classdef SoleFEM<handle
    %soleFEM class
       
    properties
        A                   % Stiffness matrix with all sole nodes
        K                   % Stiffness matrix without Dirichlet nodes
        Ks                  % Stiffness matrix of surface nodes
        Cs                  % Compliance matrix of surface nodes
        coor                % Node coordinates
        elements_surf       % Surface elements
        elements_vol        % Volume elements
        trasl               % displacement to move the sole in the correct position during update sole shape
        nTot                % number of nodes         
        nodesDirichlet      % index of Dirichlet nodes in the set of all nodes; [indx1 indy1 indz1; ... indxnDirichlet indynDirichlet indznDirichlet;]' format; 
        nodesDirichlet3     % index of Dirichlet nodes in the set of all nodes; [indx1 indy1 indz1 ... indxnDirichlet indynDirichlet indznDirichlet]' format;
        nDirichlet          % number of Dirichlet nodes
        nodesSurf           % index of Surface nodes in the set of all nodes;
        nodesFreeSurf       % index of Free surface nodes in the set of all nodes; [indx1 indy1 indz1; ... indxnFreeSurf indynFreeSurf indznFreeSurf;]' format; 
        nodesFreeSurf3      % index of Free surface nodes in the set of all nodes; [indx1 indy1 indz1 ... indxnFreeSurf indynFreeSurf indznFreeSurf]' format;
        nFreeSurf           % number of Free surface nodes
        nodesInt            % index of Internal surface nodes in the set of all nodes; [indx1 indy1 indz1; ... indxnInt indynInt indznInt;]' format; 
        nodesInt3           % index of Internal surface nodes in the set of all nodes; [indx1 indy1 indz1 ... indxnInt indynInt indznInt]' format;
        nInt                % number of Internal surface nodes
        nodesFree           % index of Free (FreeSurf + Int) nodes in the set of all nodes; [indx1 indy1 indz1; ... indxnFree indynFree indznFree;]' format; 
        nodesFree3          % index of Free (FreeSurf + Int) nodes in the set of all nodes; [indx1 indy1 indz1 ... indxnFree indynFree indznFree]' format; 
        nFree               % number of Free nodes
        m_invKii_Kis        
        dof                 % dof table
        E                   % Young's modulus
        nu                  % Poisson's modulus
        Kid                 
        Kii
        Kis
        der_A               % derivative of A with respect to coor
        der_Cs_dp           % derivative of Cc with respect to p
        mu
        lambda
        ABe                 % Preparation of Stress Von Mises computation 
        stressVM            % Stress Von Mises 
        vol                 % Sole Volume
        lx_foot             % width of the foot
        ly_foot             % length of the foot
        lz_foot             % thickness of the foot
        pos_ankle_loc       % ankle position respect to sole frame
    end     
    methods
        function obj = SoleFEM(pname,fname,l,L,e)
            [obj.elements_surf,obj.elements_vol,obj.coor] = input_mesh(pname,fname);     % reads nodes and elements from GMSH file
            obj.lx_foot = L;
            obj.ly_foot = l;
            obj.lz_foot = e;
            %obj.nEle = size(obj.elements_surf,1) + size(obj.elements_vol,1);
            nodes = (1:1:size(obj.coor,1))';
            obj.nTot = length(nodes);
            obj.trasl = [obj.lx_foot/2,0,obj.lz_foot];
            % Dirichlet nodes
            obj.nodesDirichlet = find(obj.coor(:,3)==obj.lz_foot); 
            obj.nDirichlet = length(obj.nodesDirichlet);
            % Surface nodes
            obj.nodesSurf = unique(obj.elements_surf,'first');
            % Free surface nodes
            obj.nodesFreeSurf = setdiff(obj.nodesSurf,obj.nodesDirichlet);
            obj.nFreeSurf = length(obj.nodesFreeSurf);
            % Nodes of volume - Internal nodes
            obj.nodesInt = nodes;
            obj.nodesInt([obj.nodesFreeSurf; obj.nodesDirichlet])=[];
            obj.nInt = length(obj.nodesInt);
            % Nodes free - Free surface nodes + Internal nodes
            obj.nodesFree = nodes;
            obj.nodesFree(obj.nodesDirichlet)=[];
            obj.nFree = length(obj.nodesFree);
            % Compute dof matrix
            cont = 0;
            obj.dof = zeros(size(obj.coor));
            contDir = 1;
            for i=1:obj.nTot
                if contDir <= obj.nDirichlet
                    if i~=obj.nodesDirichlet(contDir)
                        for j=1:3
                            cont = cont + 1;
                            obj.dof(i,j) = cont;
                        end
                    else
                        contDir = contDir + 1;
                    end
                else
                     for j=1:3
                        cont = cont + 1;
                        obj.dof(i,j) = cont;
                    end                   
                end
            end
            % Nodes in x;y;z;
            obj.nodesFree3(1:3:3*obj.nFree-2) = 3*obj.nodesFree-2;
            obj.nodesFree3(2:3:3*obj.nFree-1) = 3*obj.nodesFree-1;
            obj.nodesFree3(3:3:3*obj.nFree-0) = 3*obj.nodesFree-0;            
            obj.nodesFreeSurf3(1:3:3*obj.nFreeSurf-2) = 3*obj.nodesFreeSurf-2;
            obj.nodesFreeSurf3(2:3:3*obj.nFreeSurf-1) = 3*obj.nodesFreeSurf-1;
            obj.nodesFreeSurf3(3:3:3*obj.nFreeSurf-0) = 3*obj.nodesFreeSurf-0;
            obj.nodesInt3(1:3:3*obj.nInt-2) = 3*obj.nodesInt-2;
            obj.nodesInt3(2:3:3*obj.nInt-1) = 3*obj.nodesInt-1;
            obj.nodesInt3(3:3:3*obj.nInt-0) = 3*obj.nodesInt-0;
            obj.nodesDirichlet3(1:3:3*obj.nDirichlet-2) = 3*obj.nodesDirichlet-2;
            obj.nodesDirichlet3(2:3:3*obj.nDirichlet-1) = 3*obj.nodesDirichlet-1;
            obj.nodesDirichlet3(3:3:3*obj.nDirichlet-0) = 3*obj.nodesDirichlet-0;          
        end
        function [obj] = setMaterial(obj,Young,Poisson)
            obj.E = Young;
            obj.nu = Poisson;
        end
        function [obj] = stiffness(obj)
            obj.mu=obj.E/(2*(1+obj.nu));
            obj.lambda=obj.E*obj.nu/((1+obj.nu)*(1-2*obj.nu));
            obj.A = sparse(3*obj.nTot,3*obj.nTot);
            % element data file for tetraeder: 1-node/2-node/3-node/4-node
            for j = 1:size(obj.elements_vol,1)
              I = 3*obj.elements_vol(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
              obj.A(I,I) = obj.A(I,I) + stima(obj.coor(obj.elements_vol(j,:),:),obj.lambda,obj.mu); % this is the stiffness matrix
            end 
            obj.K = obj.A(obj.nodesFree3,obj.nodesFree3);
            obj.Kid = obj.A(obj.nodesInt3,obj.nodesDirichlet3);
            obj.stiffnessSurface();
        end
        function [obj] = derStiff(obj,lp,dPFree_dp)
            obj.der_A = cell(3*obj.nTot,1);
            for i = 1:(3*obj.nTot)
                obj.der_A{i}=sparse(3*obj.nTot,3*obj.nTot);
            end
            
            for j = 1:size(obj.elements_vol,1)
                vertices = obj.coor(obj.elements_vol(j,:),:);
                detJ = det([1,1,1,1;vertices']);
                if detJ > 0
                    S = 1;
                else
                    S = -1;
                end
                I = 3*obj.elements_vol(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
                obj.der_A{I(1),1}(I,I) = obj.der_A{I(1),1}(I,I) + S * der_stima_dx1(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(2),1}(I,I) = obj.der_A{I(2),1}(I,I) + S * der_stima_dy1(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(3),1}(I,I) = obj.der_A{I(3),1}(I,I) + S * der_stima_dz1(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(4),1}(I,I) = obj.der_A{I(4),1}(I,I) + S * der_stima_dx2(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(5),1}(I,I) = obj.der_A{I(5),1}(I,I) + S * der_stima_dy2(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(6),1}(I,I) = obj.der_A{I(6),1}(I,I) + S * der_stima_dz2(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(7),1}(I,I) = obj.der_A{I(7),1}(I,I) + S * der_stima_dx3(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(8),1}(I,I) = obj.der_A{I(8),1}(I,I) + S * der_stima_dy3(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(9),1}(I,I) = obj.der_A{I(9),1}(I,I) + S * der_stima_dz3(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(10),1}(I,I) = obj.der_A{I(10),1}(I,I) + S * der_stima_dx4(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(11),1}(I,I) = obj.der_A{I(11),1}(I,I) + S * der_stima_dy4(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(12),1}(I,I) = obj.der_A{I(12),1}(I,I) + S * der_stima_dz4(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
            end
            obj.der_Cs_dp = cell(lp,1);
            F = full(obj.m_invKii_Kis);
            for i = 1:lp
                der_Kss_dpi = sparse(zeros(3*obj.nFreeSurf,3*obj.nFreeSurf));
                der_Kis_dpi = sparse(zeros(3*obj.nInt,3*obj.nFreeSurf));
                der_Ksi_dpi = sparse(zeros(3*obj.nFreeSurf,3*obj.nInt));
                der_Kii_dpi = sparse(zeros(3*obj.nInt,3*obj.nInt));
                s = find(dPFree_dp(:,i) ~=0);
                for j = s'
                    der_Kss     = obj.der_A{j}(obj.nodesFreeSurf3,obj.nodesFreeSurf3);
                    der_Kss_dpi = der_Kss_dpi + der_Kss * dPFree_dp(j,i);
                    der_Ksi     = obj.der_A{j}(obj.nodesFreeSurf3,obj.nodesInt3);
                    der_Ksi_dpi = der_Ksi_dpi + der_Ksi * dPFree_dp(j,i);
                    der_Kis     = obj.der_A{j}(obj.nodesInt3,obj.nodesFreeSurf3);
                    der_Kis_dpi = der_Kis_dpi + der_Kis * dPFree_dp(j,i);                    
                    der_Kii     = obj.der_A{j}(obj.nodesInt3,obj.nodesInt3);
                    der_Kii_dpi = der_Kii_dpi + der_Kii * dPFree_dp(j,i);
                end
                tmp = der_Ksi_dpi * obj.m_invKii_Kis;
                O = der_Kss_dpi - tmp - tmp';
                O1 = F' * full(der_Kii_dpi) * F;
                O = full(O) + O1;
                obj.der_Cs_dp{i} = -obj.Cs * O * obj.Cs;
            end
            
        end     
        
        function [obj] = stiffnessSurface(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the surface Stiffness Ks [Fi;Fs]=[Kii Kis;Ksi Kss]*[Xi;Xs]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            KInternalNodes = reshape(obj.dof(obj.nodesInt,:)',[1,3*obj.nInt]);
            KNodesFreeSurf = reshape(obj.dof(obj.nodesFreeSurf,:)',[1,3*obj.nFreeSurf]);
            
            obj.Kii = obj.K(KInternalNodes,KInternalNodes);
            obj.Kis = obj.K(KInternalNodes,KNodesFreeSurf);
            Ksi = obj.K(KNodesFreeSurf,KInternalNodes);
            Kss = obj.K(KNodesFreeSurf,KNodesFreeSurf);
            
            obj.m_invKii_Kis = obj.Kii\obj.Kis;           
            obj.Cs = full(inv(Kss - Ksi * obj.m_invKii_Kis));
        end
        
        function [dP,obj] = getdP(obj,dPFreeSurf,dPDir)
            dPSurf3 = reshape(dPFreeSurf([1:obj.nFreeSurf],:)',[1, 3*obj.nFreeSurf]);
            dP = zeros(obj.nTot,3);
            dP(obj.nodesFreeSurf,:) = dPFreeSurf;
            dP(obj.nodesDirichlet,:) = dPDir;
            if ~isempty(obj.m_invKii_Kis)
                dPInt = -obj.m_invKii_Kis*dPSurf3';
                dP(obj.nodesInt,:) = reshape(dPInt,3,length(obj.nodesInt))';
            end
        end        
        function [dP,obj] = getdPTot(obj,dPFreeSurf,dPDir)
            dPFreeSurf3 = reshape(dPFreeSurf([1:obj.nFreeSurf],:)',[1, 3*obj.nFreeSurf]);
            dPDir3 = reshape(dPDir([1:obj.nDirichlet],:)',[1, 3*obj.nDirichlet]);
            dP = zeros(obj.nTot,3);
            dP(obj.nodesFreeSurf,:) = dPFreeSurf;
            dP(obj.nodesDirichlet,:) = dPDir;
            dPInt3 = -inv(obj.Kii)*(obj.Kid*dPDir3' + obj.Kis*dPFreeSurf3');
            dP(obj.nodesInt,:) = reshape(dPInt3,3,length(obj.nodesInt))';
        end
        function obj = prep_stressVonMises(obj)
            % Preparation of VonMises' Stress computations
            obj.ABe = cell(size(obj.elements_vol,1),1);
            for j = 1:size(obj.elements_vol,1) % computes nodal stresses
                vertices = obj.coor(obj.elements_vol(j,:),:);
                PhiGrad = [1,1,1,1;vertices']\[zeros(1,3);eye(3)];
                Be = zeros(6,12);
                Be([1,4,5],1:3:10) = PhiGrad';
                Be([4,2,6],2:3:11) = PhiGrad';
                Be([5,6,3],3:3:12) = PhiGrad';
                Ce(1:3,1:3) = obj.lambda*ones(3,3)+2*obj.mu*eye(3);
                Ce(4:6,4:6) = obj.mu*eye(3);
                obj.ABe{j} = Ce*Be;    
            end
        end
        function obj = stressVonMises(obj,dPabstot)
            % Computation of the VonMises' Stress
            displEle = zeros(obj.nTot,3);
            dPabstot(obj.nodesDirichlet3) = [];  
            ic = find(obj.dof>0); % finds active dofs

            displEle(ic) = dPabstot(obj.dof(ic));

            stress = zeros(obj.nTot,6);                             % initializes "stress"
            stressG2 = zeros(4,6);
            counter = zeros(obj.nTot,1);                            % initializes "counter"
            for e = 1:size(obj.elements_vol,1)                      % computes nodal stresses
                nodes = obj.elements_vol(e,:);
                stressG2(1,:) = [displEle(nodes(1),:) displEle(nodes(2),:) displEle(nodes(3),:) displEle(nodes(4),:)]*obj.ABe{e}';
                stressG2(2,:) = stressG2(1,:);
                stressG2(3:4,:) = stressG2(1:2,:);
                stress(nodes,:) = stress(nodes,:) + stressG2;   % adds stresses to triangle nodes
                counter(nodes) = counter(nodes) + 1;        
            end
            for icomp=1:6
                stress(:,icomp)=stress(:,icomp)./counter;  % naive average of stresses
            end

            obj.stressVM=sqrt((3/2)*(((2*stress(:,1)-stress(:,2)-stress(:,3))/3).^2 ...
                    +((2*stress(:,2)-stress(:,1)-stress(:,3))/3).^2 ...
                    +((2*stress(:,3)-stress(:,2)-stress(:,1))/3).^2 ...
                    +(stress(:,4)).^2+(stress(:,5)).^2+(stress(:,6)).^2));
        end        
        function obj = computeVol(obj)
            detJ = zeros(size(obj.elements_vol,1),1);
            obj.vol = 0;
            for j = 1:size(obj.elements_vol,1)
                vertices = obj.coor(obj.elements_vol(j,:),:);
                detJ(j) = det([1,1,1,1;vertices']);
                obj.vol = obj.vol + abs(detJ(j)/6);
            end
        end
    end
end