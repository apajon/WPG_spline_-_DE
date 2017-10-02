classdef Simu<handle
    % Simu class
       
    properties
        displ       % Ol in the paper, displacement of the sole
        angleact    % Upsilon in the paper, angles defing the sole rotation
        theta       % angle around x axis
        phi         % angle around y axis
        psi         % angle around z axis
        R           % Sole rotation matrix
        Rini        % Sole rotation matrix for the first time step
        Rtheta      % Sole rotation matrix around x axis
        Rphi        % Sole rotation matrix around y axis
        Rpsi        % Sole rotation matrix around z axis
        RpsiFirst   % Sole rotation matrix around z axis for the first time step
        displ_first % displacement of the sole for the first time step
        psi_first   % roll angle for the first step
        P           % node positions
        PfreeSurf   % surface free node position; [x1 y1 z1; ... xm ym zm;]' format   
        PfreeSurf3  % surface free node position; [x1 y1 z1 ... xm ym zm]' format 
        PfreeCont   % contact free node position; [x1 y1 z1; ... xnc ync nc;]' format 
        PfreeCont3  % contact free node position; [x1 y1 z1 ... xm ym zm]' format      
        PSurfOld    % surface node position at previous time step
        PContOld    % contact node position at previous time step
        PSurf3      % Surface Position Vector; [x1 y1 z1 ... xm ym zm]' format 
        Pc3         % Vector of Contact Node Positions; [x1 y1 z1 ... xm ym zm]' format
        Pc          % Vector of Contact Node Positions; [x1 y1 z1; ... xm ym zm;]' format
        no_conv     % Time step when the algorithm it does not converge
        FSurf       % Surface Forces Vector; [x1 y1 z1; ...; xm ym zm;]' format
        FSurf3      % Surface Forces Vector; [x1 y1 z1 ... xm ym zm]' format
        Fc          % Contact Forces Vector; [x1 y1 z1; ... x_nc y_nc z_nc;]' format
        Fc3         % Contact Forces Vector; [x1 y1 z1 ... x_nc y_nc z_nc]' format
        Ftot        % Vector containing the total force in the 3 directions
        Z           % Vector containing the zmp position in the 2 directions [Zx;Zy]
        Fdes        % Fdes of the actual time step        
        Zdes        % Zdes of the actual time step        
        FtotZMP     % Vector Containt force tot and zmp position [Ftot;Z]        
        ind_cont    % Index of contact nodes in the surface node set
        ind_cont3   % Index of contact nodes in the surface node set; [x1 y1 z1 ... x_nc y_nc z_nc]' format
        ind_slip    % Index of contact nodes in slip condition in the surface node set
        nc          % Number of contact nodes
        ns          % Number of contact nodes in slip condition
        dPabstot    % displacement of every node deriving by the contact with the ground
        Kcart       % cartesian stiffness matrix
    end     
    methods
        function obj=Simu()

        end
        function obj=updateContact(obj)
            % Update the contact properties
            obj.nc = length(obj.ind_cont);
            obj.ind_cont3 = sort(([3*obj.ind_cont-2; 3*obj.ind_cont-1; 3*obj.ind_cont]));
            obj.ind_cont3 = reshape(obj.ind_cont3,3*obj.nc,1);
            obj.ns = length(obj.ind_slip);
            obj.Fc = reshape(obj.Fc3,3,obj.nc);
            obj.Pc3 = obj.PSurf3(obj.ind_cont3);
            obj.Pc = reshape(obj.Pc3,3,obj.nc);
            obj.PfreeCont = obj.PfreeSurf(:,obj.ind_cont);
            obj.PfreeCont3 = reshape(obj.PfreeCont,3*obj.nc,1); 
            % Update Rotation matrix
            obj.theta = obj.angleact(1);
            obj.phi = obj.angleact(2);
            obj.psi = obj.angleact(3);
            obj.Rtheta = [1 0 0; 0 cos(obj.theta) -sin(obj.theta); 0 sin(obj.theta) cos(obj.theta)];
            obj.Rphi = [cos(obj.phi) 0 sin(obj.phi); 0 1 0; -sin(obj.phi) 0 cos(obj.phi)];
            obj.Rpsi = [cos(obj.psi) -sin(obj.psi) 0; sin(obj.psi) cos(obj.psi) 0; 0 0 1];
            obj.R = obj.Rtheta*obj.Rphi*obj.Rpsi;
            % Update RotationIni matrix
            obj.RpsiFirst = [cos(obj.psi_first) -sin(obj.psi_first) 0; sin(obj.psi_first) cos(obj.psi_first) 0; 0 0 1];
            obj.Rini = obj.Rtheta*obj.Rphi*obj.RpsiFirst;
            obj.PContOld = obj.PSurfOld(:,obj.ind_cont);
        end
        function obj = updatePosition(obj,sole)
            %%% Find the displacement of internal nodes and position of all nodes
            obj.FSurf3 = zeros(3*sole.nFreeSurf,1);
            obj.FSurf3(obj.ind_cont3) = obj.Fc3;
            obj.FSurf = reshape(obj.FSurf3,3,sole.nFreeSurf)';
            Fcloc = obj.R'*obj.FSurf';
            Fcloc3 = reshape(Fcloc,3*sole.nFreeSurf,1);
            dPlocSurf3 = sole.Cs * Fcloc3;
            % Displacement of all nodes
            dPInt3 = -(sole.m_invKii_Kis) * dPlocSurf3;
            dPloc3 = zeros(3*sole.nTot,1);
            dPloc3(sole.nodesFreeSurf3) = dPlocSurf3;
            dPloc3(sole.nodesInt3) = dPInt3;
            dPloc = reshape(dPloc3,3,sole.nTot)';
            obj.dPabstot = obj.R*dPloc';
            displacement = repmat(obj.displ',sole.nTot, 1);
            PRot = obj.R*sole.coor';
            obj.P = PRot' + obj.dPabstot' + displacement;
        end             
    end
end