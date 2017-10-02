classdef GradientShape<handle
    % GradientShape class. This class is explained in the technical report
    % at the following link
       
    properties
        A                   % see the journal paper 
        inv_A               % Inverse of the matrix A
        B                   % see the journal paper
        G                   % see the journal paper
        Bw                  % see the journal paper
        Xi                  % see the journal paper
        E                   % see the journal paper
        N                   % see the journal paper
        M                   % see the journal paper
        inv_M               % Inverse of the matrix M
        Phi                 % see the journal paper
        Phi_wpg             % see the journal paper
        Theta               % see the journal paper
        Y2                  % see the journal paper
        Y2_wpg              % see the journal paper
        bi                  % see the journal paper
        bi_wpg              % see the journal paper
        Cc                  % compliance matrix of contact nodes
        Ws                  % RCsR'
        Wc                  % RCcR'
        dPFrees_dp          % derivative of PfreeSurf with respect to polynomial coefficients
        dPFreec_dp          % derivative of PfreeCont with respect to polynomial coefficients
        Fc_slip             % Forces for contact nodes in slip conditions
        Pc_slip             % Forces for contact nodes in slip conditions
        PContOld_slip       % slip contact node position at previous time step
        Islip               % index of slip nodes in the set of surface nodes
        lp                  % number of polynomial coefficients
        n_wpg               % number of wpg parameters
        deltaSlip_t         % contact displacement vector for the contact in slip conditions
        Fc_hat              % skew matrix of force contact vector
        dPcprevious3_dp     % gradient of contact node positions at the previous time step with respect to p
                            % format [dPc3_dp1;dPc3_dp2;...dPc3_dplp]
        dPs3_dp             % gradient of surface node positions with respect to p
                            % format [dPs3_dp1;dPs3_dp2;...dPs3_dplp]
        dPc3_dp             % gradient of contact node positions with respect to p
                            % format [dPc3_dp1;dPc3_dp2;...dPc3_dplp]        
        dPs3_dwpg           % gradient of surface node positions with respect to wpg parameters
                            % format [dPs3_dwpg1;dPs3_dwpg2;...dPs3_dwpgn_wpg]
        dPc3_dwpg           % gradient of contact node positions with respect to wpg parameters
                            % format [dPc3_dwpg1;dPc3_dwpg2;...dPc3_dwpgn_wpg]
        dPcprevious3_dwpg   % gradient of contact node positions at the previous time step with respect to wpg parameters
                            % format [dPc3_dwpg1;dPc3_dwpg2;...dPc3_dwpgn_wpg]                            
        der_HCcHt_dp        % derivative of HCcH' with respect to p
        der_Cc_dp           % derivative of Cc with respect to p
        t                   % matrix to select just tangential components
        t1                  % vector to select just the first tangential components
        t2                  % vector to select just the second tangential components
        dFc3_dp             % gradient of surface node forces with respect to p
        ddeltat2_dp         % gradient of deltat with respect to p (just for contact in slip condition)
        dOl_dp              % gradient of Ol(displ) with respect to p
        dY_dp               % gradient of Upsilon(vector of Angles) with respect to p
        dFc3_dwpg           % gradient of surface node forces with respect to wpg parameters
        ddeltat2_dwpg       % gradient of deltat with respect to wpg parameters (just for contact in slip condition)
        dOl_dwpg            % gradient of Ol(displ) with respect to wpg parameters
        dY_dwpg             % gradient of Upsilon(vector of Angles) with respect to wpg parameters          
        RdPFrees_dp         % R * dPFrees_dp
        RdPFreec_dp         % R * dPFreec_dp
        dmu                 % [dmu_dtheta dmu_dphi dmu_dpsi], dmu_dtheta = [dmu_dtheta_1...dmu_dtheta_m], m are the number of surface nodes
                            % with dmu_dtheta_j = dR_dtheta*P_j + sum_j^m(dR_dtheta * Cs * R' * Fj + R * Cs * dR_dtheta' * Fj)
        dR_dtheta           % derivative of rotation matrix with respect to theta
        dR_dphi             % derivative of rotation matrix with respect to phi
        dR_dpsi             % derivative of rotation matrix with respect to psi
        dWs_dFs             % Ws(:,ind_cont3)*dFc3_dp;
        dWs_dFs_wpg         % Ws(:,ind_cont3)*dFc3_dwpg;
        dWc_dFc             % Wc*dFc3_dp;
        dWc_dFc_wpg         % Wc*dFc3_dwpg;
        invABwXi            % inv_A * Bw * Xi;
        G_inv_A             % G * inv_A;
        G_inv_ABw           % G * inv_A * Bw;    
        u_max               % eigen vector of the max eigenvalue(Kc^tr)
        u_min               % eigen vector of the min eigenvalue(Kc^rot)
        derG_pi             % derivative of G with respect to pi
        derG_wpgi           % derivative of G with respect to wpgi
        derA_pi             % derivative of A with respect to pi
        derA_wpgi           % derivative of A with respect to wpgi
        derW_pi             % derivative of W with respect to pi
        derW_wpgi           % derivative of W with respect to wpgi
        dR_dpi              % derivative of R with respect to pi
        dR_dwpgi            % derivative of R with respect to wpgi
        derBw_pi            % derivative of Bw with respect to pi
        derBw_wpgi          % derivative of Bw with respect to wpgi
        derXi_pi            % derivative of Xi with respect to pi
        derXi_wpgi          % derivative of Xi with respect to wpgi
        derKc_pi            % derivative of Kc with respect to pi
        derKc_wpgi          % derivative of Kc with respect to wpgi
        der_costS_tr_dp     % derivative of costShape^tr with respect to p
        der_costS_rot_dp    % derivative of costShape^rot with respect to p
        der_costS_tr_dwpg   % derivative of costShape^tr with respect to wpg
        der_costS_rot_dwpg  % derivative of costShape^rot with respect to wpg   
    end     
    methods
        function obj = GradientShape()

        end
        function obj = updContVar(obj,sole,simu,friction,contAngle)
            %%% With this function, we update the contact variables to
            %%% compute gradient
            % We compute this part here and not in simu class because we
            % just need for gradient computation
            obj.Cc = sole.Cs(simu.ind_cont3,simu.ind_cont3);
            % RCsR' with optimized code
            nf = sole.nFreeSurf;
            tmp = reshape(simu.R*reshape(sole.Cs,3,3*nf*nf),3*nf,3*nf);
            tmp2 = tmp';
            obj.Ws = reshape(simu.R*reshape(tmp2,3,3*nf*nf),3*nf,3*nf);
            %%%
            obj.Wc = obj.Ws(simu.ind_cont3,simu.ind_cont3);
            % Update Slip variables
            obj.Islip = zeros(simu.ns,1);
            for i=1:simu.ns
                obj.Islip(i) = find(simu.ind_cont==simu.ind_slip(i));
            end
            obj.Fc_slip = simu.Fc(:,obj.Islip);
            obj.Pc_slip = simu.Pc(:,obj.Islip);
            obj.PContOld_slip = simu.PContOld(:,obj.Islip);
            obj.t = [1 0 0; 0 1 0];
            obj.t1 = [1 0 0];
            obj.t2 = [0 1 0];
            obj.deltaSlip_t = zeros(2,simu.ns);
            for i=1:simu.ns
                if contAngle==1
                    delta = obj.Pc_slip(:,i) - [obj.t * simu.displ_first + obj.t * simu.Rini * obj.PContOld_slip(:,i);0];
                else
                    delta = obj.Pc_slip(:,i) - [obj.t * obj.PContOld_slip(:,i);0];
                end   
                obj.deltaSlip_t(:,i) = delta(1:2,1);
            end
            obj.dRot(simu);
            obj.compute_G(simu);
            obj.compute_A(simu,friction);
            obj.compute_B(simu,contAngle);
            obj.inv_A = inv(obj.A);
            obj.compute_Bw(simu);
            obj.compute_Xi(simu);
            obj.compute_E(simu);  
            if contAngle==1
                obj.compute_N(simu);
            end
            obj.compute_M(contAngle);
            obj.inv_M = inv(obj.M);
            Kc_squared_tr = simu.Kcart(1:3,1:3)'*simu.Kcart(1:3,1:3);
            [u,DE] = eig(Kc_squared_tr);
            [~,I_lamba_max] = max(diag(DE));
            obj.u_max = u(:,I_lamba_max);
            Kc_squared_rot = simu.Kcart(4:6,4:6)'*simu.Kcart(4:6,4:6);
            [u,DE] = eig(Kc_squared_rot);
            [~,I_lamba_min] = min(diag(DE));
            obj.u_min = u(:,I_lamba_min); 
            obj.invABwXi = obj.inv_A * obj.Bw * obj.Xi;
            obj.G_inv_A = obj.G * obj.inv_A;
            obj.G_inv_ABw = obj.G * obj.inv_A * obj.Bw;
        end
        function obj = derShape_pol(obj,sole,simu,dPFree_dp,friction,contAngle)
            % Compute gradient of cost with respect to the polynomial
            % coefficients
            obj.lp = size(dPFree_dp,2);
            obj.dPFrees_dp = dPFree_dp(sole.nodesFreeSurf3,:);
            obj.dPFreec_dp = obj.dPFrees_dp(simu.ind_cont3,:);       
            obj.der_Cc_dp = cell(obj.lp,1);
            for i = 1:obj.lp
                obj.der_Cc_dp{i} = sole.der_Cs_dp{i}(simu.ind_cont3,simu.ind_cont3);
            end
            if contAngle==1
                obj.dPcprevious3_dp = zeros(3*simu.nc,obj.lp);
            else
                obj.dPcprevious3_dp = obj.dPs3_dp(simu.ind_cont3,:);
                obj.Phi = obj.dPcprevious3_dp;
                obj.Phi(3:3:end,:) = 0;
            end
            
            obj.compute_Theta(simu,contAngle);
            obj.compute_Y2(simu,contAngle);
            obj = compute_bi(obj,simu);
            dXi_dp = obj.inv_M * obj.bi;
            obj.dFc3_dp = dXi_dp(1:3*simu.nc,:);
            obj.ddeltat2_dp = dXi_dp((3*simu.nc+1):((3*simu.nc)+2*simu.ns),:);
            obj.dOl_dp = dXi_dp(((3*simu.nc)+2*simu.ns)+1:((3*simu.nc)+2*simu.ns)+3,:);
            obj.dY_dp = dXi_dp(((3*simu.nc)+2*simu.ns)+4:end,:);
            obj.derPs3_dp(sole,simu); 
            obj.RdPFreec_dp = obj.RdPFrees_dp(simu.ind_cont3);
            obj.dWc_dFc = obj.dWs_dFs(simu.ind_cont3,:);
            obj.dPc3_dp =  obj.dPs3_dp(simu.ind_cont3,:);
            obj.der_costS_tr_dp = zeros(obj.lp,1);
            obj.der_costS_rot_dp = zeros(obj.lp,1);
            for i=1:obj.lp
                obj.compute_derG_pi(simu,i);
                obj.compute_derA_pi(simu,friction,i);
                obj.compute_derBw_pi(simu,i);
                obj.compute_derXi_pi(i);
                obj.derKc_pi = obj.derG_pi * obj.invABwXi - obj.G_inv_A * obj.derA_pi * obj.invABwXi +...
                               obj.G_inv_A * obj.derBw_pi * obj.Xi + obj.G_inv_ABw * obj.derXi_pi;
                obj.compute_dercost_dp(simu,i);
            end
        end
        function obj = derShape_dwpg(obj,sole,simu,friction,gradZdesx,gradZdesy,gradFdesx,gradFdesy,gradFdesz,contAngle)
            % Compute gradient of the shape cost function with respect to the WPG parameters
            obj.n_wpg = length(gradFdesx);         
            
            if contAngle==1
                obj.Phi_wpg = zeros(3*simu.nc,obj.n_wpg);
            else           
                obj.dPcprevious3_dwpg = obj.dPs3_dwpg(simu.ind_cont3,:);
                obj.Phi_wpg = obj.dPcprevious3_dwpg;
            end
            %%% gradient of the shape cost function with respect to wpg parameters
            obj.compute_bi_wpg(simu,gradFdesx,gradFdesy,gradFdesz,gradZdesx,gradZdesy,contAngle);
            dXi_dwpg = obj.inv_M * obj.bi_wpg;
            obj.dFc3_dwpg = dXi_dwpg(1:3*simu.nc,:);
            obj.ddeltat2_dwpg = dXi_dwpg((3*simu.nc+1):((3*simu.nc)+2*simu.ns),:);
            obj.dOl_dwpg = dXi_dwpg(((3*simu.nc)+2*simu.ns)+1:((3*simu.nc)+2*simu.ns)+3,:);
            obj.dY_dwpg = dXi_dwpg(((3*simu.nc)+2*simu.ns)+4:end,:);
            obj.derPs3_dwpg(sole,simu);
            obj.dWc_dFc_wpg = obj.dWs_dFs_wpg(simu.ind_cont3,:);
            obj.dPc3_dwpg = obj.dPs3_dwpg(simu.ind_cont3,:);
            obj.der_costS_tr_dwpg = zeros(obj.n_wpg,1);
            obj.der_costS_rot_dwpg = zeros(obj.n_wpg,1);
            for i=1:obj.n_wpg
                obj.compute_derG_wpgi(simu,i,gradZdesx(i),gradZdesy(i));
                obj.compute_derA_wpgi(simu,friction,i);
                obj.compute_derBw_wpgi(simu,i);
                obj.compute_derXi_wpgi(i,gradZdesx(i),gradZdesy(i));
                obj.derKc_wpgi = obj.derG_wpgi * obj.invABwXi - obj.G_inv_A * obj.derA_wpgi * obj.invABwXi +...
                                  obj.G_inv_A * obj.derBw_wpgi * obj.Xi + obj.G_inv_ABw * obj.derXi_wpgi;
                obj.compute_dercost_dwpg(simu,i);
            end            
            
        end
        function obj = compute_A(obj,simu,friction)
            A11 = -obj.Wc;
            A12 = zeros(3*simu.nc,2*simu.ns);
            A21 = zeros(2*simu.ns,3*simu.nc);
            A22 = zeros(2*simu.ns,2*simu.ns);
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    A21((2*j)-1:2*j,(3*obj.Islip(j))-2:3*obj.Islip(j)) = [norm(obj.deltaSlip_t(:,j))*eye(2) friction*obj.deltaSlip_t(:,j)];
                    A12((3*obj.Islip(j))-2:3*obj.Islip(j),(2*j)-1:2*j) = [1 0; 0 1; 0 0];
                    A22((2*j)-1:2*j,(2*j)-1:2*j) =  obj.Fc_slip(1:2,j)*(obj.deltaSlip_t(:,j)'./norm(obj.deltaSlip_t(:,j)))+friction*obj.Fc_slip(3,j)*eye(2);
                end
            end
            obj.A = [A11 A12;A21 A22];           
        end
        function obj = compute_B(obj,simu,contAngle)
            if contAngle==1            
                dRinidtheta = [0 0 0; 0 -sin(simu.theta) -cos(simu.theta); 0 cos(simu.theta) -sin(simu.theta)]*simu.Rphi*simu.RpsiFirst;        
                dRinidphi = simu.Rtheta*[-sin(simu.phi) 0 cos(simu.phi); 0 0 0; -cos(simu.phi) 0 -sin(simu.phi)]*simu.RpsiFirst;
            end
            for i=1:simu.nc
                if contAngle==1
                    btheta = obj.dR_dtheta*simu.PfreeCont3((3*i-2):3*i);
                    bphi = obj.dR_dphi*simu.PfreeCont3((3*i-2):3*i);
                    bpsi = obj.dR_dpsi*simu.PfreeCont3((3*i-2):3*i);
                    btheta(1) = btheta(1) - dRinidtheta(1,:)*simu.PfreeCont3((3*i-2):3*i);
                    btheta(2) = btheta(2) - dRinidtheta(2,:)*simu.PfreeCont3((3*i-2):3*i);
                    bphi(1) = bphi(1) - dRinidphi(1,:)*simu.PfreeCont3((3*i-2):3*i);
                    bphi(2) = bphi(2) - dRinidphi(2,:)*simu.PfreeCont3((3*i-2):3*i);
                else
                    btheta = obj.dR_dtheta*simu.PfreeCont3((3*i-2):3*i);
                    bphi = obj.dR_dphi*simu.PfreeCont3((3*i-2):3*i);
                    bpsi = obj.dR_dpsi*simu.PfreeCont3((3*i-2):3*i);           
                end
                for j=1:simu.nc
                    F = simu.Fc3(3*j-2:3*j);
                    C = obj.Cc(3*i-2:3*i,3*j-2:3*j);
                    RtF = simu.R'*F;
                    CRtF = C*RtF;
                    a = obj.dR_dtheta * CRtF + simu.R * (C * (obj.dR_dtheta' * F));
                    b = obj.dR_dphi * CRtF + simu.R * (C * (obj.dR_dphi' * F));
                    c = obj.dR_dpsi  * CRtF + simu.R * (C * (obj.dR_dpsi' * F));
                    btheta = btheta + a;
                    bphi = bphi + b; 
                    bpsi = bpsi + c;
                end
                B12((3*i-2):3*i,:) = [btheta bphi bpsi];
            end
            B11 = repmat(eye(3),simu.nc,1);
            obj.B = [B11 B12;zeros(2*simu.ns,3) zeros(2*simu.ns,3)];
        end        
        function obj = compute_G(obj,simu)
            G11 = repmat(eye(3),1,simu.nc);
            G21 = zeros(3,3*simu.nc);
            for j=1:simu.nc
                G21(:,(3*j)-2:3*j) = [0 0 simu.Pc(2,j)-simu.Zdes(2);0 0 -(simu.Pc(1,j)-simu.Zdes(1));-simu.Pc(2,j)+simu.Zdes(2) simu.Pc(1,j)-simu.Zdes(1) 0];
            end
            G22 = zeros(3,2*simu.ns);
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    G22(:,(2*j)-1:2*j) = [0 obj.Fc_slip(3,j);-obj.Fc_slip(3,j) 0;obj.Fc_slip(2,j) -obj.Fc_slip(1,j)];
                end
            end
            obj.G = [G11 zeros(3,2*simu.ns);G21 G22];        
        end
        function obj = compute_Bw(obj,simu)
            delta_mc3 = obj.Wc * simu.Fc3;

            obj.Fc_hat = multiSkew(simu.Fc3);
            delta_mc_hat = multiSkew(delta_mc3);
            
            Bw12 = obj.Wc * obj.Fc_hat - delta_mc_hat;
            for j=1:simu.nc
                Bw12(3*j-2:3*j,:) = Bw12(3*j-2:3*j,:)-skew(simu.R*simu.PfreeCont(:,j));
            end
            Bw11 = repmat(eye(3),simu.nc,1);
            obj.Bw = [Bw11 Bw12;zeros(2*simu.ns,3) zeros(2*simu.ns,3)];
        end
        function obj = compute_Xi(obj,simu)
            obj.Xi = [1, 0, 0, 0,                           simu.displ(3),                  -simu.displ(2)+simu.Zdes(2);
                      0, 1, 0, -simu.displ(3),              0,                              simu.displ(1)-simu.Zdes(1);
                      0, 0, 1, simu.displ(2)-simu.Zdes(2),  -simu.displ(1)+simu.Zdes(1),    0;
                      0, 0, 0, 1,                           0,                              0;
                      0, 0, 0, 0,                           1,                              0;
                      0, 0, 0, 0,                           0,                              1];
        end        
        function obj = compute_E(obj,simu)
            E11 = repmat(eye(3),1,simu.nc);
            E21 = zeros(3,3*simu.nc);
            for j=1:simu.nc
                E21(:,(3*j)-2:3*j) = [0 0 -(simu.Pc(1,j)-simu.Zdes(1));0 0 (simu.Pc(2,j)-simu.Zdes(2));-(simu.Pc(2,j)-simu.Zdes(2)) (simu.Pc(1,j)-simu.Zdes(1)) 0];
            end
            E22 = zeros(3,2*simu.ns);
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    E22(:,(2*j)-1:2*j) = [-obj.Fc_slip(3,j) 0;0 obj.Fc_slip(3,j);obj.Fc_slip(2,j) -obj.Fc_slip(1,j)];
                end
            end
            obj.E = [E11 zeros(3,2*simu.ns);E21 E22];
        end
        function obj = compute_M(obj,contAngle)
            if contAngle==1
                obj.M = [obj.A -obj.B;obj.E obj.N];
            else
                obj.M = [obj.A -obj.B;obj.E zeros(6,6)];
            end
        end        
        function obj = compute_Theta(obj,simu,contAngle)
            obj.Theta = zeros(3*simu.nc,obj.lp);
            obj.der_HCcHt_dp = cell(obj.lp,1);
            if contAngle==1
                Rini1 = [obj.t * simu.Rini;0 0 0];
                for i=1:obj.lp
                    tmp = reshape(simu.R*reshape(obj.der_Cc_dp{i},3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
                    tmp2 = tmp';
                    obj.der_HCcHt_dp{i} = reshape(simu.R*reshape(tmp2,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc); %obj.der_HCcHt_dp{i} = Hc * der_Cc_dp{i} * Hc';                    
                    der_HCcHF = obj.der_HCcHt_dp{i} * simu.Fc3;            
                    obj.Theta(:,i) = der_HCcHF + reshape((simu.R - Rini1) * reshape(obj.dPFreec_dp(:,i),3,simu.nc),3*simu.nc,1); %obj.Theta(:,i) = Hc * der_Cc_dp{i} * Hc' * Fc_mc3 + Hc * dPFreec_dp(:,i) - Hc_ini * dPFreec_dp(:,i);           
                end
            else
                for i=1:obj.lp
                    tmp = reshape(simu.R*reshape(obj.der_Cc_dp{i},3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
                    tmp2 = tmp';
                    obj.der_HCcHt_dp{i} = reshape(simu.R*reshape(tmp2,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc); %obj.der_HCcHt_dp{i} = Hc * der_Cc_dp{i} * Hc';
                    der_HCcHF = obj.der_HCcHt_dp{i} * simu.Fc3;
                    obj.Theta(:,i) = der_HCcHF + reshape(simu.R* reshape(obj.dPFreec_dp(:,i),3,simu.nc),3*simu.nc,1) - obj.Phi(:,i);%obj.Theta(:,i) = der_HCcHF + Hc * dPFreec_dp(:,i) - dPc_dp(:,i);
                end
            end
        end
        function obj = compute_N(obj,simu)
            obj.N = zeros(6,6);          
            dRinidtheta = [0 0 0; 0 -sin(simu.theta) -cos(simu.theta); 0 cos(simu.theta) -sin(simu.theta)]*simu.Rphi*simu.RpsiFirst;        
            dRinidphi = simu.Rtheta*[-sin(simu.phi) 0 cos(simu.phi); 0 0 0; -cos(simu.phi) 0 -sin(simu.phi)]*simu.RpsiFirst;    
            PfreeFt1 = zeros(3,1);
            PfreeFt2 = zeros(3,1);
            PfreeFn = zeros(3,1);
            for j=1:simu.nc
                PfreeFt1 = PfreeFt1 + simu.PfreeCont3((3*j-2):3*j) * simu.Fc3(3*j-2);
                PfreeFt2 = PfreeFt2 + simu.PfreeCont3((3*j-2):3*j) * simu.Fc3(3*j-1);
                PfreeFn = PfreeFn + simu.PfreeCont3((3*j-2):3*j) * simu.Fc3(3*j);     
            end
            obj.N(4,4) = -obj.t1 * dRinidtheta * PfreeFn;
            obj.N(4,5) = -obj.t1 * dRinidphi * PfreeFn;
            obj.N(5,4) = obj.t2 * dRinidtheta * PfreeFn;
            obj.N(5,5) = obj.t2 * dRinidphi * PfreeFn;
            obj.N(6,4) = obj.t1 * dRinidtheta * PfreeFt2 - obj.t2 * dRinidtheta * PfreeFt1;
            obj.N(6,5) = obj.t1 * dRinidphi * PfreeFt2 - obj.t2 * dRinidphi * PfreeFt1;
        end
      
        function obj = compute_bi(obj,simu)
            obj.bi = [obj.Theta;zeros(2*simu.ns,obj.lp);zeros(3,obj.lp);obj.Y2];
        end
        function obj = compute_bi_wpg(obj,simu,gradFdesx,gradFdesy,gradFdesz,gradZdesx,gradZdesy,contAngle)
            obj.compute_Y2_wpg(simu,contAngle,gradZdesx,gradZdesy);
            obj.bi_wpg = [-obj.Phi_wpg;zeros(2*simu.ns,obj.n_wpg);[gradFdesx;gradFdesy;gradFdesz];obj.Y2_wpg];
        end
        function obj = compute_Y2(obj,simu,contAngle)  
            obj.Y2 = zeros(3,obj.lp);
            if contAngle==1
                for i=1:obj.lp
                    At1 = zeros(3,1);
                    At2 = zeros(3,1);
                    An = zeros(3,1);
                    for j=1:simu.nc
                        At1 = At1 + simu.Rini * obj.dPFreec_dp(3*j-2:3*j,i) * simu.Fc3(3*j-2);
                        At2 = At2 + simu.Rini * obj.dPFreec_dp(3*j-2:3*j,i) * simu.Fc3(3*j-1);
                        An = An + simu.Rini * obj.dPFreec_dp(3*j-2:3*j,i) * simu.Fc3(3*j);
                    end
                    obj.Y2(:,i) = [obj.t1*An; -obj.t2*An; -obj.t1*At2+obj.t2*At1];
                end
            else 
                for i=1:obj.lp
                    %for j=1:length(ind_cont)
                    %    At1 = At1 + dPcq_1_dp(3*j-2:3*j,i) * Fc_mc3(3*j-2);
                    %    At2 = At2 + dPcq_1_dp(3*j-2:3*j,i) * Fc_mc3(3*j-1);
                    %    An = An + dPcq_1_dp(3*j-2:3*j,i) * Fc_mc3(3*j);
                    %end
                    At1 = reshape(obj.dPcprevious3_dp(:,i),3,simu.nc)*simu.Fc3(1:3:3*simu.nc);
                    At2 = reshape(obj.dPcprevious3_dp(:,i),3,simu.nc)*simu.Fc3(2:3:3*simu.nc);
                    An = reshape(obj.dPcprevious3_dp(:,i),3,simu.nc)*simu.Fc3(3:3:3*simu.nc);
                    obj.Y2(:,i) = [obj.t1*An; -obj.t2*An; -obj.t1*At2+obj.t2*At1];
                end
            end
        end
        function obj = compute_Y2_wpg(obj,simu,contAngle,gradZdesx,gradZdesy)  
            if contAngle==1
                Ftotn = sum(simu.Fc3(3:3:end));
                Ftott1 = sum(simu.Fc3(1:3:end));
                Ftott2 = sum(simu.Fc3(2:3:end));
                obj.Y2_wpg = zeros(3,obj.n_wpg);
                obj.Y2_wpg(1,:) = -gradZdesx*Ftotn;
                obj.Y2_wpg(2,:) = gradZdesy*Ftotn;
                obj.Y2_wpg(3,:) = gradZdesx*Ftott2-gradZdesy*Ftott1; 
            else
                Ftotn = sum(simu.Fc3(3:3:end));
                Ftott1 = sum(simu.Fc3(1:3:end));
                Ftott2 = sum(simu.Fc3(2:3:end));
                Yz2 = zeros(3,obj.n_wpg);
                Yz2(1,:) = -gradZdesx*Ftotn;
                Yz2(2,:) = gradZdesy*Ftotn;
                Yz2(3,:) = gradZdesx*Ftott2-gradZdesy*Ftott1;
                for i=1:obj.n_wpg
                    At1 = reshape(obj.dPcprevious3_dwpg(:,i),3,simu.nc)*simu.Fc3(1:3:3*simu.nc);
                    At2 = reshape(obj.dPcprevious3_dwpg(:,i),3,simu.nc)*simu.Fc3(2:3:3*simu.nc);
                    An = reshape(obj.dPcprevious3_dwpg(:,i),3,simu.nc)*simu.Fc3(3:3:3*simu.nc);
                    obj.Y2_wpg(:,i) = Yz2(:,i) + [obj.t1*An; -obj.t2*An; -obj.t1*At2+obj.t2*At1];
                end
            end
        end          
        function obj = derPs3_dp(obj,sole,simu)  
            obj.dPs3_dp = zeros(3*sole.nFreeSurf,obj.lp);
            obj.dmu = zeros(3*sole.nFreeSurf,3);
            RtF = simu.R'  * reshape(simu.Fc3,3,simu.nc);
            dRthetatF = obj.dR_dtheta' * reshape(simu.Fc3,3,simu.nc);
            dRphitF = obj.dR_dphi' * reshape(simu.Fc3,3,simu.nc);
            dRpsitF = obj.dR_dpsi' * reshape(simu.Fc3,3,simu.nc);
            for j=1:sole.nFreeSurf
                P = simu.PfreeSurf3((3*j-2):3*j);
                dmu_theta = obj.dR_dtheta*P + simu.R * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(dRthetatF,3*simu.nc,1))...
                                 + obj.dR_dtheta * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(RtF,3*simu.nc,1));
                dmu_phi = obj.dR_dphi*P + simu.R * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(dRphitF,3*simu.nc,1))...
                             + obj.dR_dphi * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(RtF,3*simu.nc,1));
                dmu_psi = obj.dR_dpsi*P + simu.R * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(dRpsitF,3*simu.nc,1))...
                             + obj.dR_dpsi * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(RtF,3*simu.nc,1));
                obj.dmu(3*j-2:3*j,:) = [dmu_theta dmu_phi dmu_psi];
            end

            obj.dWs_dFs = obj.Ws(:,simu.ind_cont3)*obj.dFc3_dp;
            for i=1:obj.lp
                tmp = sole.der_Cs_dp{i}(:,simu.ind_cont3)*reshape(RtF,3*simu.nc,1);
                RdCRt = reshape(simu.R*reshape(tmp,3,sole.nFreeSurf),3*sole.nFreeSurf,1);
                obj.RdPFrees_dp = reshape(simu.R*reshape(obj.dPFrees_dp(:,i),3,sole.nFreeSurf),3*sole.nFreeSurf,1);
                obj.dPs3_dp(:,i) = repmat(obj.dOl_dp(:,i),sole.nFreeSurf,1) + obj.dmu * obj.dY_dp(:,i) + obj.RdPFrees_dp + obj.dWs_dFs(:,i) + RdCRt;
            end
        end
        function obj = derPs3_dwpg(obj,sole,simu)  
            obj.dPs3_dwpg = zeros(3*sole.nFreeSurf,obj.n_wpg);
            obj.dRot(simu);
            obj.dWs_dFs_wpg = obj.Ws(:,simu.ind_cont3)*obj.dFc3_dwpg;
            if isempty(obj.dmu)
                obj.dmu = zeros(3*sole.nFreeSurf,3);
                RtF = simu.R'  * reshape(simu.Fc3,3,simu.nc);
                dRthetatF = obj.dR_dtheta' * reshape(simu.Fc3,3,simu.nc);
                dRphitF = obj.dR_dphi' * reshape(simu.Fc3,3,simu.nc);
                dRpsitF = obj.dR_dpsi' * reshape(simu.Fc3,3,simu.nc);
                for j=1:sole.nFreeSurf
                    P = simu.PfreeSurf3((3*j-2):3*j);
                    dmu_theta = obj.dR_dtheta*P + simu.R * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(dRthetatF,3*simu.nc,1))...
                                     + obj.dR_dtheta * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(RtF,3*simu.nc,1));
                    dmu_phi = obj.dR_dphi*P + simu.R * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(dRphitF,3*simu.nc,1))...
                                 + obj.dR_dphi * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(RtF,3*simu.nc,1));
                    dmu_psi = obj.dR_dpsi*P + simu.R * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(dRpsitF,3*simu.nc,1))...
                                 + obj.dR_dpsi * (sole.Cs(3*j-2:3*j,simu.ind_cont3)*reshape(RtF,3*simu.nc,1));
                    obj.dmu(3*j-2:3*j,:) = [dmu_theta dmu_phi dmu_psi];
                end
            end
            for i=1:obj.n_wpg
                obj.dPs3_dwpg(:,i) = repmat(obj.dOl_dwpg(:,i),sole.nFreeSurf,1) + obj.dmu * obj.dY_dwpg(:,i) + obj.dWs_dFs_wpg(:,i);
            end
        end      
        function obj = dRot(obj,simu)
            obj.dR_dtheta = [0 0 0; 0 -sin(simu.theta) -cos(simu.theta); 0 cos(simu.theta) -sin(simu.theta)]* simu.Rphi * simu.Rpsi;
            obj.dR_dphi = simu.Rtheta * [-sin(simu.phi) 0 cos(simu.phi); 0 0 0; -cos(simu.phi) 0 -sin(simu.phi)] * simu.Rpsi;
            obj.dR_dpsi = simu.Rtheta * simu.Rphi * [-sin(simu.psi) -cos(simu.psi) 0; cos(simu.psi) -sin(simu.psi) 0; 0 0 0];
        end
        function obj = compute_derG_pi(obj,simu,ind_pi)
            derG21_pi = zeros(3,3*simu.nc);
            for j=1:simu.nc
                derG21_pi(:,3*j-2:3*j) = [0 0 obj.dPc3_dp(3*j-1,ind_pi);0 0 -obj.dPc3_dp(3*j-2,ind_pi); -obj.dPc3_dp(3*j-1,ind_pi) obj.dPc3_dp(3*j-2,ind_pi) 0];
            end
            derG22_pi = [];
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derG22_pi(:,2*j-1:2*j) = [0 obj.dFc3_dp(3*obj.Islip(j),ind_pi);-obj.dFc3_dp(3*obj.Islip(j),ind_pi) 0; obj.dFc3_dp(3*obj.Islip(j)-1,ind_pi) -obj.dFc3_dp(3*obj.Islip(j)-2,ind_pi)]; 
                end    
            end
            obj.derG_pi = [zeros(3,3*simu.nc) zeros(3,2*simu.ns);derG21_pi derG22_pi];
        end
        function obj = compute_derG_wpgi(obj,simu,ind_wpgi,gradZdesx_i,gradZdesy_i)
            derG21_wpgi = zeros(3,3*simu.nc);
            for j=1:simu.nc
                derG21_wpgi(:,3*j-2:3*j) = [0 0 obj.dPc3_dwpg(3*j-1,ind_wpgi)-gradZdesy_i;0 0 -obj.dPc3_dwpg(3*j-2,ind_wpgi)+gradZdesx_i; -obj.dPc3_dwpg(3*j-1,ind_wpgi)+gradZdesy_i obj.dPc3_dwpg(3*j-2,ind_wpgi)-gradZdesx_i 0];
            end
            derG22_wpgi = [];
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derG22_wpgi(:,2*j-1:2*j) = [0 obj.dFc3_dwpg(3*obj.Islip(j),ind_wpgi);-obj.dFc3_dwpg(3*obj.Islip(j),ind_wpgi) 0; obj.dFc3_dwpg(3*obj.Islip(j)-1,ind_wpgi) -obj.dFc3_dwpg(3*obj.Islip(j)-2,ind_wpgi)]; 
                end    
            end
            obj.derG_wpgi = [zeros(3,3*simu.nc) zeros(3,2*simu.ns);derG21_wpgi derG22_wpgi];
        end      
        function obj  = compute_derA_pi(obj,simu,friction,ind_pi)
            derA12_pi = zeros(3*simu.nc,2*simu.ns);
            derA21_pi = zeros(2*simu.ns,3*simu.nc);
            derA22_pi = zeros(2*simu.ns,2*simu.ns);
            obj.dR_dpi = obj.dR_dtheta * obj.dY_dp(1,ind_pi) + obj.dR_dphi * obj.dY_dp(2,ind_pi) + obj.dR_dpsi * obj.dY_dp(3,ind_pi);
            tmp = reshape(simu.R*reshape(obj.Cc,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            tmp2 = tmp';
            dR_1 = reshape(obj.dR_dpi*reshape(tmp2,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            derA11_pi = dR_1 + obj.der_HCcHt_dp{ind_pi} + dR_1';
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derA21_pi((2*j)-1:2*j,(3*obj.Islip(j))-2:3*obj.Islip(j)) = [(obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j))*obj.ddeltat2_dp(2*j-1:2*j,ind_pi)*eye(2) friction.*obj.ddeltat2_dp(2*j-1:2*j,ind_pi)];
                    a1 = obj.dFc3_dp(3*obj.Islip(j)-2:3*obj.Islip(j)-1,ind_pi)*((obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j)));
                    a2 = obj.Fc_slip(1:2,j) *(((eye(2).*norm(obj.deltaSlip_t(:,j))-obj.deltaSlip_t(:,j)*(obj.deltaSlip_t(:,j)'./norm(obj.deltaSlip_t(:,j))))./(norm(obj.deltaSlip_t(:,j))^2))*obj.ddeltat2_dp(2*j-1:2*j,ind_pi))';
                    derA22_pi((2*j)-1:2*j,(2*j)-1:2*j) = a1 + a2 + friction * obj.dFc3_dp(3*obj.Islip(j),ind_pi)*eye(2);
                end
            end
            obj.derW_pi = derA11_pi;
            obj.derA_pi = [-derA11_pi derA12_pi;derA21_pi derA22_pi];
        end
        function obj  = compute_derA_wpgi(obj,simu,friction,ind_wpgi)
            derA12_wpgi = zeros(3*simu.nc,2*simu.ns);
            derA21_wpgi = zeros(2*simu.ns,3*simu.nc);
            derA22_wpgi = zeros(2*simu.ns,2*simu.ns);
            obj.dR_dwpgi = obj.dR_dtheta * obj.dY_dwpg(1,ind_wpgi) + obj.dR_dphi * obj.dY_dwpg(2,ind_wpgi) + obj.dR_dpsi * obj.dY_dwpg(3,ind_wpgi);
            tmp = reshape(simu.R*reshape(obj.Cc,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            tmp2 = tmp';
            dR_1 = reshape(obj.dR_dwpgi*reshape(tmp2,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            derA11_wpgi = dR_1 + dR_1';
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derA21_wpgi((2*j)-1:2*j,(3*obj.Islip(j))-2:3*obj.Islip(j)) = [(obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j))*obj.ddeltat2_dwpg(2*j-1:2*j,ind_wpgi)*eye(2) friction.*obj.ddeltat2_dwpg(2*j-1:2*j,ind_wpgi)];
                    a1 = obj.dFc3_dwpg(3*obj.Islip(j)-2:3*obj.Islip(j)-1,ind_wpgi)*((obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j)));
                    a2 = obj.Fc_slip(1:2,j) *(((eye(2).*norm(obj.deltaSlip_t(:,j))-obj.deltaSlip_t(:,j)*(obj.deltaSlip_t(:,j)'./norm(obj.deltaSlip_t(:,j))))./(norm(obj.deltaSlip_t(:,j))^2))*obj.ddeltat2_dwpg(2*j-1:2*j,ind_wpgi))';
                    derA22_wpgi((2*j)-1:2*j,(2*j)-1:2*j) = a1 + a2 + friction * obj.dFc3_dwpg(3*obj.Islip(j),ind_wpgi)*eye(2);
                end
            end
            obj.derW_wpgi = derA11_wpgi;
            obj.derA_wpgi = [-derA11_wpgi derA12_wpgi;derA21_wpgi derA22_wpgi];
        end
        function obj = compute_derBw_pi(obj,simu,ind_pi)
            dFc_hat_dpi = multiSkew(obj.dFc3_dp(:,ind_pi));
            dWFh = obj.derW_pi * obj.Fc_hat;
            WdF = obj.Wc * dFc_hat_dpi;
            dWf3 = obj.derW_pi * simu.Fc3;
            a = reshape(obj.dR_dpi*simu.PfreeCont,3*simu.nc,1) + obj.RdPFreec_dp;
            ad = a + dWf3 + obj.dWc_dFc(:,ind_pi);
            derBw_pi12 = dWFh + WdF - multiSkew(ad);
            obj.derBw_pi = [zeros(3*simu.nc,3) derBw_pi12;zeros(2*simu.ns,3) zeros(2*simu.ns,3)];    
        end
        function obj = compute_derBw_wpgi(obj,simu,ind_wpgi)
            dFc_hat_dwpgi = multiSkew(obj.dFc3_dwpg(:,ind_wpgi));
            dWFh = obj.derW_wpgi * obj.Fc_hat;
            WdF = obj.Wc * dFc_hat_dwpgi;
            dWf3 = obj.derW_wpgi * simu.Fc3;
            a = reshape(obj.dR_dwpgi*simu.PfreeCont,3*simu.nc,1);
            ad = a + dWf3 + obj.dWc_dFc_wpg(:,ind_wpgi);
            derBw_wpgi12 = dWFh + WdF - multiSkew(ad);
            obj.derBw_wpgi = [zeros(3*simu.nc,3) derBw_wpgi12;zeros(2*simu.ns,3) zeros(2*simu.ns,3)];    
        end     
        function obj = compute_derXi_pi(obj,ind_pi)
            obj.derXi_pi = [0, 0, 0, 0,                       obj.dOl_dp(3,ind_pi),    -obj.dOl_dp(2,ind_pi);
                            0, 0, 0, -obj.dOl_dp(3,ind_pi),   0,                       obj.dOl_dp(1,ind_pi);
                            0, 0, 0, obj.dOl_dp(2,ind_pi),    -obj.dOl_dp(1,ind_pi),   0;
                            0, 0, 0, 0,                       0,                       0;
                            0, 0, 0, 0,                       0,                       0;
                            0, 0, 0, 0,                       0,                       0];   
        end
        function obj = compute_derXi_wpgi(obj,ind_wpgi,gradZdesx_i,gradZdesy_i)
            obj.derXi_wpgi = [0, 0, 0, 0,                                      obj.dOl_dwpg(3,ind_wpgi),                     -obj.dOl_dwpg(2,ind_wpgi)+gradZdesy_i;
                              0, 0, 0, -obj.dOl_dwpg(3,ind_wpgi),             0,                                              obj.dOl_dwpg(1,ind_wpgi)-gradZdesx_i;
                              0, 0, 0, obj.dOl_dwpg(2,ind_wpgi)-gradZdesy_i,  -obj.dOl_dwpg(1,ind_wpgi)+gradZdesx_i,      0;
                              0, 0, 0, 0,                                       0,                                              0;
                              0, 0, 0, 0,                                       0,                                              0;
                              0, 0, 0, 0,                                       0,                                              0];   
        end         
        function obj = compute_dercost_dp(obj,simu,ind_pi)
            %%% For the translational part
            derKc_pi_tr = obj.derKc_pi(1:3,1:3);
            derKc_pi_tr_sym = derKc_pi_tr' * simu.Kcart(1:3,1:3) + simu.Kcart(1:3,1:3)' * derKc_pi_tr;
            % derivative of max singular value
            dlambda1 = obj.u_max'*[derKc_pi_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda2 = obj.u_max'*[0 derKc_pi_tr_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda3 = obj.u_max'*[0 0 derKc_pi_tr_sym(1,3); 0 0 0; 0 0 0] *obj.u_max;
            dlambda4 = obj.u_max'*[0 0 0; derKc_pi_tr_sym(2,1) 0 0; 0 0 0] *obj.u_max;
            dlambda5 = obj.u_max'*[0 0 0; 0 derKc_pi_tr_sym(2,2) 0; 0 0 0] *obj.u_max;
            dlambda6 = obj.u_max'*[0 0 0; 0 0 derKc_pi_tr_sym(2,3); 0 0 0] *obj.u_max;
            dlambda7 = obj.u_max'*[0 0 0; 0 0 0; derKc_pi_tr_sym(3,1) 0 0] *obj.u_max;
            dlambda8 = obj.u_max'*[0 0 0; 0 0 0; 0 derKc_pi_tr_sym(3,2) 0] *obj.u_max;
            dlambda9 = obj.u_max'*[0 0 0; 0 0 0; 0 0 derKc_pi_tr_sym(3,3)] *obj.u_max;
            obj.der_costS_tr_dp(ind_pi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(simu.Kcart(1:3,1:3)'*simu.Kcart(1:3,1:3)))));

            %%% For the rotational part
            derKc_pi_rot = obj.derKc_pi(4:6,4:6);
            derKc_pi_rot_sym = derKc_pi_rot' * simu.Kcart(4:6,4:6) + simu.Kcart(4:6,4:6)' * derKc_pi_rot; 
            % derivative of min singular value
            dlambda1 = obj.u_min'*[derKc_pi_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda2 = obj.u_min'*[0 derKc_pi_rot_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda3 = obj.u_min'*[0 0 derKc_pi_rot_sym(1,3); 0 0 0; 0 0 0] *obj.u_min;
            dlambda4 = obj.u_min'*[0 0 0; derKc_pi_rot_sym(2,1) 0 0; 0 0 0] *obj.u_min;
            dlambda5 = obj.u_min'*[0 0 0; 0 derKc_pi_rot_sym(2,2) 0; 0 0 0] *obj.u_min;
            dlambda6 = obj.u_min'*[0 0 0; 0 0 derKc_pi_rot_sym(2,3); 0 0 0] *obj.u_min;
            dlambda7 = obj.u_min'*[0 0 0; 0 0 0; derKc_pi_rot_sym(3,1) 0 0] *obj.u_min;
            dlambda8 = obj.u_min'*[0 0 0; 0 0 0; 0 derKc_pi_rot_sym(3,2) 0] *obj.u_min;
            dlambda9 = obj.u_min'*[0 0 0; 0 0 0; 0 0 derKc_pi_rot_sym(3,3)] *obj.u_min;
            obj.der_costS_rot_dp(ind_pi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(simu.Kcart(4:6,4:6)'*simu.Kcart(4:6,4:6)))));
        end 
        function obj = compute_dercost_dwpg(obj,simu,ind_wpgi)
            %%% For the translational part
            derKc_wpgi_tr = obj.derKc_wpgi(1:3,1:3);
            derKc_wpgi_tr_sym = derKc_wpgi_tr' * simu.Kcart(1:3,1:3) + simu.Kcart(1:3,1:3)' * derKc_wpgi_tr;
            % derivative of max singular value
            dlambda1 = obj.u_max'*[derKc_wpgi_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda2 = obj.u_max'*[0 derKc_wpgi_tr_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda3 = obj.u_max'*[0 0 derKc_wpgi_tr_sym(1,3); 0 0 0; 0 0 0] *obj.u_max;
            dlambda4 = obj.u_max'*[0 0 0; derKc_wpgi_tr_sym(2,1) 0 0; 0 0 0] *obj.u_max;
            dlambda5 = obj.u_max'*[0 0 0; 0 derKc_wpgi_tr_sym(2,2) 0; 0 0 0] *obj.u_max;
            dlambda6 = obj.u_max'*[0 0 0; 0 0 derKc_wpgi_tr_sym(2,3); 0 0 0] *obj.u_max;
            dlambda7 = obj.u_max'*[0 0 0; 0 0 0; derKc_wpgi_tr_sym(3,1) 0 0] *obj.u_max;
            dlambda8 = obj.u_max'*[0 0 0; 0 0 0; 0 derKc_wpgi_tr_sym(3,2) 0] *obj.u_max;
            dlambda9 = obj.u_max'*[0 0 0; 0 0 0; 0 0 derKc_wpgi_tr_sym(3,3)] *obj.u_max;
            obj.der_costS_tr_dwpg(ind_wpgi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(simu.Kcart(1:3,1:3)'*simu.Kcart(1:3,1:3)))));

            %%% For the rotational part
            derKc_wpgi_rot = obj.derKc_wpgi(4:6,4:6);
            derKc_wpgi_rot_sym = derKc_wpgi_rot' * simu.Kcart(4:6,4:6) + simu.Kcart(4:6,4:6)' * derKc_wpgi_rot; 
            % derivative of min singular value
            dlambda1 = obj.u_min'*[derKc_wpgi_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda2 = obj.u_min'*[0 derKc_wpgi_rot_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda3 = obj.u_min'*[0 0 derKc_wpgi_rot_sym(1,3); 0 0 0; 0 0 0] *obj.u_min;
            dlambda4 = obj.u_min'*[0 0 0; derKc_wpgi_rot_sym(2,1) 0 0; 0 0 0] *obj.u_min;
            dlambda5 = obj.u_min'*[0 0 0; 0 derKc_wpgi_rot_sym(2,2) 0; 0 0 0] *obj.u_min;
            dlambda6 = obj.u_min'*[0 0 0; 0 0 derKc_wpgi_rot_sym(2,3); 0 0 0] *obj.u_min;
            dlambda7 = obj.u_min'*[0 0 0; 0 0 0; derKc_wpgi_rot_sym(3,1) 0 0] *obj.u_min;
            dlambda8 = obj.u_min'*[0 0 0; 0 0 0; 0 derKc_wpgi_rot_sym(3,2) 0] *obj.u_min;
            dlambda9 = obj.u_min'*[0 0 0; 0 0 0; 0 0 derKc_wpgi_rot_sym(3,3)] *obj.u_min;
            obj.der_costS_rot_dwpg(ind_wpgi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(simu.Kcart(4:6,4:6)'*simu.Kcart(4:6,4:6)))));
        end    
    end
end