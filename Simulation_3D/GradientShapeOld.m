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
        Phi_fwpg            % see the journal paper
        Phi_zwpg            % see the journal paper
        Theta               % see the journal paper
        Y2                  % see the journal paper
        Y2_fwpg             % see the journal paper
        Y2_zwpg             % see the journal paper
        bi                  % see the journal paper
        bi_fwpg             % see the journal paper
        bi_zwpg             % see the journal paper
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
        dPs3_dfwpg          % gradient of surface node positions with respect to force wpg parameters
                            % format [dPs3_dfwpg1;dPs3_dfwpg2;...dPs3_dfwpgn_wpg]
        dPc3_dfwpg          % gradient of contact node positions with respect to force wpg parameters
                            % format [dPc3_dfwpg1;dPc3_dfwpg2;...dPc3_dfwpgn_wpg]
        dPcprevious3_dfwpg  % gradient of contact node positions at the previous time step with respect to force wpg parameters
                            % format [dPc3_dfwpg1;dPc3_dfwpg2;...dPc3_dfwpgn_wpg]                            
        dPs3_dzwpg          % gradient of surface node positions with respect to zmp position wpg parameters
                            % format [dPs3_dzwpg1;dPs3_dzwpg2;...dPs3_dzwpgn_wpg]
        dPc3_dzwpg          % gradient of contact node positions with respect to zmp position wpg parameters
                            % format [dPc3_dzwpg1;dPc3_dzwpg2;...dPc3_dzwpgn_wpg]  
        dPcprevious3_dzwpg  % gradient of contact node positions at the previous time step with respect to zmp position wpg parameters
                            % format [dPc3_dzwpg1;dPc3_dzwpg2;...dPc3_dzwpgn_wpg]                             
        der_HCcHt_dp        % derivative of HCcH' with respect to p
        der_Cc_dp           % derivative of Cc with respect to p
        t                   % matrix to select just tangential components
        t1                  % vector to select just the first tangential components
        t2                  % vector to select just the second tangential components
        dFc3_dp             % gradient of surface node forces with respect to p
        ddeltat2_dp         % gradient of deltat with respect to p (just for contact in slip condition)
        dOl_dp              % gradient of Ol(displ) with respect to p
        dY_dp               % gradient of Upsilon(vector of Angles) with respect to p
        dFc3_dfwpg          % gradient of surface node forces with respect to force wpg parameters
        ddeltat2_dfwpg      % gradient of deltat with respect to force wpg parameters (just for contact in slip condition)
        dOl_dfwpg           % gradient of Ol(displ) with respect to force wpg parameters
        dY_dfwpg            % gradient of Upsilon(vector of Angles) with respect to force wpg parameters
        dFc3_dzwpg          % gradient of surface node forces with respect to force zmp position parameters
        ddeltat2_dzwpg      % gradient of deltat with respect to zmp position wpg parameters (just for contact in slip condition)
        dOl_dzwpg           % gradient of Ol(displ) with respect to zmp position wpg parameters
        dY_dzwpg            % gradient of Upsilon(vector of Angles) with respect to zmp position wpg parameters           
        RdPFrees_dp         % R * dPFrees_dp
        RdPFreec_dp         % R * dPFreec_dp
        dmu                 % [dmu_dtheta dmu_dphi dmu_dpsi], dmu_dtheta = [dmu_dtheta_1...dmu_dtheta_m], m are the number of surface nodes
                            % with dmu_dtheta_j = dR_dtheta*P_j + sum_j^m(dR_dtheta * Cs * R' * Fj + R * Cs * dR_dtheta' * Fj)
        dR_dtheta           % derivative of rotation matrix with respect to theta
        dR_dphi             % derivative of rotation matrix with respect to phi
        dR_dpsi             % derivative of rotation matrix with respect to psi
        dWs_dFs             % Ws(:,ind_cont3)*dFc3_dp;
        dWs_dFs_fwpg        % Ws(:,ind_cont3)*dFc3_dfwpg;
        dWs_dFs_zwpg        % Ws(:,ind_cont3)*dFc3_dzwpg;
        dWc_dFc             % Wc*dFc3_dp;
        dWc_dFc_fwpg        % Wc*dFc3_dfwpg;
        dWc_dFc_zwpg        % Wc*dFc3_dzwpg;
        invABwXi            % inv_A * Bw * Xi;
        G_inv_A             % G * inv_A;
        G_inv_ABw           % G * inv_A * Bw;    
        u_max               % eigen vector of the max eigenvalue(Kc^tr)
        u_min               % eigen vector of the min eigenvalue(Kc^rot)
        derG_pi             % derivative of G with respect to pi
        derG_fwpgi          % derivative of G with respect to fwpgi
        derG_zwpgi          % derivative of G with respect to zwpgi
        derA_pi             % derivative of A with respect to pi
        derA_fwpgi          % derivative of A with respect to fwpgi
        derA_zwpgi          % derivative of A with respect to zwpgi
        derW_pi             % derivative of W with respect to pi
        derW_fwpgi          % derivative of W with respect to fwpgi
        derW_zwpgi          % derivative of W with respect to zwpgi
        dR_dpi              % derivative of R with respect to pi
        dR_dfwpgi           % derivative of R with respect to fwpgi
        dR_dzwpgi           % derivative of R with respect to zwpgi
        derBw_pi            % derivative of Bw with respect to pi
        derBw_fwpgi         % derivative of Bw with respect to fwpgi
        derBw_zwpgi         % derivative of Bw with respect to zwpgi
        derXi_pi            % derivative of Xi with respect to pi
        derXi_fwpgi         % derivative of Xi with respect to fwpgi
        derXi_zwpgi         % derivative of Xi with respect to zwpgi
        derKc_pi            % derivative of Kc with respect to pi
        derKc_fwpgi         % derivative of Kc with respect to fwpgi
        derKc_zwpgi         % derivative of Kc with respect to zwpgi
        der_costS_tr_dp     % derivative of costShape^tr with respect to p
        der_costS_rot_dp    % derivative of costShape^rot with respect to p
        der_costS_tr_dfwpg  % derivative of costShape^tr with respect to fwpg
        der_costS_rot_dfwpg % derivative of costShape^rot with respect to fwpg
        der_costS_tr_dzwpg  % derivative of costShape^tr with respect to zwpg
        der_costS_rot_dzwpg % derivative of costShape^rot with respect to zwpg         
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
                obj.Phi_fwpg = zeros(3*simu.nc,obj.n_wpg);
                obj.Phi_zwpg = zeros(3*simu.nc,obj.n_wpg);
            else           
                obj.dPcprevious3_dfwpg = obj.dPs3_dfwpg(simu.ind_cont3,:);
                obj.Phi_fwpg = obj.dPcprevious3_dfwpg;
                obj.dPcprevious3_dzwpg = obj.dPs3_dzwpg(simu.ind_cont3,:);
                obj.Phi_zwpg = obj.dPcprevious3_dzwpg;
            end
            %%% gradient of the shape cost function with respect to force wpg parameters
            obj.compute_bi_fwpg(simu,gradFdesx,gradFdesy,gradFdesz,contAngle);
            dXi_dfwpg_Fdes = obj.inv_M * obj.bi_fwpg;
            obj.dFc3_dfwpg = dXi_dfwpg_Fdes(1:3*simu.nc,:);
            obj.ddeltat2_dfwpg = dXi_dfwpg_Fdes((3*simu.nc+1):((3*simu.nc)+2*simu.ns),:);
            obj.dOl_dfwpg = dXi_dfwpg_Fdes(((3*simu.nc)+2*simu.ns)+1:((3*simu.nc)+2*simu.ns)+3,:);
            obj.dY_dfwpg = dXi_dfwpg_Fdes(((3*simu.nc)+2*simu.ns)+4:end,:);
            obj.derPs3_dfwpg(sole,simu);
            obj.dWc_dFc_fwpg = obj.dWs_dFs_fwpg(simu.ind_cont3,:);
            obj.dPc3_dfwpg = obj.dPs3_dfwpg(simu.ind_cont3,:);
            obj.der_costS_tr_dfwpg = zeros(obj.n_wpg,1);
            obj.der_costS_rot_dfwpg = zeros(obj.n_wpg,1);
            for i=1:obj.n_wpg
                obj.compute_derG_fwpgi(simu,i);
                obj.compute_derA_fwpgi(simu,friction,i);
                obj.compute_derBw_fwpgi(simu,i);
                obj.compute_derXi_fwpgi(i);
                obj.derKc_fwpgi = obj.derG_fwpgi * obj.invABwXi - obj.G_inv_A * obj.derA_fwpgi * obj.invABwXi +...
                                  obj.G_inv_A * obj.derBw_fwpgi * obj.Xi + obj.G_inv_ABw * obj.derXi_fwpgi;
                obj.compute_dercost_dfwpg(simu,i);
            end
            %%% gradient of the shape cost function with respect to zmp wpg parameters
            obj.compute_bi_zwpg(simu,gradZdesx,gradZdesy,contAngle);
            dXi_dzwpg_Zdes = obj.inv_M * obj.bi_zwpg;
            obj.dFc3_dzwpg = dXi_dzwpg_Zdes(1:3*simu.nc,:);
            obj.ddeltat2_dzwpg = dXi_dzwpg_Zdes((3*simu.nc+1):((3*simu.nc)+2*simu.ns),:);
            obj.dOl_dzwpg = dXi_dzwpg_Zdes(((3*simu.nc)+2*simu.ns)+1:((3*simu.nc)+2*simu.ns)+3,:);
            obj.dY_dzwpg = dXi_dzwpg_Zdes(((3*simu.nc)+2*simu.ns)+4:end,:);
            obj.derPs3_dzwpg(sole,simu);
            obj.dWc_dFc_zwpg = obj.dWs_dFs_zwpg(simu.ind_cont3,:);
            obj.dPc3_dzwpg = obj.dPs3_dzwpg(simu.ind_cont3,:);
            obj.der_costS_tr_dzwpg = zeros(obj.n_wpg,1);
            obj.der_costS_rot_dzwpg = zeros(obj.n_wpg,1);
            for i=1:obj.n_wpg
                obj.compute_derG_zwpgi(simu,i,gradZdesx(i),gradZdesy(i));
                obj.compute_derA_zwpgi(simu,friction,i);
                obj.compute_derBw_zwpgi(simu,i);
                obj.compute_derXi_zwpgi(i,gradZdesx(i),gradZdesy(i));
                obj.derKc_zwpgi = obj.derG_zwpgi * obj.invABwXi - obj.G_inv_A * obj.derA_zwpgi * obj.invABwXi +...
                                  obj.G_inv_A * obj.derBw_zwpgi * obj.Xi + obj.G_inv_ABw * obj.derXi_zwpgi;
                obj.compute_dercost_dzwpg(simu,i);
            end  
            
            
        end
%         function obj = derCankle_dwpg(obj,sole,simu,friction,gradZdesx,gradZdesy,gradFdesx,gradFdesy,gradFdesz,contAngle)        
%             dP_dwpg = [obj.dR_dtheta*sole.pos_ankle_loc obj.dR_dphi*sole.pos_ankle_loc obj.dR_dpsi*sole.pos_ankle_loc] * [obj.dY_dfwpg(:,1)*[gradFdesx(1);gradFdesy(1);gradFdesz(1)] obj.dY_dfwpg obj.dY_dfwpg]' * gradFdesx;
%         end
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
        function obj = compute_Y2_fwpg(obj,simu,contAngle)  
            obj.Y2_fwpg = zeros(3,obj.n_wpg);
            if contAngle>1
                for i=1:obj.n_wpg
                    At1 = reshape(obj.dPcprevious3_dfwpg(:,i),3,simu.nc)*simu.Fc3(1:3:3*simu.nc);
                    At2 = reshape(obj.dPcprevious3_dfwpg(:,i),3,simu.nc)*simu.Fc3(2:3:3*simu.nc);
                    An = reshape(obj.dPcprevious3_dfwpg(:,i),3,simu.nc)*simu.Fc3(3:3:3*simu.nc);
                    obj.Y2_fwpg(:,i) = [obj.t1*An; -obj.t2*An; -obj.t1*At2+obj.t2*At1];
                end
            end
        end
        function obj = compute_Y2_zwpg(obj,simu,contAngle,gradZdesx,gradZdesy)  
            if contAngle==1
                Ftotn = sum(simu.Fc3(3:3:end));
                Ftott1 = sum(simu.Fc3(1:3:end));
                Ftott2 = sum(simu.Fc3(2:3:end));
                obj.Y2_zwpg = zeros(3,obj.n_wpg);
                obj.Y2_zwpg(1,:) = -gradZdesx*Ftotn;
                obj.Y2_zwpg(2,:) = gradZdesy*Ftotn;
                obj.Y2_zwpg(3,:) = gradZdesx*Ftott2-gradZdesy*Ftott1; 
            else
                Ftotn = sum(simu.Fc3(3:3:end));
                Ftott1 = sum(simu.Fc3(1:3:end));
                Ftott2 = sum(simu.Fc3(2:3:end));
                Yz2 = zeros(3,obj.n_wpg);
                Yz2(1,:) = -gradZdesx*Ftotn;
                Yz2(2,:) = gradZdesy*Ftotn;
                Yz2(3,:) = gradZdesx*Ftott2-gradZdesy*Ftott1;
                for i=1:obj.n_wpg
                    At1 = reshape(obj.dPcprevious3_dzwpg(:,i),3,simu.nc)*simu.Fc3(1:3:3*simu.nc);
                    At2 = reshape(obj.dPcprevious3_dzwpg(:,i),3,simu.nc)*simu.Fc3(2:3:3*simu.nc);
                    An = reshape(obj.dPcprevious3_dzwpg(:,i),3,simu.nc)*simu.Fc3(3:3:3*simu.nc);
                    obj.Y2_zwpg(:,i) = Yz2(:,i) + [obj.t1*An; -obj.t2*An; -obj.t1*At2+obj.t2*At1];
                end
            end
        end        
        function obj = compute_bi(obj,simu)
            obj.bi = [obj.Theta;zeros(2*simu.ns,obj.lp);zeros(3,obj.lp);obj.Y2];
        end
        function obj = compute_bi_fwpg(obj,simu,gradFdesx,gradFdesy,gradFdesz,contAngle)
            if contAngle==1
                obj.bi_fwpg = [zeros(3*simu.nc,obj.n_wpg);zeros(2*simu.ns,obj.n_wpg);[gradFdesx;gradFdesy;gradFdesz];zeros(3,obj.n_wpg)];
            else
                obj.compute_Y2_fwpg(simu,contAngle);
                obj.bi_fwpg = [-obj.Phi_fwpg;zeros(2*simu.ns,obj.n_wpg);[gradFdesx;gradFdesy;gradFdesz];obj.Y2_fwpg];                
            end
        end
    
        function obj = compute_bi_zwpg(obj,simu,gradZdesx,gradZdesy,contAngle)
            obj.compute_Y2_zwpg(simu,contAngle,gradZdesx,gradZdesy);
            if contAngle==1
                obj.bi_zwpg = [zeros(3*simu.nc,obj.n_wpg);zeros(2*simu.ns,obj.n_wpg);zeros(3,obj.n_wpg);obj.Y2_zwpg];
            else            
                obj.bi_zwpg = [-obj.Phi_zwpg;zeros(2*simu.ns,obj.n_wpg);zeros(3,obj.n_wpg);obj.Y2_zwpg];                
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
        function obj = derPs3_dfwpg(obj,sole,simu)  
            obj.dPs3_dfwpg = zeros(3*sole.nFreeSurf,obj.n_wpg);
            obj.dRot(simu);
            obj.dWs_dFs_fwpg = obj.Ws(:,simu.ind_cont3)*obj.dFc3_dfwpg;
            for i=1:obj.n_wpg
                obj.dPs3_dfwpg(:,i) = repmat(obj.dOl_dfwpg(:,i),sole.nFreeSurf,1) + obj.dmu * obj.dY_dfwpg(:,i) + obj.dWs_dFs_fwpg(:,i);
            end
        end
        function obj = derPs3_dzwpg(obj,sole,simu)  
            obj.dPs3_dzwpg = zeros(3*sole.nFreeSurf,obj.n_wpg);
            obj.dRot(simu);
            obj.dWs_dFs_zwpg = obj.Ws(:,simu.ind_cont3)*obj.dFc3_dzwpg;
            for i=1:obj.n_wpg
                obj.dPs3_dzwpg(:,i) = repmat(obj.dOl_dzwpg(:,i),sole.nFreeSurf,1) + obj.dmu * obj.dY_dzwpg(:,i) + obj.dWs_dFs_zwpg(:,i);
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
        function obj = compute_derG_fwpgi(obj,simu,ind_fwpgi)
            derG21_fwpgi = zeros(3,3*simu.nc);
            for j=1:simu.nc
                derG21_fwpgi(:,3*j-2:3*j) = [0 0 obj.dPc3_dfwpg(3*j-1,ind_fwpgi);0 0 -obj.dPc3_dfwpg(3*j-2,ind_fwpgi); -obj.dPc3_dfwpg(3*j-1,ind_fwpgi) obj.dPc3_dfwpg(3*j-2,ind_fwpgi) 0];
            end
            derG22_fwpgi = [];
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derG22_fwpgi(:,2*j-1:2*j) = [0 obj.dFc3_dfwpg(3*obj.Islip(j),ind_fwpgi);-obj.dFc3_dfwpg(3*obj.Islip(j),ind_fwpgi) 0; obj.dFc3_dfwpg(3*obj.Islip(j)-1,ind_fwpgi) -obj.dFc3_dfwpg(3*obj.Islip(j)-2,ind_fwpgi)]; 
                end    
            end
            obj.derG_fwpgi = [zeros(3,3*simu.nc) zeros(3,2*simu.ns);derG21_fwpgi derG22_fwpgi];
        end
        function obj = compute_derG_zwpgi(obj,simu,ind_zwpgi,gradZdesx_i,gradZdesy_i)
            derG21_zwpgi = zeros(3,3*simu.nc);
            for j=1:simu.nc
                derG21_zwpgi(:,3*j-2:3*j) = [0 0 obj.dPc3_dzwpg(3*j-1,ind_zwpgi)-gradZdesy_i;0 0 -obj.dPc3_dzwpg(3*j-2,ind_zwpgi)+gradZdesx_i; -obj.dPc3_dzwpg(3*j-1,ind_zwpgi)+gradZdesy_i obj.dPc3_dzwpg(3*j-2,ind_zwpgi)-gradZdesx_i 0];
            end
            derG22_zwpgi = [];
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derG22_zwpgi(:,2*j-1:2*j) = [0 obj.dFc3_dzwpg(3*obj.Islip(j),ind_zwpgi);-obj.dFc3_dzwpg(3*obj.Islip(j),ind_zwpgi) 0; obj.dFc3_dzwpg(3*obj.Islip(j)-1,ind_zwpgi) -obj.dFc3_dzwpg(3*obj.Islip(j)-2,ind_zwpgi)]; 
                end    
            end
            obj.derG_zwpgi = [zeros(3,3*simu.nc) zeros(3,2*simu.ns);derG21_zwpgi derG22_zwpgi];
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
        function obj  = compute_derA_fwpgi(obj,simu,friction,ind_fwpgi)
            derA12_fwpgi = zeros(3*simu.nc,2*simu.ns);
            derA21_fwpgi = zeros(2*simu.ns,3*simu.nc);
            derA22_fwpgi = zeros(2*simu.ns,2*simu.ns);
            obj.dR_dfwpgi = obj.dR_dtheta * obj.dY_dfwpg(1,ind_fwpgi) + obj.dR_dphi * obj.dY_dfwpg(2,ind_fwpgi) + obj.dR_dpsi * obj.dY_dfwpg(3,ind_fwpgi);
            tmp = reshape(simu.R*reshape(obj.Cc,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            tmp2 = tmp';
            dR_1 = reshape(obj.dR_dfwpgi*reshape(tmp2,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            derA11_fwpgi = dR_1 + dR_1';
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derA21_fwpgi((2*j)-1:2*j,(3*obj.Islip(j))-2:3*obj.Islip(j)) = [(obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j))*obj.ddeltat2_dfwpg(2*j-1:2*j,ind_fwpgi)*eye(2) friction.*obj.ddeltat2_dfwpg(2*j-1:2*j,ind_fwpgi)];
                    a1 = obj.dFc3_dfwpg(3*obj.Islip(j)-2:3*obj.Islip(j)-1,ind_fwpgi)*((obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j)));
                    a2 = obj.Fc_slip(1:2,j) *(((eye(2).*norm(obj.deltaSlip_t(:,j))-obj.deltaSlip_t(:,j)*(obj.deltaSlip_t(:,j)'./norm(obj.deltaSlip_t(:,j))))./(norm(obj.deltaSlip_t(:,j))^2))*obj.ddeltat2_dfwpg(2*j-1:2*j,ind_fwpgi))';
                    derA22_fwpgi((2*j)-1:2*j,(2*j)-1:2*j) = a1 + a2 + friction * obj.dFc3_dfwpg(3*obj.Islip(j),ind_fwpgi)*eye(2);
                end
            end
            obj.derW_fwpgi = derA11_fwpgi;
            obj.derA_fwpgi = [-derA11_fwpgi derA12_fwpgi;derA21_fwpgi derA22_fwpgi];
        end
        function obj  = compute_derA_zwpgi(obj,simu,friction,ind_zwpgi)
            derA12_zwpgi = zeros(3*simu.nc,2*simu.ns);
            derA21_zwpgi = zeros(2*simu.ns,3*simu.nc);
            derA22_zwpgi = zeros(2*simu.ns,2*simu.ns);
            obj.dR_dzwpgi = obj.dR_dtheta * obj.dY_dzwpg(1,ind_zwpgi) + obj.dR_dphi * obj.dY_dzwpg(2,ind_zwpgi) + obj.dR_dpsi * obj.dY_dzwpg(3,ind_zwpgi);
            tmp = reshape(simu.R*reshape(obj.Cc,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            tmp2 = tmp';
            dR_1 = reshape(obj.dR_dzwpgi*reshape(tmp2,3,3*simu.nc*simu.nc),3*simu.nc,3*simu.nc);
            derA11_zwpgi = dR_1 + dR_1';
            if ~isempty(simu.ind_slip)
                for j=1:simu.ns
                    derA21_zwpgi((2*j)-1:2*j,(3*obj.Islip(j))-2:3*obj.Islip(j)) = [(obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j))*obj.ddeltat2_dzwpg(2*j-1:2*j,ind_zwpgi)*eye(2) friction.*obj.ddeltat2_dzwpg(2*j-1:2*j,ind_zwpgi)];
                    a1 = obj.dFc3_dzwpg(3*obj.Islip(j)-2:3*obj.Islip(j)-1,ind_zwpgi)*((obj.deltaSlip_t(:,j))'./norm(obj.deltaSlip_t(:,j)));
                    a2 = obj.Fc_slip(1:2,j) *(((eye(2).*norm(obj.deltaSlip_t(:,j))-obj.deltaSlip_t(:,j)*(obj.deltaSlip_t(:,j)'./norm(obj.deltaSlip_t(:,j))))./(norm(obj.deltaSlip_t(:,j))^2))*obj.ddeltat2_dzwpg(2*j-1:2*j,ind_zwpgi))';
                    derA22_zwpgi((2*j)-1:2*j,(2*j)-1:2*j) = a1 + a2 + friction * obj.dFc3_dzwpg(3*obj.Islip(j),ind_zwpgi)*eye(2);
                end
            end
            obj.derW_zwpgi = derA11_zwpgi;
            obj.derA_zwpgi = [-derA11_zwpgi derA12_zwpgi;derA21_zwpgi derA22_zwpgi];
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
        function obj = compute_derBw_fwpgi(obj,simu,ind_fwpgi)
            dFc_hat_dfwpgi = multiSkew(obj.dFc3_dfwpg(:,ind_fwpgi));
            dWFh = obj.derW_fwpgi * obj.Fc_hat;
            WdF = obj.Wc * dFc_hat_dfwpgi;
            dWf3 = obj.derW_fwpgi * simu.Fc3;
            a = reshape(obj.dR_dfwpgi*simu.PfreeCont,3*simu.nc,1);
            ad = a + dWf3 + obj.dWc_dFc_fwpg(:,ind_fwpgi);
            derBw_fwpgi12 = dWFh + WdF - multiSkew(ad);
            obj.derBw_fwpgi = [zeros(3*simu.nc,3) derBw_fwpgi12;zeros(2*simu.ns,3) zeros(2*simu.ns,3)];    
        end
        function obj = compute_derBw_zwpgi(obj,simu,ind_zwpgi)
            dFc_hat_dzwpgi = multiSkew(obj.dFc3_dzwpg(:,ind_zwpgi));
            dWFh = obj.derW_zwpgi * obj.Fc_hat;
            WdF = obj.Wc * dFc_hat_dzwpgi;
            dWf3 = obj.derW_zwpgi * simu.Fc3;
            a = reshape(obj.dR_dzwpgi*simu.PfreeCont,3*simu.nc,1);
            ad = a + dWf3 + obj.dWc_dFc_zwpg(:,ind_zwpgi);
            derBw_zwpgi12 = dWFh + WdF - multiSkew(ad);
            obj.derBw_zwpgi = [zeros(3*simu.nc,3) derBw_zwpgi12;zeros(2*simu.ns,3) zeros(2*simu.ns,3)];    
        end        
        function obj = compute_derXi_pi(obj,ind_pi)
            obj.derXi_pi = [0, 0, 0, 0,                       obj.dOl_dp(3,ind_pi),    -obj.dOl_dp(2,ind_pi);
                            0, 0, 0, -obj.dOl_dp(3,ind_pi),   0,                       obj.dOl_dp(1,ind_pi);
                            0, 0, 0, obj.dOl_dp(2,ind_pi),    -obj.dOl_dp(1,ind_pi),   0;
                            0, 0, 0, 0,                       0,                       0;
                            0, 0, 0, 0,                       0,                       0;
                            0, 0, 0, 0,                       0,                       0];   
        end
        function obj = compute_derXi_fwpgi(obj,ind_fwpgi)
            obj.derXi_fwpgi = [0, 0, 0, 0,                          obj.dOl_dfwpg(3,ind_fwpgi),  -obj.dOl_dfwpg(2,ind_fwpgi);
                              0, 0, 0, -obj.dOl_dfwpg(3,ind_fwpgi), 0,                           obj.dOl_dfwpg(1,ind_fwpgi);
                              0, 0, 0, obj.dOl_dfwpg(2,ind_fwpgi),  -obj.dOl_dfwpg(1,ind_fwpgi), 0;
                              0, 0, 0, 0,                           0,                           0;
                              0, 0, 0, 0,                           0,                           0;
                              0, 0, 0, 0,                           0,                           0];   
        end 
        function obj = compute_derXi_zwpgi(obj,ind_zwpgi,gradZdesx_i,gradZdesy_i)
            obj.derXi_zwpgi = [0, 0, 0, 0,                                      obj.dOl_dzwpg(3,ind_zwpgi),                     -obj.dOl_dzwpg(2,ind_zwpgi)+gradZdesy_i;
                              0, 0, 0, -obj.dOl_dzwpg(3,ind_zwpgi),             0,                                              obj.dOl_dzwpg(1,ind_zwpgi)-gradZdesx_i;
                              0, 0, 0, obj.dOl_dzwpg(2,ind_zwpgi)-gradZdesy_i,  -obj.dOl_dzwpg(1,ind_zwpgi)+gradZdesx_i,      0;
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
        function obj = compute_dercost_dfwpg(obj,simu,ind_fwpgi)
            %%% For the translational part
            derKc_fwpgi_tr = obj.derKc_fwpgi(1:3,1:3);
            derKc_fwpgi_tr_sym = derKc_fwpgi_tr' * simu.Kcart(1:3,1:3) + simu.Kcart(1:3,1:3)' * derKc_fwpgi_tr;
            % derivative of max singular value
            dlambda1 = obj.u_max'*[derKc_fwpgi_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda2 = obj.u_max'*[0 derKc_fwpgi_tr_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda3 = obj.u_max'*[0 0 derKc_fwpgi_tr_sym(1,3); 0 0 0; 0 0 0] *obj.u_max;
            dlambda4 = obj.u_max'*[0 0 0; derKc_fwpgi_tr_sym(2,1) 0 0; 0 0 0] *obj.u_max;
            dlambda5 = obj.u_max'*[0 0 0; 0 derKc_fwpgi_tr_sym(2,2) 0; 0 0 0] *obj.u_max;
            dlambda6 = obj.u_max'*[0 0 0; 0 0 derKc_fwpgi_tr_sym(2,3); 0 0 0] *obj.u_max;
            dlambda7 = obj.u_max'*[0 0 0; 0 0 0; derKc_fwpgi_tr_sym(3,1) 0 0] *obj.u_max;
            dlambda8 = obj.u_max'*[0 0 0; 0 0 0; 0 derKc_fwpgi_tr_sym(3,2) 0] *obj.u_max;
            dlambda9 = obj.u_max'*[0 0 0; 0 0 0; 0 0 derKc_fwpgi_tr_sym(3,3)] *obj.u_max;
            obj.der_costS_tr_dfwpg(ind_fwpgi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(simu.Kcart(1:3,1:3)'*simu.Kcart(1:3,1:3)))));

            %%% For the rotational part
            derKc_fwpgi_rot = obj.derKc_fwpgi(4:6,4:6);
            derKc_fwpgi_rot_sym = derKc_fwpgi_rot' * simu.Kcart(4:6,4:6) + simu.Kcart(4:6,4:6)' * derKc_fwpgi_rot; 
            % derivative of min singular value
            dlambda1 = obj.u_min'*[derKc_fwpgi_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda2 = obj.u_min'*[0 derKc_fwpgi_rot_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda3 = obj.u_min'*[0 0 derKc_fwpgi_rot_sym(1,3); 0 0 0; 0 0 0] *obj.u_min;
            dlambda4 = obj.u_min'*[0 0 0; derKc_fwpgi_rot_sym(2,1) 0 0; 0 0 0] *obj.u_min;
            dlambda5 = obj.u_min'*[0 0 0; 0 derKc_fwpgi_rot_sym(2,2) 0; 0 0 0] *obj.u_min;
            dlambda6 = obj.u_min'*[0 0 0; 0 0 derKc_fwpgi_rot_sym(2,3); 0 0 0] *obj.u_min;
            dlambda7 = obj.u_min'*[0 0 0; 0 0 0; derKc_fwpgi_rot_sym(3,1) 0 0] *obj.u_min;
            dlambda8 = obj.u_min'*[0 0 0; 0 0 0; 0 derKc_fwpgi_rot_sym(3,2) 0] *obj.u_min;
            dlambda9 = obj.u_min'*[0 0 0; 0 0 0; 0 0 derKc_fwpgi_rot_sym(3,3)] *obj.u_min;
            obj.der_costS_rot_dfwpg(ind_fwpgi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(simu.Kcart(4:6,4:6)'*simu.Kcart(4:6,4:6)))));
        end
        function obj = compute_dercost_dzwpg(obj,simu,ind_zwpgi)
            %%% For the translational part
            derKc_zwpgi_tr = obj.derKc_zwpgi(1:3,1:3);
            derKc_zwpgi_tr_sym = derKc_zwpgi_tr' * simu.Kcart(1:3,1:3) + simu.Kcart(1:3,1:3)' * derKc_zwpgi_tr;
            % derivative of max singular value
            dlambda1 = obj.u_max'*[derKc_zwpgi_tr_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda2 = obj.u_max'*[0 derKc_zwpgi_tr_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_max;
            dlambda3 = obj.u_max'*[0 0 derKc_zwpgi_tr_sym(1,3); 0 0 0; 0 0 0] *obj.u_max;
            dlambda4 = obj.u_max'*[0 0 0; derKc_zwpgi_tr_sym(2,1) 0 0; 0 0 0] *obj.u_max;
            dlambda5 = obj.u_max'*[0 0 0; 0 derKc_zwpgi_tr_sym(2,2) 0; 0 0 0] *obj.u_max;
            dlambda6 = obj.u_max'*[0 0 0; 0 0 derKc_zwpgi_tr_sym(2,3); 0 0 0] *obj.u_max;
            dlambda7 = obj.u_max'*[0 0 0; 0 0 0; derKc_zwpgi_tr_sym(3,1) 0 0] *obj.u_max;
            dlambda8 = obj.u_max'*[0 0 0; 0 0 0; 0 derKc_zwpgi_tr_sym(3,2) 0] *obj.u_max;
            dlambda9 = obj.u_max'*[0 0 0; 0 0 0; 0 0 derKc_zwpgi_tr_sym(3,3)] *obj.u_max;
            obj.der_costS_tr_dzwpg(ind_zwpgi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(max(eig(simu.Kcart(1:3,1:3)'*simu.Kcart(1:3,1:3)))));

            %%% For the rotational part
            derKc_zwpgi_rot = obj.derKc_zwpgi(4:6,4:6);
            derKc_zwpgi_rot_sym = derKc_zwpgi_rot' * simu.Kcart(4:6,4:6) + simu.Kcart(4:6,4:6)' * derKc_zwpgi_rot; 
            % derivative of min singular value
            dlambda1 = obj.u_min'*[derKc_zwpgi_rot_sym(1,1) 0 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda2 = obj.u_min'*[0 derKc_zwpgi_rot_sym(1,2) 0; 0 0 0; 0 0 0] *obj.u_min;
            dlambda3 = obj.u_min'*[0 0 derKc_zwpgi_rot_sym(1,3); 0 0 0; 0 0 0] *obj.u_min;
            dlambda4 = obj.u_min'*[0 0 0; derKc_zwpgi_rot_sym(2,1) 0 0; 0 0 0] *obj.u_min;
            dlambda5 = obj.u_min'*[0 0 0; 0 derKc_zwpgi_rot_sym(2,2) 0; 0 0 0] *obj.u_min;
            dlambda6 = obj.u_min'*[0 0 0; 0 0 derKc_zwpgi_rot_sym(2,3); 0 0 0] *obj.u_min;
            dlambda7 = obj.u_min'*[0 0 0; 0 0 0; derKc_zwpgi_rot_sym(3,1) 0 0] *obj.u_min;
            dlambda8 = obj.u_min'*[0 0 0; 0 0 0; 0 derKc_zwpgi_rot_sym(3,2) 0] *obj.u_min;
            dlambda9 = obj.u_min'*[0 0 0; 0 0 0; 0 0 derKc_zwpgi_rot_sym(3,3)] *obj.u_min;
            obj.der_costS_rot_dzwpg(ind_zwpgi) = (dlambda1+dlambda2+dlambda3+dlambda4+dlambda5+dlambda6+dlambda7+dlambda8+dlambda9)/(2*sqrt(min(eig(simu.Kcart(4:6,4:6)'*simu.Kcart(4:6,4:6)))));
        end         
    end
end