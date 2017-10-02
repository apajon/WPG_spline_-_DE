% format long
% clc
% clear all
% close all
% 
% addpath ./input
% addpath ./FEM
% addpath ./iso2mesh
% addpath ./gradient
% 
% % Sole FEM
% fname = 'semelle1.msh';
% pname = 'input/semelle1 L=0.23, l=0.13, e=0.03 m new centre/';
% l=0.13; % width of the foot
% L=0.23; % length of the foot
% e=0.03; % thickness of the foot
% sole = SoleFEM(pname,fname,l,L,e); % create sole object with FEM properties
% 
% % B-spline to define the sole shape - it is a 2D spline
% spline_res = -pi/2:pi/10:pi/2; % spline resolution
% l_spli = length(spline_res);
% spl = Spline(sole,spline_res); % create Spline object
% 
% % Defining the type of shape Optimization as in the paper of Journal of
% % Optimization and Engineering
% type = 1;
% [move_dirichlet,first_opt,trajectory,material,w,s1,s2,s3,p_ini_v] = setConfiguration(type,l_spli);
% % Choose the sole Material defining Young's and Poisson's coefficients
% % Elastomer
% % Butyl Rubber 0.001-0.002 GPa
% % Sylicon Elastomers 0.005-0.02 GPa
% % Neoprene (CR) 0.0007-0.002 GPa
% % Neoprene (CR)
% if material == 1
%     Young = 1000000;
% else
%     Young = 1400000;    
% end
% Poisson = 0.3;
% sole.setMaterial(Young,Poisson);
% 
% % Load the ZMP trajectory comes from WPG
% % Linux
% % [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('/media/63DFD35C3FD8A3D0/Giappone/Code/sim3d/Simulation 3D - Desired ZMP - Position and Force - 3 Angles - meters - Linux/trajectory/exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
% % Windows
% if exist('trajectory/exemple_trajectoire.txt','file')~=0
%     % Linux
%     [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('trajectory/exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
% elseif exist('trajectory\exemple_trajectoire.txt','file')~=0
%     % Windows
%     [t, ZMPx, ZMPy, COMx, COMy, Fx, Fy, Fz] = textread('trajectory\exemple_trajectoire.txt', '%f %f %f %f %f %f %f %f', 'headerlines', 1);
% end
% % Change the reference of the WPG trajectory to be in the good format of
% % shape optimization movement simulation
% backtoankle=0.098; %from back to ankle of foot
% fronttoankle=0.128; %from  front to ankle of foot
% exttoankle=0.076; %from exterior to ankle of foot
% inttoankle=0.054; %from interior to ankle of foot
% xpankle=1.257728354176543; %x coordinate of ankle position
% ypankle=-0.045000000000611; %y coordinate of ankle position
% rightorleft = 1;
% [Zdes,Fdes] = changeRef(ZMPx, ZMPy, Fx, Fy, Fz,backtoankle,fronttoankle,exttoankle,inttoankle,xpankle,ypankle,rightorleft,trajectory);
% Zdes(:,end-10:end)=[];
% Zdes(:,1:10)=[];
% Zdes = Zdes(:,1:1:length(Zdes));
% Fdes(:,end-10:end)=[];
% Fdes(:,1:10)=[];
% Fdes = Fdes(:,1:1:length(Fdes));
% 
% % Create structure for the shape optimization parameters
% param_sopt = struct;
% param_sopt.friction = 0.8;
% param_sopt.Fdes = Fdes;
% param_sopt.Zdes = Zdes;
% param_sopt.move_dirichlet = move_dirichlet;
% param_sopt.w = w;
% 
% simu = Simu(); %create Simulation object 
% gradS = GradientShape(); %create GradientShape object 
% 
% pbc = opt_var(1:walking_param.nbparamtotal*2,1);
% der_cost_dopt_var = zeros(size(opt_var));
% walking_param.psa_abcdDSP=pbc;
% costZMP=1/2*pbc'*H_zmp*pbc+f_zmp'*pbc;
% 
% friction = param_sopt.friction;
% Fdes = param_sopt.Fdes;
% Zdes = param_sopt.Zdes;
% w = param_sopt.w;
% 
% %sole.coor = deformation_moveDiri(sole,p_ini_v,spl,param_sopt.move_dirichlet); % Moving Dirichlet
% st1 = 1e-04;
% % p_ini_v(5) = p_ini_v(5) + st1;
% [sole.coor,dPFree_dp] = deformation_moveDiri(sole,p_ini_v,spl,move_dirichlet);
% 
% sole.stiffness();
% sole.derStiff(length(p_ini_v),dPFree_dp);
% 
% sole.prep_stressVonMises(); % Preparation of VonMises' Stress computations
% 
% % Variable Initialisation
% simu.P = [];
% simu.no_conv = [];
% P0 = sole.coor';
% simu.PfreeSurf = P0(:,sole.nodesFreeSurf);
% simu.PfreeSurf3 = reshape(simu.PfreeSurf,3*sole.nFreeSurf,1);
% simu.displ = 1.0e-04 *[0.200672417450334;0.412701716712745;-0.813000335911768];
% simu.angleact = [0.002718435798210;-0.001201953074406;-0.000033898454071];
% simu.displ_first = [0;0;0];
% simu.psi_first = 0;
% cost_lin = 0;
% cost_rot = 0;
% 
% der_Kctr_dptot = zeros(size(dPFree_dp,2),1);
% der_Kcrot_dptot = zeros(size(dPFree_dp,2),1);
% der_costS_tr_dp_tot = zeros(size(dPFree_dp,2),1);
% der_costS_rot_dp_tot = zeros(size(dPFree_dp,2),1);
% der_costS_tr_dwpg_tot = zeros(size(gradfzmp_upd,2),1);
% der_costS_rot_dwpg_tot = zeros(size(gradfzmp_upd,2),1);
% for s=1:size(Fdes,2)
%     simu.Fdes = Fdes(:,s);
%     simu.Zdes = Zdes(:,s);
%     [simu,B_alg] = contacts(sole,simu,friction,s);
%     gradS.B = B_alg;
%     gradS.updContVar(sole,simu,friction,s);
%     gradS.derShape_pol(sole,simu,dPFree_dp,friction,s);
%     
%     sole.stressVonMises(simu.dPabstot);
%     plotsole(2,sole.elements_surf,simu.P,sole.stressVM,simu.Pc,simu.Fc,simu.Z,simu.Ftot,-37.5,30);
% 
% 
% %     cost_rot = cost_rot + sqrt(min(eig(Kcart(4:6,4:6)'*Kcart(4:6,4:6)))); % maximize the minimum stiffness. Be carefull, this is not differentiable when changing of minimum eigenvalue !!!!!!
% %     if comp_grad
% %         if s==1
% %            [der_Kctr_dp,der_Kcrot_dp,dPs_dp] = compute_gradient_step(sole,size(dPFree_dp,2),Kcart,displ,angleact,Zdes,friction,PfreeSurf,PSurf3,PabsOld,Fc_mc3,sole.Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,s,w,B_alg,[]);
% %         else
% %             [der_Kctr_dp,der_Kcrot_dp,dPs_dp] = compute_gradient_step(sole,size(dPFree_dp,2),Kcart,displ,angleact,Zdes,friction,PfreeSurf,PSurf3,PabsOld,Fc_mc3,sole.Cs,dPFree_dp,ind_cont,ind_slip,displ_first,psi_first,s,w,B_alg,dPs_dp);
% %         end
% %     end
% % 
%     if s<(size(Fdes,2)/10) %minimize linear stiffness for the first 10% of the step
% %         cost_lin = cost_lin + sqrt(max(eig(Kcart(1:3,1:3)'*Kcart(1:3,1:3))));
% %         if comp_grad
%             der_Kctr_dptot = der_Kctr_dptot + gradS.der_costS_tr_dp;
% %         end
%     end
% %     
% %     if comp_grad
%         der_Kcrot_dptot = der_Kcrot_dptot + gradS.der_costS_rot_dp;
% %     end    
% end
% derCostShape_dp = der_Kctr_dptot - w * der_Kcrot_dptot;
% 
% simu.displ = 1.0e-04 *[0.200672417450334;0.412701716712745;-0.813000335911768];
% simu.angleact = [0.002718435798210;-0.001201953074406;-0.000033898454071];
% 
% f = @(x)testder(x,sole,Fdes,Zdes,friction,spl,move_dirichlet,w);
% JOl = diff5points(f,1,1e-9,p_ini_v);     


% clear all
close all
clear walking_param
clc
format long

addpath ./cost_viapoint
addpath ./generator_zmp
addpath ./generator_com
addpath ./f_com
addpath ./torque_ankle
addpath ./divers
addpath ./Simulation_3D
addpath ./Simulation_3D/input
addpath ./Simulation_3D/results
addpath ./Simulation_3D/FEM
addpath ./Simulation_3D/iso2mesh
addpath ./Simulation_3D/gradient
addpath ./generalized_functions


%%%%%%%%%%%%%%%%%%%%%%%%%% WPG parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create walking parameters
robot=2;
type_traj=1;
firstSS=0;
frequency=200;
lambda=0.6;
epsilon=10;
%e=0.21;
e=0.3;
walking_param=new_walking_param(robot,type_traj,firstSS,frequency,lambda,epsilon,e);
%% %fix ankle positions
% load('step_number_pankle_fixed.mat');
% walking_param.choose_fixed_ankle_position(step_number_pankle_fixed)
%%
% COM init and final position
xpcominit=0.125;
xpcomfin=0.250;
ypcominit=-0.0405;
ypcomfin=0.0405;

%COM init and final position
xscominit=0.0617;
xscomfin=0.0617;
yscominit=0;
yscomfin=0;
walking_param=walking_param.choose_com_init_fin(xpcominit,ypcominit,xpcomfin,ypcomfin,xscominit,yscominit,xscomfin,yscomfin);
% generate Qp problem
[H_zmp f_zmp A_zmp b_zmp Aeq_zmp beq_zmp]=generating_problem(walking_param);

%optimization quadratic
opt = optimset('Algorithm', 'interior-point-convex','Display','iter');%,'MaxIter',1000,'TolFun',1e-6);
switch(walking_param.optim_type)
    case {1,2,3}
        [psa_abcdDSP toto]=quadprog(H_zmp,f_zmp,A_zmp,b_zmp,Aeq_zmp,-beq_zmp,[],[],[],opt);
    case 4
%         psa_abcdDSP=quadprog(AviapointDSP,BviapointDSP,AconstraintDSP_anal,BconstraintDSP_anal,AscomeqDSP_path,-Bscomeq_path,[],[],[],opt);
end
walking_param.psa_abcdDSP=psa_abcdDSP;
walking_param.nbparamABCD=walking_param.nbparamABCD+4;
%draw trajectories
drawing(walking_param,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% Shape parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = 'semelle1.msh';
%pname = 'Simulation_3D/input/semelle1 L=0.23, l=0.13, e=0.021 m new centre/';
pname = 'Simulation_3D/input/semelle1 L=0.23, l=0.13, e=0.03 m new centre/';
%%% foot size %%%
l=0.13;
L=0.23;
e=0.03;
sole = SoleFEM(pname,fname,l,L,e);

% B-spline
spline_res = -pi/2:pi/10:pi/2;
l_spli = length(spline_res);
spl = Spline(sole,spline_res);

% Defining the type of shape Optimization as in the paper of Journal of
% Optimization and Engineering
type = 1;
[move_dirichlet,first_opt,trajectory,material,w,s1,s2,s3,p_ini_v] = setConfiguration(type,l_spli);
% a = load('resultsOpt/opt_var1.mat');
% opt_var = a.opt_var;
% p_ini_v = opt_var((walking_param.nbparamtotal*2+1):end,1);

% Choose the sole Material defining Young's and Poisson's coefficients
% Elastomer
% Butyl Rubber 0.001-0.002 GPa
% Sylicon Elastomers 0.005-0.02 GPa
% Neoprene (CR) 0.0007-0.002 GPa
% Neoprene (CR)
if material == 1
    Young = 1000000;
else
    Young = 1400000;    
end

Poisson = 0.3;
sole.setMaterial(Young,Poisson);

%%% Create shape parameters struct
param_sopt = struct;
param_sopt.friction = 0.8;
param_sopt.move_dirichlet = move_dirichlet;
param_sopt.trajectory = trajectory;
param_sopt.w = w;

simu = Simu(); %create Simulation object 
gradS = GradientShape(); %create GradientShape object 

% Remesh
if first_opt==0
    [sole,spl] = remesh_auto(sole,p_ini_v,spl,param_sopt.move_dirichlet,l,L,e);
end

%%% Initialization of the optimization parameters 2
% We transform the problem to Cx = s with s < d.
pbc = walking_param.psa_abcdDSP;
s = A_zmp * pbc;

s1 = [pbc;s];

Aeq = [Aeq_zmp zeros(size(Aeq_zmp,1),length(s)) zeros(size(Aeq_zmp,1),length(p_ini_v)); A_zmp -eye(length(s),length(s)) zeros(size(A_zmp,1),length(p_ini_v))];
beq = [-beq_zmp;zeros(length(s),1)];

opt_var0 = [s1;p_ini_v];



pbc = opt_var(1:walking_param.nbparamtotal*2,1);
der_cost_dopt_var = zeros(size(opt_var));
walking_param.psa_abcdDSP=pbc;
costZMP=1/2*pbc'*H_zmp*pbc+f_zmp'*pbc;
