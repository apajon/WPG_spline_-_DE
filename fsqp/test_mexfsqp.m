%programme de test de mexfsqp
clear all

   mode=100;
   iprint=1;
   miter=500;  
   bigbnd=1.e10;
   eps=1.e-8;
   epsneq=0.e0;
   udelta=0.e0;
   nparam=3;
   nf=1;
   neqn=0;
   nineqn=1;
   nineq=1;
   neq=1;
   ncsrl=0;
   ncsrn=0;
   nfsr=0;
   mesh_pts(1)=0;

   bl=zeros(1,3);
   bu=bigbnd*ones(1,3);

   x(1)=0.1e0;
   x(2)=0.7e0;
   x(3)=0.2e0;
   
   cd=[0];

  
% [inform,x,f,g,lambda] = ...
%      mexfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts, ...
%               mode,iprint,miter,bigbnd,eps,epsneq,udelta,bl,bu,x, ...
%               'obj32','cntr32_bis','grob32','grcn32_bis',cd);
%               
% inform
% x
% f
% g
% lambda
ciao=1;
cd = cell(2,1);
cd{1} =  testClass(5);
cd{2} =  testClass(5);
[inform,x,f,g,lambda] = ...
     mexfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts, ...
              mode,iprint,miter,bigbnd,eps,epsneq,udelta,bl,bu,x, ...
              'obj32','cntr32_bis','grob32','grcn32_bis',cd);
              
inform
x
f
g
lambda


% T = testClass(5);
% fobj = @(j,x,cd) T.eval(x);
% fgrad = @(j,x,dummy,cd) T.diff(x);
% %fgrad = @(j,x,dummy,cd) diff(T,x);
% [inform,x,f,g,lambda] = ...
%      mexfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts, ...
%               mode,iprint,miter,bigbnd,eps,epsneq,udelta,bl,bu,x, ...
%               'fobj','cntr32_bis','fgrad','grcn32_bis',cd);
%               
% inform
% x
% f
% g
% lambda



   
   
