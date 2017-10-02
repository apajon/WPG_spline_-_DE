%---------------------------------------------------------------------
% Use of mexfsqp
%
% [inform,x,f,g,lambda] =
%    mexfsqp(nparam,nf,nfsr,nineqn,nineq,neqn,neq,ncsrl,ncsrn,mesh_pts,
%            mode,iprint,miter,bigbnd,eps,epseqn,udelta,bl,bu,x,
%            obj,constr,gradob,gradcn,cd);
%
%
%
%---------------------------------------------------------------------
% Brief specification of various arrays and parameters in the calling
% sequence. See manual for a more detailed description.
%
% nparam : number of variables
% nf     : number of objective functions (count each set of sequentially
%          related objective functions once)
% nfsr   : number of sets of sequentially related objectives (possibly
%          zero)
% nineqn : number of nonlinear inequality constraints
% nineq  : total number of inequality constraints
% neqn   : number of nonlinear equality constraints
% neq    : total number of equality constraints
% ncsrl  : number of sets of linear sequentially related inequality
%          constraints
% ncsrn  : number of sets of nonlinear sequentially related inequality
%          constraints
% mesh_pts : array of integers giving the number of actual objectives/
%            constraints in each sequentially related objective or
%            constraint set. The order is as follows:
%            (i) objective sets, (ii) nonlinear constraint sets,
%            (iii) linear constraint sets. If one or no sequentially
%            related constraint or objectives sets are present, the
%            user may simply pass the address of an integer variable
%            containing the appropriate number (possibly zero).
% mode   : mode=CBA specifies job options as described below:
%          A = 0 : ordinary minimax problems
%            = 1 : ordinary minimax problems with each individual
%                  function replaced by its absolute value, ie,
%                  an L_infty problem
%          B = 0 : monotone decrease of objective function
%                  after each iteration
%            = 1 : monotone decrease of objective function after
%                  at most four iterations
%          C = 1 : default operation.
%            = 2 : requires that constraints always be evaluated
%                  before objectives during the line search.
% iprint : print level indicator with the following options-
%          iprint=0: no normal output, only error information
%                    (this option is imposed during phase 1)
%          iprint=1: a final printout at a local solution
%          iprint=2: a brief printout at the end of each iteration
%          iprint=3: detailed infomation is printed out at the end
%                    of each iteration (for debugging purposes)
%          For iprint=2 or 3, the information may be printed at
%          iterations that are multiples of 10, instead of every
%          iteration. This may be done by adding the desired number
%          of iterations to skip printing to the desired iprint value
%          as specified above. e.g., sending iprint=23 would give
%          the iprint=3 information once every 20 iterations.
% miter  : maximum number of iterations allowed by the user to solve
%          the problem
% inform : status report at the end of execution
%          inform= 0:normal termination
%          inform= 1:no feasible point found for linear constraints
%          inform= 2:no feasible point found for nonlinear constraints
%          inform= 3:no solution has been found in miter iterations
%          inform= 4:stepsize smaller than machine precision before
%                    a successful new iterate is found
%          inform= 5:failure in attempting to construct d0
%          inform= 6:failure in attempting to construct d1
%          inform= 7:inconsistent input data
%          inform= 8:new iterate essentially identical to previous
%                    iterate, though stopping criterion not satisfied.
%          inform= 9:penalty parameter too large, unable to satisfy
%                    nonlinear equality constraint
% bigbnd : plus infinity
% eps    : stopping criterion. Execution stopped when the norm of the
%          Newton direction vector is smaller than eps
% epseqn : tolerance of the violation of nonlinear equality constraints
%          allowed by the user at an optimal solution
% udelta : perturbation size in computing gradients by finite
%          difference. The actual perturbation is determined by
%          sign(x_i) X max{udelta, rteps X max{1, |x_i|}} for each
%          component of x, where rteps is the square root of machine
%          precision.
% bl     : array of dimension nparam,containing lower bound of x
% bu     : array of dimension nparam,containing upper bound of x
% x      : array of dimension nparam,containing initial guess in input
%          and final iterate at the end of execution
% f      : array of dimension sufficient enough to hold the value of
%          all regular objective functions and the value of all
%          members of the sequentially related objective sets.
%          (dimension must be at least 1)
% g      : array of dimension sufficient enough to hold the value of
%          all regular constraint functions and the value of all
%          members of the sequentially related constraint sets.
%          (dimension must be at least 1)
% lambda : array of dimension nparam+dim(f)+dim(g), containing
%          Lagrange multiplier values at x in output. (A concerns the
%          mode, see above). The first nparam positions contain the
%          multipliers associated with the simple bounds, the next
%          dim(g) positions contain the multipliers associated with
%          the constraints. The final dim(f) positions contain the
%          multipliers associated with the objective functions. The
%          multipliers are in the order they were specified in the
%          user-defined objective and constraint functions.
% obj    : Pointer to function that returns the value of objective
%          functions, one upon each call
% constr : Pointer to function that returns the value of constraints
%          one upon each call
% gradob : Pointer to function that computes gradients of f,
%          alternatively it can be replaced by grobfd to compute
%          finite difference approximations
% gradcn : Pointer to function that computes gradients of g,
%          alternatively it can be replaced by grcnfd to compute
%          finite difference approximations
% cd     : Void pointer that may be used by the user for the passing of
%          "client data" (untouched by CFSQP)
%
%----------------------------------------------------------------------
%
%
%                       CFSQP  Version 2.5d
%
%                  Craig Lawrence, Jian L. Zhou
%                         and Andre Tits
%                  Institute for Systems Research
%                               and
%                Electrical Engineering Department
%                     University of Maryland
%                     College Park, Md 20742
%
%                         February, 1998
%
%
%  The purpose of CFSQP is to solve general nonlinear constrained
%  minimax optimization problems of the form
%
%   (A=0 in mode)     minimize    max_i f_i(x)   for i=1,...,n_f
%                        or
%   (A=1 in mode)     minimize    max_j |f_i(x)|   for i=1,...,n_f
%                       s.t.      bl   <= x <=  bu
%                                 g_j(x) <= 0,   for j=1,...,nineqn
%                                 A_1 x - B_1 <= 0
%
%                                 h_i(x)  = 0,   for i=1,...,neqn
%                                 A_2 x - B_2  = 0
%
% CFSQP is also able to efficiently handle problems with large sets of
% sequentially related objectives or constraints, see the manual for
% details.
%
%
%                  Conditions for External Use
%                  ===========================
%
%   1. The CFSQP routines may not be distributed to third parties.
%      Interested parties shall contact AEM Design directly.
%   2. If modifications are performed on the routines, these
%      modifications shall be communicated to AEM Design.  The
%      modified routines will remain the sole property of the authors.
%   3. Due acknowledgment shall be made of the use of the CFSQP
%      routines in research reports or publications. Whenever
%      such reports are released for public access, a copy shall
%      be forwarded to AEM Design.
%   4. The CFSQP routines may only be used for research and
%      development, unless it has been agreed otherwise with AEM
%      Design in writing.
%
% Copyright (c) 1993-1998 by Craig T. Lawrence, Jian L. Zhou, and
%                         Andre L. Tits
% All Rights Reserved.
