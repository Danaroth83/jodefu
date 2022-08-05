function x=PDFP2O(y,varargin)
% PDFP2O is a solver the problem
% min (1/2)||Ax-y||_2^2 + g(Lx)
% where g(x) is a lower semicontinuous convex function
% It implements algorithm 4 and 5 of the article:
% 'A primal–dual fixed point algorithm for convex separable minimization 
% with applications to image restoration' by Chen et al.
% The script tries the invert the input:
% y: observation variable
% These parameter are used for the inversion
% 'Nbiter': The number of iterations
% 'opA': handle for calculating the function Ax
% 'opL': handle for calculating the function Lx
% 'prox_gamma_g': handle for the proximity operator of the function \gamma*g
% If A or L are not matrices, the operator can take as input
% 'opAadj': the adjoint operator handle for opA
% 'opLadj': the adjoint operator handle for opL
% 'opnormA': The operator norm ||A||^2
% 'opnormL': The operator norm ||L||^2
% If underrelaxation is needed, it can also take as input:
% 'rho': Parameter of under relaxation
% Optional parameters:
% 'prox_gamma_g_conj': proximity operator of the conjugate of \gamma*g
% 'initx': Initialization for the parameter x

if rem(numel(varargin),2)~=0
    fprintf(['Usage: PDFP2O(y,''Nbiter'',Nbiter,''opA'',opA,''opL'',opL,',...
        '''prox_gamma_g'',prox_gamma_g,''opLadj'',opLadj,''opnormL'',opnormL)\n']);
    error('Wrong input format');
end

rho=1;
flag_opA=0;
flag_opL=0;
flag_opAadj=0;
flag_opLadj=0;
flag_opnormA=0;
flag_opnormL=0;
flag_prox_gamma_g=0;
flag_prox_gamma_g_conj=0;
flag_initx=0;
verbose=0;
opA=@(x) x;
opAadj=@(x) x;
opnormA=1;
opL=@(x) x;
opLadj=@(x) x;
opnormL=1;
Nbiter=250;

for ii=1:2:numel(varargin)
   pname=varargin{ii};
   pval=varargin{ii+1};
   if strcmpi(pname,'rho')         % Underrelaxation parameter
      rho=pval;
      if rho>1 || rho<=0, error('Rho has to be between 1 and 0'); end
   elseif strcmpi(pname,'Nbiter')  % Number of iterations
      Nbiter=pval;   
   elseif strcmpi(pname,'opA')     % Handle of operator A
      flag_opA=1;
      opA=pval;
   elseif strcmpi(pname,'opL')     % Handle of operator L
      flag_opL=1;
      opL=pval;
   elseif strcmpi(pname,'opAadj')  % Adjoint of operator A
      flag_opAadj=1;
      opAadj=pval;
   elseif strcmpi(pname,'opLadj')  % Adjoint of operator L
      flag_opLadj=1;
      opLadj=pval;
   elseif strcmpi(pname,'opnormA') % Operator norm ||A||^2
      flag_opnormA=1;
      opnormA=pval; 
   elseif strcmpi(pname,'opnormL') % Operator norm ||L||^2
      flag_opnormL=1;
      opnormL=pval;
   elseif strcmpi(pname,'prox_gamma_g') % Proximity operator of regularization function
      flag_prox_gamma_g=1;
      prox_gamma_g=pval;
   elseif strcmpi(pname,'prox_gamma_g_conj') % Proximity operator of conjugate regularization function
      flag_prox_gamma_g_conj=1;
      prox_gamma_g_conj=pval;   
   elseif strcmpi(pname,'initx')   % Initialization of x
      flag_initx=1;
      initx=pval;
   elseif strcmpi(pname,'verbose') % Verbose
      verbose=pval;
   end
end

%% Generating linear operators
if ismatrix(opA)
    if flag_opAadj==0, opAadj=opA'; end
    if flag_opnormA==0, opnormA=abs(eig(opAadj*opA,1)); end
else
    if flag_opA==0 && flag_opAadj==0 && flag_opnormA==0
        fprintf('Operator A was set to identity.\n');
    elseif flag_opAadj==0 || flag_opnormA==0
        fprintf('Information needed on adjoint operator for A.\n');
    end
end

if ismatrix(opL)
    if flag_opLadj==0, opLadj=opL'; end
    if flag_opnormL==0, opnormL=max_eig(opLadj*opL); end
else
    if flag_opL==0 && flag_opLadj==0 && flag_opnormL==0
        fprintf('Operator L was set to identity.\n');
    elseif flag_opLadj==0 || flag_opnormL==0
        fprintf('Information needed on adjoint operator for L.\n');
    end
end

%% Generating proximity operators
if flag_prox_gamma_g==0
    if flag_prox_gamma_g_conj==1
        prox_gamma_g= @(u,gamma) u-gamma*prox_gamma_g_conj(u/gamma,1/gamma);
    else
        error('Unassigned proximity operator for regularization function');
    end
end
    
I_minus_prox_gamma_g=@(u,gamma) u-prox_gamma_g(u,gamma);

%% Convergence parameters
sigma=1/opnormL;  % sigma=0.99/opnormL;
tau=1.99/opnormA;

%% Initialization
if flag_initx==0, x=opAadj(y); else, x=initx; end
u = I_minus_prox_gamma_g(opL(x),tau/sigma);

%% Iterative algorithm
for iter = 1:Nbiter
    if rho<1, x_old=x; u_old=u; end
    x= x-tau*opAadj(opA(x)-y);
    u = I_minus_prox_gamma_g(opL(x-sigma*opLadj(u))+u,tau/sigma);
    x = x-sigma*opLadj(u);
    if rho<1
        x=x_old+rho*(x-x_old);
        u=u_old+rho*(u-u_old);
    end
    if verbose==1 && mod(iter,25)==0
        fprintf('Current iteration: %d\n',iter)
    end
end