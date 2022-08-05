function x=ChambollePock(varargin)
% ChambollePock is a solver the problem
% min f(x) + g(Lx)
% where f(x) and g(x) and lower semicontinuous convex function
% The over-relaxed Chambolle-Pock algorithm is described in
% L. Condat, "A primal-dual splitting method for convex optimization 
% involving Lipschitzian, proximable and linear composite terms", 
% J. Optimization Theory and Applications, vol. 158, no. 2, 
% pp. 460-479, 2013.
% These parameter are used for the inversion:
% 'Nbiter': The number of iterations
% 'opL': handle for calculating the function Lx
% 'prox_gamma_f': handle for the proximity operator of the function \gamma*f
% 'prox_gamma_g_conj': handle for the proximity operator of the conjugate of the function \gamma*g
% If A or L are not matrices, the operator can take as input
% 'opLadj': the adjoint operator handle for opL
% 'opnormL': The operator norm ||L||^2
% The following parameter influence the convergence:
% 'tau': Convergence parameter
% 'rho': Parameter of over-relaxation
% Optional parameters:
% 'prox_gamma_f_conj': proximity operator of the conjugate of \gamma*f
% 'prox_gamma_g': proximity operator of \gamma*g
% 'initx': Initialization for the parameter x

if rem(numel(varargin),2)~=0
    fprintf(['Usage: ChambollePock(''Nbiter'',Nbiter,''opA'',opA,''opL'',opL,',...
        '''prox_gamma_f'',prox_gamma_f,''prox_gamma_g_conj'',prox_gamma_g_conj,',...
        '''opLadj'',opLadj,''opnormL'',opnormL)\n']);
    error('Wrong input format');
end

tau=0.01;
rho=1.99;
flag_opL=0;
flag_opLadj=0;
flag_opnormL=0;
flag_prox_gamma_f=0;
flag_prox_gamma_f_conj=0;
flag_prox_gamma_g=0;
flag_prox_gamma_g_conj=0;
flag_initx=0;
verbose=0;
opL=@(x) x;
opLadj=@(x) x;
opnormL=1;
Nbiter=250;

for ii=1:2:numel(varargin)
   pname=varargin{ii};
   pval=varargin{ii+1};
   if strcmpi(pname,'rho')         % Over-relaxation parameter
      rho=pval;
      if rho>2 || rho<1, error('Rho (over-relaxation parameter) has to be between 1 and 2'); end
   elseif strcmpi(pname,'tau')     % Convergence parameter
      tau=pval;
   elseif strcmpi(pname,'Nbiter')  % Number of iterations
      Nbiter=pval;
   elseif strcmpi(pname,'opL')     % Handle of operator L
      flag_opL=1;
      opL=pval;
   elseif strcmpi(pname,'opLadj')  % Adjoint of operator L
      flag_opLadj=1;
      opLadj=pval;
   elseif strcmpi(pname,'opnormL') % Operator norm ||L||^2
      flag_opnormL=1;
      opnormL=pval;
   elseif strcmpi(pname,'prox_gamma_f') % Proximity operator of error function
      flag_prox_gamma_f=1;
      prox_gamma_f=pval;
   elseif strcmpi(pname,'prox_gamma_f_conj') % Proximity operator of conjugate error function
      flag_prox_gamma_f_conj=1;
      prox_gamma_f_conj=pval; 
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
if flag_prox_gamma_g_conj==0
    if flag_prox_gamma_g==1
        prox_gamma_g_conj= @(u,gamma) u-gamma*prox_gamma_g(u/gamma,1/gamma);
    else
        error('Unassigned proximity operator for regularization function');
    end
end
if flag_prox_gamma_f==0
    if flag_prox_gamma_f_conj==1
        prox_gamma_f= @(u,gamma) u-gamma*prox_gamma_f_conj(u/gamma,1/gamma);
    else
        error('Unassigned proximity operator for error function');
    end
end

%% Convergence parameters
sigma = 1/tau/opnormL;

%% Initialization
if flag_initx==0, error('Initialization needed for x'); else, x2=initx; end
u2 = prox_gamma_g_conj(opL(x2));

%% Iterative algorithm
for iter = 1:Nbiter
    x = prox_gamma_f(x2-tau*opLadj(u2),tau);
    u = prox_gamma_g_conj(u2+sigma*opL(2*x-x2),sigma);
	x2 = x2+rho*(x-x2);
	u2 = u2+rho*(u-u2);
    if verbose==1 && mod(iter,25)==0
        fprintf('Current iteration: %d\n',iter)
    end
end