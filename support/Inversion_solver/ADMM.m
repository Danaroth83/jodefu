function x=ADMM(varargin)
% ADMM is a solver the problem
% min f(x) + g(x)
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
    fprintf(['Usage: ADMM(''Nbiter'',Nbiter,''tau'',tau,''rho'',rho',...
        '''prox_gamma_f'',prox_gamma_f,''prox_gamma_g_conj'',prox_gamma_g_conj)\n']);
    error('Wrong input format');
end

rho=0.99;
tau=0.01;
flag_prox_gamma_f=0;
flag_prox_gamma_f_conj=0;
flag_prox_gamma_g=0;
flag_prox_gamma_g_conj=0;
flag_initx=0;
verbose=0;
Nbiter=250;

for ii=1:2:numel(varargin)
   pname=varargin{ii};
   pval=varargin{ii+1};
   if strcmpi(pname,'rho')         % Under-relaxation parameter
      rho=pval;
      if rho>1 || rho<=0, error('Rho (under-relaxation parameter) has to be between 0 and 1'); end
   elseif strcmpi(pname,'tau')     % Convergence parameter
      tau=pval;
      if tau<=0, error('Tau (Convergence parameter) has to be positive'); end
   elseif strcmpi(pname,'Nbiter')  % Number of iterations
      Nbiter=pval;
   elseif strcmpi(pname,'prox_gamma_f')      % Proximity operator of error function
      flag_prox_gamma_f=1;
      prox_gamma_f=pval;
   elseif strcmpi(pname,'prox_gamma_f_conj') % Proximity operator of conjugate error function
      flag_prox_gamma_f_conj=1;
      prox_gamma_f_conj=pval; 
   elseif strcmpi(pname,'prox_gamma_g')      % Proximity operator of regularization function
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

%% Initialization
if flag_initx==0, error('Initialization needed for x'); else x=initx; end
u = 2*rho*(prox_gamma_g_conj(2*x,tau)-x);

%% Iterative algorithm
for iter = 1:Nbiter
    x = prox_gamma_f(u2,tau);
    u = u + 2*rho*(prox_gamma_g_conj(2*x-u,tau)-x);
    if verbose==1 && mod(iter,25)==0
        fprintf('Current iteration: %d\n',iter)
    end
end