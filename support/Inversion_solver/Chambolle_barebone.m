%% CHAMBOLLE-POCK INVERSION SOLVER
%
% (copyright) Daniele Picone
% version 8.0
% 
% Description:
% Loris is a solver the problem
%       argmin_x   f(Ax) + <x,c> + lambda g(Lx)
% where: A and L are bounded linear operators
%        f(x) and g(x) are lower semicontinuous convex function
%        lambda>0 is a regularization parameter
% It implements  the generalized Chambolle-Pock algorithm (22) in [1] and
% the over-relaxation of eq.53 in [2]. Also known as PDFP2O [3]. This
% implementation uses regularization paths [5], using the results of
% previous convergences for a certain lambda as warm start for following
% ones
% 
% Usage:
% x=Chambolle('Field',Value);
%
% Example:
% im=double(imread('peppers.png'))/255;
% y=im+0.1*randn(size(im));
% opf.direct= @(x) sum(x(:).^2)/2;
% opf.proxconj=@(x,gamma) 1/(1+gamma)*x;
% opA.direct= @(x) x;
% opA.adjoint= @(x) x;
% opA.norm=1;
% opL.direct = @(x) cat(4,[diff(x,1,1);zeros(1,size(x,2),size(x,3))],[diff(x,1,2) zeros(size(x,1),1,size(x,3))]);
% opL.adjoint = @(u) -[u(1,:,:,1);diff(u(:,:,:,1),1,1)]-[u(:,1,:,2) diff(u(:,:,:,2),1,2)];
% opL.norm=8;
% opy = @(x) x ./ max(sqrt(sum(x.^2,2:3)),1);
% opg.proxconj= @(x,~) reshape(opy(reshape(x,[],size(x,3),size(x,4))),size(x));
% opg.direct = @(x) sum(sqrt(sum(reshape(x,[],size(x,3),size(x,4)).^2,2:3)),1);
% 
% [x,cost]=Loris('y',y,'Nbiter',250,'opf',opf,'opA',opA,'opL',opL,'reg',opg,...
%                'lambda',0.1:0.01:0.2,'rho',1,'tol',10E-6);
%
% Fields:
% 'c': Scalar product term in the cost function
% 'y': If this is provided in place of above, this problem is solved instead:
%          argmin_x   f(Ax) - <Ax,y> + lambda g(Lx)
% 'lambda': Array of regularization parameter (Array of scalars >0 - default: 0.1)
% 'opA': struct for calculating the Direct Tranfer Function Ax, whose fields are:
%      - direct:  Handle for the direct tranfer model Ax
%      - adjoint: Handle for the adjoint operator of the above
%      - norm:    Operator norm ||A||^2
% 'opL': struct for calculating Lx, whose fields are:
%      - direct:  Handle for calculating Lx
%      - adjoint: Handle for the adjoint operator of the above
%      - norm:    Operator norm ||L||^2
% 'opW': struct for preprocessing of data, whose fields are
%      - direct:  Handle for calculating the preprocessing
%      - adjoint: Handle for the adjoint operator of the above
%      - norm:    Operator norm ||W||^2
%      If this operator is provided, the algorithm solves the following
%      minimization problem:
%          W^(-1) { argmin_x 1/(2 lambda)||(AWx-y)||_2^2 + g(Lx)) }
% 'opg': struct relative to the operators relative to the function g(x), whose fields are:
%      - proxconj: handle for the proximal operator; it has to have three 
%                  inputs (x,gamma,lambda) to calculate the proximal
%                  operator of the function gamma multiplied by the convex
%                  conjugate of lambda*g(x)
%      - prox:     handle for the proximal operator of gamma*g(x); it has 
%                  two inputs (x,gamma) and generates the handle proxconj
%                  via Moreau identity if the above is not provided
%      - direct:   handle to calculate g(x)
% 'opf': struct relative to the operators relative to the function f(x), whose fields are:
%      - proxconj: handle for the proximal operator; it has to have three 
%                  inputs (x,gamma,lambda) to calculate the proximal
%                  operator of the function gamma multiplied by the convex
%                  conjugate of lambda*f(x)
%      - prox:     handle for the proximal operator of gamma*g(x); it has 
%                  two inputs (x,gamma) and generates the handle proxconj
%                  via Moreau identity if the above is not provided
%      - direct:   handle to calculate f(x)
% 'Nbiter': Max number of iterations (default: 1000).
% 'init': Initialization for the solution x (Default: adjoint transformation applied on y)
% 'rho': Over-relaxation parameter (Must be 1<=rho<2, default: 1)
% 'tol': difference in cost value to stop iterations, scaled to the
%              norm of y (If 0, check is ignored. Default: 1E-7)
% 'verbose': If 1, prints the iteration status (Default: 1)
% 
% References:
% [1] Loris I. and Verhoeven C., "On a generalization of the iterative
% soft-thresholding algorithm for the case of non-separable penalty."
% Inverse Problems 27.12 (2011): 125007.
% [2] Condat L.,  Kitahara D., Contreras A. and Hirabayashi A., "Proximal
% splitting algorithms: Relax them all!" (preprint)
% [3] Chen P., Huang J., and Zhang X., "A primal–dual fixed point algorithm
% for convex separable minimization with applications to image 
% restoration." Inverse Problems 29.2 (2013): 025011.
% [4] Borwein J. and Lewis A. S., "Convex analysis and nonlinear 
% optimization: theory and examples." Springer Science & Business Media, 
% 2010.
% [5] Jerome F., Hastie T., and Tibshirani R., "Regularization paths for 
% generalized linear models via coordinate descent." Journal of statistical
% software 33.1 (2010): 1.

function [x_out,cost_out]=Chambolle_barebone(varargin)

%% Parse inputs

if rem(numel(varargin),2)~=0
    fprintf(['Usage: [x,cost]=Chambolle(''y'',y,''Nbiter'',Nbiter,''opA'',opA,''opL'',opL,','''prox'',prox)\n']);
    error('Wrong input format');
end

iter_save=1000;
lambda=0.1;
verbose=1;
opA=[];
opL=[];
opW=[];
rho=1;
prox_gamma_f_conj=[];
prox_gamma_g_conj=[];
initx=[];
tau_init=[];
deltacost_max=10E-7;
Npos=3; % Number of consecutive correct checks to stop the iterations
c=[];
y=[];

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if any(strcmpi(pname,{'c'})) 
        c=pval;
    elseif any(strcmpi(pname,{'y'}))
        y=pval;
    elseif any(strcmpi(pname,{'Nbiter','iter_save'})) % Number of iterations
        iter_save=pval;
    elseif strcmpi(pname,'lambda')  % Regularization Parameter
        lambda=pval;   
    elseif any(strcmpi(pname,{'A','opA'}))     % Handle of operator A
        if isnumeric(pval), opA=pval; end
        if isfield(pval,'direct')
            opA=pval.direct;
            if isfield(pval,'adjoint'), opAadj=pval.adjoint; else, error('Missing information on adjoint operator'); end
            if isfield(pval,'norm'),    opAnorm=pval.norm; else, error('Missing information on operator norm'); end
        end
        if isnumeric(opA)
            opAnorm=max(svd(opA));
            opAadj=@(x) opA'*x;
            opA=@(x) opA*x; 
        end
    elseif any(strcmpi(pname,{'L','D','opL','opD'}))     % Handle of operator L
        if isnumeric(pval), opL=pval; end
        if isfield(pval,'direct')
            opL=pval.direct;
            if isfield(pval,'adjoint'), opLadj=pval.adjoint; else, error('Missing information on adjoint operator'); end
            if isfield(pval,'norm'), opLnorm=pval.norm; else, error('Missing information on operator norm'); end
        end
        if isnumeric(opL)
            opLnorm=max(svd(opL));
            opLadj=@(x) opL'*x;
            opL=@(x) opL*x; 
        end
        elseif any(strcmpi(pname,{'W','opW'}))     % Handle of operator W
        if isnumeric(pval), opW=pval; end
        if isfield(pval,'direct')
            opW=pval.direct;
            if isfield(pval,'adjoint'), opWadj=pval.adjoint; else, error('Missing information on adjoint operator'); end
            if isfield(pval,'norm'), opWnorm=pval.norm; else, error('Missing information on operator norm'); end
            if isfield(pval,'inverse'), opWinv=pval.inverse; else, opWinv=opWadj; warning('Assumption: the transfomation W is orthogonal'); end
        end
        if isnumeric(opW)
            opWnorm=max(svd(opW));
            opWadj=@(x) opW'*x;
            opWinv=@(x) opW\x;
            opW=@(x) opW*x; 
        end
    elseif any(strcmpi(pname,{'f','opf','fidelity'})) % Proximal operator for function g
        if isfield(pval,'proxconj')
            if nargin(pval.proxconj)==3
                prox_gamma_f_conj=@(x,gamma) pval.proxconj(x,gamma,1);
            elseif nargin(pval.proxconj)==2
                prox_gamma_f_conj=pval.proxconj;
            else
                error('The proximal operator is not properly formatted');
            end
        elseif isfield(pval,'prox')
            % Obtaining the conjugate proximal operator via Moreau identity
            if nargin(pval.prox)==3,prox_gamma_f_conj=@(u,gamma) u-gamma*pval.prox(u/gamma,1/gamma,1);
            elseif nargin(pval.prox)==2, prox_gamma_f_conj=@(u,gamma) u-gamma*pval.prox(u/gamma,1/gamma);
            else, error('The proximal operator is not properly formatted');
            end
        else
            error('The proximal operator is not properly formatted');
        end
        if isfield(pval,'direct')
            if nargin(pval.direct)==2, norm_handle_f=@(x) pval.direct(x,1);
            elseif nargin(pval.direct)==1, norm_handle_f=pval.direct;
            else, error('Not able to parse the norm operator');
            end
        elseif nargout>=1
            error('Regularization function handle is needed');
        end
    elseif any(strcmpi(pname,{'g','opg','prox','norm','proximal','regularizer','reg'})) % Proximal operator for function g
        if isfield(pval,'proxconj')
            if nargin(pval.proxconj)==3
                prox_gamma_g_conj= pval.proxconj;
            else
                error('The proximal operator is not properly formatted');
            end
        elseif isfield(pval,'prox')
            % Obtaining the conjugate proximal operator via Moreau identity
            if nargin(pval.prox)==3,prox_gamma_g_conj=@(u,gamma) u-gamma*pval.prox(u/gamma,1/gamma,lambda);
            elseif nargin(pval.prox)==2, prox_gamma_g_conj=@(u,gamma,lambda) u-gamma*pval.prox(u/gamma,lambda/gamma);
            else, error('The proximal operator is not properly formatted');
            end
        else
            error('The proximal operator is not properly formatted');
        end
        if isfield(pval,'direct')
            if nargin(pval.direct)==2, norm_handle_g=pval.direct;
            elseif nargin(pval.direct)==1, norm_handle_g=@(x,lambda) lambda*pval.direct(x);
            else, error('Not able to parse the norm operator');
            end
        elseif nargout>=1
            error('Regularization function handle is needed');
        end
    elseif strcmpi(pname,'rho')     % Over-relaxation parameter
        rho=pval;
    elseif any(strcmpi(pname,{'initx','init'}))   % Initialization of x
        initx=pval;
    elseif strcmpi(pname,'verbose') % Verbose
        verbose=pval;
    elseif any(strcmpi(pname,{'tol','tolerance','deltacost','deltacost_max'})) % Tolerance for the cost function precision
        deltacost_max_init=pval;
    elseif strcmpi(pname,'tau') % Tau parameter (Please avoid changing, unless the operator norm of your direct transform is too hard to evaluate)
        tau_init=pval;
    end
end

Nbiter=max(iter_save);
if any(lambda<=0), error('lambda parameter must be positive'); end
[lambda,idx_lambda]=sort(lambda,'descend');
Nlambda=length(lambda);

%% Generating linear operators

if isempty(opA)
    if verbose==1, fprintf('Direct transfer model set to identity.\n'); end
    opA=@(x) x;
	opAadj=@(x) x;
    opAnorm=1;
end
if isempty(opL)
    if verbose==1, fprintf('Regularization set to Total Variation.\n'); end
    opL = @(x) cat(4,[diff(x,1,1);zeros(1,size(x,2),size(x,3))],[diff(x,1,2) zeros(size(x,1),1,size(x,3))]);
    opLadj = @(u) -[u(1,:,:,1);diff(u(:,:,:,1),1,1)]-[u(:,1,:,2) diff(u(:,:,:,2),1,2)];
    opLnorm=sqrt(8);
end
if isempty(prox_gamma_f_conj)       
    if verbose==1, fprintf('Data fidelity set to Quadratic.\n'); end
    prox_gamma_f_conj = @(x) 1/2*x; % proximal operator for the conjugate of the auadratic function
    norm_handle_f = @(x) sum(x(:).^2)/2;
end
if isempty(prox_gamma_g_conj)     
    if verbose==1, fprintf('Norm set to l_221.\n'); end
    opy = @(x,lambda) x ./ max(sqrt(sum(x.^2,2:3)),lambda); % proximal operator for the conjugate of the l221 norm
    prox_gamma_g_conj= @(x,~,lambda) reshape(opy(reshape(x,[],size(x,3),size(x,4)),lambda),size(x));
    norm_handle_g = @(x,lambda) lambda*sum(sqrt(sum(reshape(x,[],size(x,3),size(x,4)).^2,2:3)),1);
end

if ~isempty(opW)
    opA=@(x) opA(opW(x));
    opAadj=@(x) opWadj(opAadj(x));
    opAnorm=opAnorm*opWnorm;
else
    opW    = @(x) x;
    opWinv = @(x) x;
end

if isempty(c)
    if isempty(y), error('The acquisition is needed to run the algorithm');
    else, c=-opAadj(y);
    end
end
    

%% Convergence Parameters
if rho>2 || rho<=0
    fprintf('The assigned value for the over-relaxation parameter does not assure convergence.\n');
end
if isempty(tau_init), tau=0.99/sqrt(opLnorm^2+opAnorm^2);
else, tau=tau_init/sqrt(opLnorm^2+opAnorm^2);
end
sigma=0.99/tau/opLnorm^2;
% eta  =0.99/tau/opAnorm^2;

%% Initialization
if isempty(initx), x=-c; else, x=opWinv(initx); end
x_out=[];
sizex=size(x);
cond_stop=~isempty(norm_handle_g) && ~isempty(norm_handle_f);
if nargin>=1, cost_out=cell(Nlambda,1); end

% sigmatau=sigma*tau;
% etatau=eta*tau;

for kk=1:Nlambda
    
    if ~isempty(initx), x=opWinv(initx); end % Bypass warm start
    
    % Itnitialization (Iteration 0)
    xtilde=x;
    u=opL(x);
    
    if cond_stop
        cost_old=norm_handle_f(opA(x))+sum(x(:).*c(:))+norm_handle_g(opL(x),lambda(kk));
    end
    
    iter=0;
    Ncor=0;
    cost=[];
    
    %% Iterative algorithm
    if verbose==1, fprintf('Current iteration (lambda=%.3e): %5d',lambda(kk),0); end
    while iter < Nbiter && Ncor < Npos
        
        opL=@(x) TVs_direct(x,1);
        opLadj=@(x) -TVs_adjoint(x,1);
        
        % Classical generalized Chambolle-Pock iteration
        u=prox_gamma_g_conj(u+sigma*opL(xtilde),sigma,lambda(kk));
        xtilde=x;
        x=prox_gamma_f_conj(x-tau*(opLadj(u)+c),tau);
        xtilde=2*x-xtilde;
        
        % sigma1=sigma/lambda(kk);
        % tau1=tau/lambda(kk);
        % u=prox_gamma_g_conj(u+sigma1*lambda(kk)*opL(xtilde),sigma1,1);
        % xtilde=x;
        % x=prox_gamma_f_conj(x-tau1*(lambda(kk)*opLadj(u)+c),tau1);
        % xtilde=2*x-xtilde;
        
        iter=iter+1;
        if cond_stop
            %x=xtilde*tau;
            cost_new=norm_handle_f(opA(x))+sum(x(:).*c(:))+norm_handle_g(opL(x),lambda(kk));
            if nargin>=1, cost=cat(1,cost,cost_new); end
            deltacost=abs(cost_old-cost_new);
            if iter==1, deltacost_max=deltacost_max_init*deltacost; end
            if deltacost<deltacost_max, Ncor=Ncor+1; end
            cost_old=cost_new;
        end
        if verbose==1, fprintf([repmat('\b',[1,5]),'%5d'],iter); end
    end
    if verbose==1, fprintf([repmat('\b',[1,5]),'%5d/ %5d. Done!\n'],iter,iter); end
    
    temp=opW(x);
    x_out=cat(ndims(temp)+1,x_out,temp);
    if nargin>=1, cost_out{kk}=cost; end
    
end

x_out=reshape(x_out,[],Nlambda);
x_out=reshape(x_out(:,idx_lambda),[sizex,Nlambda]);