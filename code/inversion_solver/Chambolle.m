%% GENERALIZED CHAMBOLLE-POCK INVERSION SOLVER
%
% (copyright) Daniele Picone
% version 2.0
% 
% Description:
% Chambolle is a solver the problem
% argmin_x  f(Ax) + g(Lx) + <x,x0>
% where: A and L are bounded linear operators
%        f(x) and g(x) are real valued lower semicontinuous convex function
%        c is a parameter in the dual space of x
% It implements  the Chambolle-Pock algorithm (110) in [1].
% 
% Usage:
% x=Chambolle('Field',Value);
%
% Example:
% im=double(imread('peppers.png'))/255;
% y=im+0.1*randn(size(im));
% opA.direct= @(x) x;
% opA.adjoint= @(x) x;
% opA.norm=1;
% opL.direct = @(x) cat(4,[diff(x,1,1);zeros(1,size(x,2),size(x,3))],[diff(x,1,2) zeros(size(x,1),1,size(x,3))]);
% opL.adjoint = @(u) -[u(1,:,:,1);diff(u(:,:,:,1),1,1)]-[u(:,1,:,2) diff(u(:,:,:,2),1,2)];
% opL.norm=sqrt(8);
% opy = @(x,lambda) x ./ max(sqrt(sum(x.^2,2:3)),lambda);
% opg.proxconj= @(x,~,lambda) reshape(opy(reshape(x,[],size(x,3),size(x,4)),lambda),size(x));
% opg.direct = @(x,lambda) lambda*sum(sqrt(sum(reshape(x,[],size(x,3),size(x,4)).^2,2:3)),1);
% opf.proxconj= @(x,~) x/2; 
% opf.direct = @(x) sum(x.^2,1:3)/2;
% 
% [x,cost]=Chambolle('c',-opA.adjoint(y),'Nbiter',250,'opA',opA,'opL',opL,...
%                'g',opg,'f',opf,'lambda',0.14:0.01:0.15,'rho',1,'tol',10E-6);
%
% Input:
% y: observation variable
%
% Fields:
% 'lambda': Array of regularization parameter (Array of scalars >0 - default: 0.1)
% 'opA': struct for calculating the Direct Tranfer Function Ax, whose fields are:
%      - direct:  Handle for the direct tranfer model Ax
%      - adjoint: Handle for the adjoint operator of the above
%      - norm:    Operator norm ||A||
% 'opL': struct for calculating Lx, whose fields are:
%      - direct:  Handle for calculating Lx
%      - adjoint: Handle for the adjoint operator of the above
%      - norm:    Operator norm ||L||
%
% 'g': struct relative to the operators relative to the function g(x), whose fields are:
%      - proxconj: handle for the proximal operator of \gamma multiplied
%                  by the convex conjugate of the function g(x); it has to
%                  have two inputs (x,gamma). It is allowed to have three
%                  inputs, (x,gamma,lambda) in case g(x) allows a
%                  multiplicative parameter
%      - prox:     handle for the proximal operator of \gamma*g(x); it has 
%                  two inputs (x,\gamma) and generates the handle proxconj
%                  via Moreau identity
%      - direct: handle to calculate g(x)
% 'f': struct relative to the operators relative to the function f(x), whose fields are:
%      - proxconj: handle for the proximal operator of \gamma multiplied
%                  by the convex conjugate of the function g(x); it has to
%                  have two inputs (x,gamma). It is allowed to have three
%                  inputs, (x,gamma,lambda) in case g(x) allows a
%                  multiplicative parameter
%      - prox:     handle for the proximal operator of \gamma*g(x); it has 
%                  two inputs (x,\gamma) and generates the handle proxconj
%                  via Moreau identity
%      - direct: handle to calculate g(x)
% 'Nbiter': Max number of iterations (default: 1000).
% 'init': Initialization for the solution x (Default: adjoint transformation applied on y)
% 'rho': Over-relaxation parameter (Must be 1<=rho<2, default: 1)
% 'tol': difference in cost value to stop iterations, scaled to the
%              norm of y (If 0, check is ignored. Default: 1E-7)
% 'verbose': If 1, prints the iteration status (Default: 1)
% 
% References:
% [1] Condat L.,  Kitahara D., Contreras A. and Hirabayashi A., "Proximal
% splitting algorithms: Relax them all!" (preprint)

function [x_out,cost_out]=Chambolle(varargin)

%% Parse inputs

if rem(numel(varargin),2)~=0
    fprintf('Usage: [x,cost]=Chambolle(''Nbiter'',Nbiter,''opA'',opA,''opL'',opL,''f'',proxf,''g'',''proxg'')\n');
    error('Wrong input format');
end

iter_save=1000;
lambda=0.1;
verbose=1;
opA=[];
opL=[];
opW=[];
rho=1;
prox_gamma_g_conj=[];
prox_gamma_f_conj=[];
initx=[];
tau_init=[];
deltacost_max=10E-7;
Npos=3; % Number of consecutive correct checks to stop the iterations

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if any(strcmpi(pname,{'Nbiter','iter_save'})) % Number of iterations
        iter_save=pval;
    elseif strcmpi(pname,'lambda')  % Regularization Parameter
        lambda=pval;   
    elseif strcmpi(pname,'opA')     % Handle of operator A
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
    elseif strcmpi(pname,'opL')     % Handle of operator L
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
    elseif strcmpi(pname,'opW')     % Handle of operator W
        if isnumeric(pval), opW=pval; end
        if isfield(pval,'direct')
            opW=pval.direct;
            if isfield(pval,'adjoint'), opWadj=pval.adjoint; else, error('Missing information on adjoint operator'); end
            if isfield(pval,'norm'), opWnorm=pval.norm; else, error('Missing information on operator norm'); end
        end
        if isnumeric(opW)
            opWnorm=max(svd(opW));
            opWadj=@(x) opW'*x;
            opWinv=@(x) opW\x;
            opW=@(x) opW*x; 
        end
    elseif any(strcmpi(pname,{'f','opf','prox','proximal'})) % Proximal operator for function f
        if isfield(pval,'proxconj')
            if nargin(pval.proxconj)==3, prox_gamma_f_conj=@(x,gamma) pval.proxconj(x,gamma,1);
            elseif nargin(pval.proxconj)==2, prox_gamma_f_conj=pval.proxconj;
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
    elseif any(strcmpi(pname,{'opg','g','regularizer','reg','norm'})) % Proximal operator for function g
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
    elseif any(strcmpi(pname,{'c','scalar'}))   % Initialization of c
        c=pval;    
    elseif strcmpi(pname,'rho')     % Over-relaxation parameter
        rho=pval;
    elseif any(strcmpi(pname,{'initx','init'}))   % Initialization of x
        initx=pval;
    elseif strcmpi(pname,'verbose') % Verbose
        verbose=pval;
    elseif any(strcmpi(pname,{'tol','tolerance','deltacost','deltacost_max'})) % Tolerance for the cost function precision
        deltacost_max=pval;
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
if isempty(prox_gamma_g_conj)
    if isempty(prox_gamma_g)       
        if verbose==1, fprintf('Norm set to l_221.\n'); end
        opy = @(x,lambda) x ./ max(sqrt(sum(x.^2,2:3)),lambda); % proximal operator for the conjugate of the l221 norm
        prox_gamma_g_conj= @(x,~,lambda) reshape(opy(reshape(x,[],size(x,3),size(x,4)),lambda),size(x));
        norm_handle_g = @(x,lambda) lambda*sum(sqrt(sum(reshape(x,[],size(x,3),size(x,4)).^2,2:3)),1);
    end
end
if isempty(prox_gamma_f_conj)
    if isempty(prox_gamma_f)       
        if verbose==1, fprintf('Main function set to half square l2 norm.\n'); end
        prox_gamma_f_conj= @(x,gamma) x/(1-gamma);
        norm_handle_f = @(x) sum(x.^2,1:3)/2;
    end
end

if ~isempty(opW)
    opA=@(x) opA(opW(x));
    opAadj=@(x) opWadj(opAadj(x));
    opAnorm=opAnorm*opWnorm;
else
    opW=@(x) x;
    opWinv=@(x) x;
end


%% Initialization
if isempty(initx) && isempty(c), error('The algorithm needs to be initialized'); end
if ~isempty(initx), x=opWinv(initx); else, x=-c; end

x_out=[];
sizex=size(x);
cond_stop=~isempty(norm_handle_g);
if nargin>=1, cost_out=cell(Nlambda,1); end
if isempty(tau_init), tau=min(0.99/opAnorm,0.99/opLnorm); else, tau=tau_init; end
% if tau>2/opAnorm^2, warning('Parameter tau does not assure convergence'); end
if rho<0 || rho>2, error('Parameter rho has to be between 0 and 2'); end
eta=0.99/tau/opAnorm^2;
sigma=0.99/tau/opLnorm^2;

etatau=eta*tau;
sigmatau=sigma*tau;

for kk=1:Nlambda
    
    if ~isempty(initx), x=opWinv(initx); end % Bypass warm start
    
    cost_old=norm_handle_f(opA(x))+c(:)'*x(:)+norm_handle_g(opL(x),lambda(kk));
    
    % Starting iteration
    xbar=x/tau;
    v=prox_gamma_f_conj(etatau*opA(xbar),eta);
    l=opAadj(v);
    u=prox_gamma_g_conj(sigmatau*opL(-2*l),sigma,lambda(kk));
    b=l;
    l=opLadj(u)+b+c;
    
    x=xbar*tau;
    cost_new=norm_handle_f(opA(x))+c(:)'*x(:)+norm_handle_g(opL(x),lambda(kk));
    deltacost_max=deltacost_max*abs(cost_old-cost_new); % Normalization of the condition to stop iteration
    cost_old=cost_new;

    iter=1;
    cost=cost_new;
    Ncor=0;
    
    %% Iterative algorithm
    if verbose==1, fprintf('Current iteration (lambda=%.3e): %5d',lambda(kk),0); end
    while iter < Nbiter && Ncor<Npos
        
        % Main iteration
        v1=prox_gamma_f_conj(v+etatau*opA(xbar-l),eta);
        l1=l+opAadj(v1)-b;
        u1=prox_gamma_g_conj(u+sigmatau*opL(xbar-2*l1),sigma,lambda(kk));
        xbar=xbar-rho*l1;
        u=u+rho*(u1-u);
        v=v+rho*(v1-v);
        b=b+rho*(l1-l);
        l=opLadj(u)+b+c;
        
        % Stopping condition        
        if cond_stop
            x=xbar*tau;
            cost_new=norm_handle_f(opA(x))+c(:)'*x(:)+norm_handle_g(opL(x),lambda(kk));
            if nargin>=1, cost=cat(1,cost,cost_new); end
            deltacost=abs(cost_old-cost_new);
            if deltacost<deltacost_max, Ncor=Ncor+1; end
            cost_old=cost_new;
        end
        iter=iter+1;
        if verbose==1, fprintf([repmat('\b',[1,5]),'%5d'],iter); end
    end
    if verbose==1, fprintf([repmat('\b',[1,5]),'%5d/ %5d. Done!\n'],iter,iter); end
    
    temp=opW(xbar*tau);
    x_out=cat(ndims(temp)+1,x_out,temp);
    if nargin>=1, cost_out{kk}=cost; end
    
end

x_out=reshape(x_out,[],Nlambda);
x_out=reshape(x_out(:,idx_lambda),[sizex,Nlambda]);