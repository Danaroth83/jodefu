function [u,E] = solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,varargin)
%
%       Solver for inverse problems with quadratic data-fidelity
%              and Shannon total variation regularity
%
% Authors: Rémy Abergel & Lionel Moisan.
%
% This program is freely available on the web page
%
% http://www.math-info.univ-paris5.fr/~rabergel/
%
% We hope that you will find it useful. If you use it for a publication,
% please mention this web page and the paper [1]. If you modify this
% program, please indicate in the code that you did so and  leave this
% message. You can also report bugs or suggestions to Remy.Abergel [AT]
% gmail.com and Lionel.Moisan [AT] parisdescartes.fr
%
% [1] R. Abergel and L. Moisan, ``The Shannon Total Variation'', Preprint MAP5,
% 2016.
%
%
% compute a minimizer of
%
%                 E(u) := ||A(u)-u0||^2 + lambda*STVn(u)
%
% using the Chambolle-Pock Algorithm.
%
%
% Input description:
% ------------------
%
% + u0 : input image
% + n : oversampling factor (use STVn regularizer)
% + lambda : regularity parameter (STVn weight)
% + A (Matlab function or macro with prototype v = A(u)) : linear operator to be inverted
% + adj_A (Matlab function or macro with prototype u = adj_A(v)) : adjoint of A
% + L_A : an upper bound for |||A||| (the the induced l2 norm of the operator A)
% + niter : number of iterations for the Chambolle-Pock Algorithm
% + init (default = zeros(size(u))): initialize algorithm with u = init
% + gain (default = 1): multiply default stepsize parameters (tau and sigma) by gain and 1/gain respectively
% + verbose (boolean, default = false): display energy at each iteration
% + video (boolean, default = false): imshow 'displayFcn(u)' at each iteration
% + displayFcn (Matlab function or macro, default = @(u)u): display function for the video mode.
%
% Output description:
% -------------------
%
% + u : restored image
% + E : energy (E(k) = energy computed at iteration k)
%
%
% Versions:
% ---------
%
%   + (2017/01/13) v 1.0 : initial version (Rémy Abergel)
%

%% Usage displayer
if(nargin < 7)
    fprintf('\nUsage: [u,E] = solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter[,''init'',value,''gain'',value,''verbose'',value,''video'',value,''displayFcn'',value])\n\n');
    fprintf('  + u0 : input image\n');
    fprintf('  + n : oversampling factor (use STVn regularizer)\n');
    fprintf('  + lambda : regularity parameter (STVn weight)\n');
    fprintf('  + A (Matlab function or macro with prototype v = A(u)) : linear operator to be inverted\n');
    fprintf('  + adj_A (Matlab function or macro with prototype u = adj_A(v)) : adjoint of A\n');
    fprintf('  + L_A : upper bound for the induced l2 norm of A\n');
    fprintf('  + niter : number of iterations for the Chambolle-Pock Algorithm\n');
    fprintf('  + init (optional, default = zeros(size(adj_A(u0)))): initialize algorithm with u = init\n');
    fprintf('  + gain (optional, default = 1) : multiply default stepsize parameters (tau and sigma) by gain and 1/gain respectively\n');
    fprintf('  + verbose (optional, default = false) : display energy at each iteration\n');
    fprintf('  + video (optional, default = false) : imshow ''displayFcn(u)'' at each iteration\n');
    fprintf('  + displayFcn : (optional, default = @(u)u) : display function for the video mode.\n\n');
    fprintf('  + u : ouptut image\n');
    fprintf('  + E : computed energy (E(u)) at each iteration.\n\n');
    fprintf('Description: minimize E(u) := ||A(u)-u0||_2^2 + lambda*STVn(u) using the Chambolle-Pock Algorithm.\n\n');
    fprintf('Reference: R. Abergel and L. Moisan, ``The Shannon Total Variation'''', Journal of Mathematical Imaging and Vision, 2017.\n\n');
    error('Incorrect number of input arguments.');
end

%% Input parser
p = inputParser;
p.addRequired('u0',@(u0) (~isempty(u0)) && isreal(u0));
p.addRequired('n',@(n) double(n)==floor(n) && n>0);
p.addRequired('lambda',@(lambda) lambda>0);
p.addRequired('A',@(A)isa(A,'function_handle'));
p.addRequired('adj_A',@(adj_A)isa(adj_A,'function_handle'));
p.addRequired('L_A',@(L_A) L_A>0);
p.addRequired('niter',@(niter) niter>=0);
p.addParameter('init',zeros(size(adj_A(u0))),@(init) prod(size(init)==size(adj_A(u0))) && isreal(init));
p.addParameter('gain',1,@(gain) gain > 0);
p.addParameter('verbose',false,@(verbose) islogical(verbose));
p.addParameter('video',false,@(video) islogical(video));
p.addParameter('displayFcn',@(u)u,@(A)isa(A,'function_handle'));
parse(p,u0,n,lambda,A,adj_A,L_A,niter,varargin{:});
init = p.Results.init;
verbose = p.Results.verbose;
gain = p.Results.gain;
video = p.Results.video;
displayFcn = p.Results.displayFcn;

%% Initialization

% deal with video mode
if(video)
    fig_hdl = figure('Name',sprintf('iteration 0.'));
    imshow2(displayFcn(init));
    im_hdl = fig_hdl.Children.Children;
end

% primal & dual variables initialization
[ny0,nx0] = size(u0); % size of omega := domain of u0 = domain of A(u)
[ny,nx] = size(init); % size of Omega := domain of u = domain of adj_A(u0)
Nx = n*nx; Ny = n*ny; % size of Omega_n := domain of the dual variable (px,py)
u = init; ubar = u; % primal variables
px = zeros(Ny,Nx); py = zeros(Ny,Nx);  % p = (px,py) : first dual variable
q = zeros(ny0,nx0); % second dual variable
needE = verbose || (nargout == 2);
if(needE); E = zeros(niter,1); else E = []; end

% set primal & dual steps for Chambolle-Pock Algorithm
L = sqrt(2*(pi*lambda/n)^2+L_A^2);
coef = gain * (norm(u0(:),2)/sqrt(nx0*ny0))/128;
tau = coef * 0.99/L;
sigma = (1/coef) * 0.99/L;
lambda_n = lambda/(n^2);

%% Main loop (Chambolle-Pock Algorithm)
for k = 1:niter

    % compute (gx,gy) = discrete gradient of ubar
    [gx,gy] = sgrad(ubar,n);

    % update (px,py) from ubar
    px = px + sigma*lambda_n * gx;
    py = py + sigma*lambda_n * gy;
    nrm = px.^2 + py.^2;
    id = find(nrm > 1);
    nrm = sqrt(nrm(id));
    px(id) = px(id)./nrm;
    py(id) = py(id)./nrm;

    % update q from ubar
    q = (2*q + 2*sigma*(A(ubar)-u0))/(2+sigma);

    % update u and ubar from (px,py) and q
    ubar = -u;
    u = u + tau*lambda_n*sdiv(px,py,n) - tau*adj_A(q);
    ubar = ubar + 2*u;

    % deal with verbose mode
    if(needE)
        [gy,gx] = sgrad(u,n);
        E(k) = sum(sum((A(u)-u0).^2)) + lambda_n*sum(sum(sqrt(gx.^2+gy.^2)));
        if(verbose); fprintf('  + iteration %d : E = %.16g\n',k,E(k)); end
    end

    % deal with video mode
    if(video && isvalid(fig_hdl))
        if(needE)
            fig_hdl.Name = sprintf('iteration %d: E = %.10g',k,E(k));
        else
            fig_hdl.Name = sprintf('iteration %d.',k);
        end
        im_hdl.CData = displayFcn(u);
        drawnow;
    end

    % deal with stop event (video mode only)
    if (video && ~isvalid(fig_hdl))
        if(needE); E = E(1:k); end
        break; % leave for loop & return
    end

end

end
