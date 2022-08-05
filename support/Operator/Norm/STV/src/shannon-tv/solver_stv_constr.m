function [u,E] = solver_stv_constr(nx,ny,n,proj_constr,niter,varargin)
%
%        Minimize the Shannon Total Variation under constraint
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
% Use the Chambolle-Pock Algorithm to compute
%
%                  argmin_{u in C} STVn(u)
%
% where the constraint set C is a subset of R^{nx*ny}.
%
% Input description:
% ------------------
%
% + nx : output width
% + ny : output height
% + n : oversampling factor (use STVn regularizer)
% + proj_constr (Matlab function or macro): projector on the constraint set C
% + niter : number of iterations for the Chambolle-Pock Algorithm
% + init (default = proj_constr(zeros(ny,nx))): initialize algorithm with u = init
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
if(nargin < 5)
    fprintf('\nUsage: [u,E] = solver_stv_constr(nx,ny,proj_constr,niter[,''init'',value,''gain'',value,''verbose'',value,''video'',value,''displayFcn'',value])\n\n');
    fprintf('  + nx : output width\n');
    fprintf('  + ny : output height\n');
    fprintf('  + n : oversampling factor (use STVn regularizer)\n');
    fprintf('  + proj_constr\n');
    fprintf('  + proj_constr (Matlab function or macro) : projector on the constraint set C\n');
    fprintf('  + niter : number of iterations for the Chambolle-Pock Algorithm\n');
    fprintf('  + init (optional, default = zeros(size(adj_A(u0)))): initialize algorithm with u = init\n');
    fprintf('  + gain (optional, default = 1) : multiply default stepsize parameters (tau and sigma) by gain and 1/gain respectively\n');
    fprintf('  + verbose (optional, default = false) : display energy at each iteration\n');
    fprintf('  + video (optional, default = false) : imshow ''displayFcn(u)'' at each iteration\n');
    fprintf('  + displayFcn : (optional, default = @(u)u) : display function for the video mode.\n\n');
    fprintf('  + u : ouptut image\n');
    fprintf('  + E : computed energy (E(u) := STVn(u)) at each iteration\n\n');
    fprintf('Description: Usethe Chambolle-Pock Algorithm to minimize E(u) := STVn(u) over the constraint set C.\n\n');
    fprintf('Reference: R. Abergel and L. Moisan, ``The Shannon Total Variation'''', Journal of Mathematical Imaging and Vision, 2017.\n\n');
    error('Incorrect number of input arguments.');
end

%% Input parser
p = inputParser;
p.addRequired('nx',@(n) length(n) == 1 && double(n)==floor(n) && n>0);
p.addRequired('ny',@(n) length(n) == 1 && double(n)==floor(n) && n>0);
p.addRequired('n',@(n) length(n) == 1 && double(n)==floor(n) && n>0);
p.addRequired('proj_constr',@(proj_constr)isa(proj_constr,'function_handle') && prod(size(proj_constr(zeros(ny,nx)))==[ny,nx]));
p.addRequired('niter',@(niter) niter>=0);
p.addParameter('init',proj_constr(zeros(ny,nx)),@(init) prod(size(init)==size(adj_A(u0))) && isreal(init));
p.addParameter('gain',1,@(gain) gain > 0);
p.addParameter('verbose',false,@(verbose) islogical(verbose));
p.addParameter('video',false,@(video) islogical(video));
p.addParameter('displayFcn',@(u)u,@(A)isa(A,'function_handle'));
p.parse(nx,ny,n,proj_constr,niter,varargin{:});
init = p.Results.init;
verbose = p.Results.verbose;
gain = p.Results.gain;
video = p.Results.video;
displayFcn = p.Results.displayFcn;

%% Initialization

% deal with video mode

if(video);
    fig_hdl = figure('Name',sprintf('iteration 0.'));
    imshow2(displayFcn(init));
    im_hdl = fig_hdl.Children.Children;
end

% primal & dual variables initialization
u = init; ubar = u; % primal variables
px = zeros(n*ny,n*nx); py = zeros(n*ny,n*nx);  % dual variables
needE = verbose || (nargout == 2);
if(needE); E = zeros(niter,1); else E = []; end

% set primal & dual steps for Chambolle-Pock Algorithm
L = sqrt(2)*pi*n;
coef = gain * (norm(init(:),2)/sqrt(nx*ny))/128;
tau = coef * 0.99/L;
sigma = (1/coef) * 0.99/L;
sigma_n = sigma/(n^2);
tau_n = tau/(n^2);

%% Main loop (Chambolle-Pock Algorithm)
for k = 1:niter

    % compute (gx,gy) = discrete gradient of ubar
    [gx,gy] = sgrad(ubar,n);

    % update (px,py) from ubar
    px = px + sigma_n * gx;
    py = py + sigma_n * gy;
    nrm = px.^2 + py.^2;
    id = find(nrm > 1);
    nrm = sqrt(nrm(id));
    px(id) = px(id)./nrm;
    py(id) = py(id)./nrm;

    % update u and ubar from (px,py) and q
    ubar = -u;
    u = proj_constr(u + tau_n*sdiv(px,py,n));
    ubar = ubar + 2*u;

    % deal with verbose mode
    if(needE)
        [gy,gx] = sgrad(u,n);
        E(k) = sum(sum(sqrt(gx.^2+gy.^2)))/(n^2);
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
