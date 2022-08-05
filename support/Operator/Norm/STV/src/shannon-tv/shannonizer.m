function [u,E] = shannonizer(u0,gam,n,lambda,niter,varargin)
%
%                        Image Shannonizer
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
% Description : given an input image 'u0' with size nx*ny, a frequency
% weigthing function 'gam' defined in the Fourier domain, use the
% Chambolle-Pock Algorithm to compute a minimizer of
%
%         E(u) := ||DFT(u)-DFT(u0)||_gam^2 + lambda*STVn(u)
%
% where DFT denotes the discrete Fourier transform and
%
%  ||DFT(u)-DFT(u0)||_gam^2 = 1/(nx*ny) * sum_{a,b} gam(a,b) * ||DFT(u)(a,b) - DFT(u0)(a,b) ||_2^2
%
% is a weighted l2 square distance between DFT(u) and DFT(u0).
%
%
% Input description:
% ------------------
%
% + u0 : input (noisy) image
% + gam : weighting map in the Fourier domain
% + n : oversampling factor (use STVn regularizer)
% + lambda : regularity parameter (STVn weight)
% + niter : number of iterations for the Chambolle-Pock Algorithm
% + init (default = u0): initialize algorithm with u = init
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
% v1.0 (12/2016): initial version
%
% Versions:
% ---------
%
%   + (2017/01/13) v 1.0 : initial version (Rémy Abergel)
%

%% Initialization
if(nargin < 5)
    fprintf('\nUsage: [u,E] = shannonizer(u0,gam,n,lambda,niter[,''init'',value,''gain'',value,''verbose'',value,''video'',value,''displayFcn'',value])\n\n');
    fprintf('  + u0 : input image\n');
    fprintf('  + gam : weighting map in the Fourier domain\n');
    fprintf('  + n : oversampling factor (use STVn regularizer)\n');
    fprintf('  + lambda : regularity parameter (STVn weight)\n');
    fprintf('  + niter : number of iterations for the Chambolle-Pock Algorithm\n');
    fprintf('  + init (optional, default = zeros(size(adj_A(u0)))): initialize algorithm with u = init\n');
    fprintf('  + gain (optional, default = 1) : multiply default stepsize parameters (tau and sigma) by gain and 1/gain respectively\n');
    fprintf('  + verbose (optional, default = false) : display energy at each iteration\n');
    fprintf('  + video (optional, default = false) : imshow ''displayFcn(u)'' at each iteration\n');
    fprintf('  + displayFcn : (optional, default = @(u)u) : display function for the video mode.\n\n');
    fprintf('  + u : ouptut image\n');
    fprintf('  + E : computed energy (E(u)) at each iteration.\n\n');
    fprintf('Description: minimize E(u) := ||DFT(u)-DFT(u0)||_gam^2 + lambda*STVn(u) using the Chambolle-Pock Algorithm.\n\n');
    fprintf('Reference: R. Abergel and L. Moisan, ``The Shannon Total Variation'''', Journal of Mathematical Imaging and Vision, 2017.\n\n');
    error('Incorrect number of input arguments.');
end

%% Input parser
p = inputParser;
p.addRequired('u0',@(u0) (~isempty(u0)) && isreal(u0));
p.addRequired('gam',@(gam) prod(size(gam)==size(u0)) && isreal(gam));
p.addRequired('n',@(n) double(n)==floor(n) && n>0);
p.addRequired('lambda',@(lambda) lambda>0);
p.addRequired('niter',@(niter) niter>=0);
p.addParameter('init',u0,@(init) prod(size(init)==size(u0)) && isreal(init));
p.addParameter('gain',1,@(gain) gain > 0);
p.addParameter('verbose',false,@(verbose) islogical(verbose));
p.addParameter('video',false,@(video) islogical(video));
p.addParameter('displayFcn',@(u)u,@(A)isa(A,'function_handle'));
parse(p,u0,gam,n,lambda,niter,varargin{:});
init = p.Results.init;
verbose = p.Results.verbose;
gain = p.Results.gain;
video = p.Results.video;
displayFcn = p.Results.displayFcn;

%% Initialization

% deal with video mode
function callback_keyboard_pressed(fig_hdl,evt) % keyboard event handler function
    if strcmp(evt.Key,'q') % stop process & quit
        close(fig_hdl);
    end
end
if(video);
    fig_hdl = figure('Name',sprintf('iteration 0.'),'KeyPressFcn',@callback_keyboard_pressed);
    imshow(displayFcn(init),[],'Border','tight');
    im_hdl = fig_hdl.Children.Children;
end

% primal & dual variables initialization
[ny,nx] = size(u0);
u = init;
ubar = init;
px = zeros(n*ny,n*nx);
py = zeros(n*ny,n*nx);
dft_u0 = fft2(u0);
needE = verbose || (nargout == 2);
if(verbose); E = zeros(1,niter); else E = []; end

%% set primal and dual steps
L = lambda*pi*sqrt(2)/n;
coef = gain * (norm(u(:),2)/sqrt(nx*ny))/128;
tau = coef * 0.99/L;
sigma = (1/coef) * 0.99/L;
lambda_n = lambda/(n^2);

%% Main loop
for k = 1:niter

    % update (px,py) from ubar //
    [gx,gy] = sgrad(ubar,n);
    px = px + sigma*lambda_n * gx;
    py = py + sigma*lambda_n * gy;
    nrm = px.^2 + py.^2;
    id = find(nrm > 1);
    nrm = sqrt(nrm(id));
    px(id) = px(id)./nrm;
    py(id) = py(id)./nrm;

    % update u and ubar from (px,py)
    ubar = -u;
    u = u + tau*lambda_n*sdiv(px,py,n);
    dft_u = fft2(u);
    dft_u = (dft_u + 2*tau*(gam.*dft_u0))./(1+2*tau*gam);
    u = real(ifft2(dft_u));
    ubar = ubar + 2*u;

    % deal with verbose mode
    if(verbose)
        [gy,gx] = sgrad(u,n);
        E(k) = sum(sum(gam.*abs(dft_u-dft_u0).^2))/(nx*ny) + lambda_n*sum(sum(sqrt(gx.^2+gy.^2)));
        fprintf('  + iteration %d : E = %.10g\n',k,E(k));
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
