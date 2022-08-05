function [u,E] = tvdenoise(u0,lambda,niter,varargin)
%
%            Discrete Total Variation based image denoising
%
% Author: Rémy Abergel.
%
% This program is freely available on the web page
%
% http://www.math-info.univ-paris5.fr/~rabergel/
%
% I hope that you will find it useful. If you use it for a publication,
% please mention this web page. If you modify this program, please indicate
% in the code that you did so and  leave this message. You can also report
% bugs or suggestions to Remy.Abergel [AT] gmail.com
%
% Description : compute a minimizer of 
%
%                 E(u) := ||u-u0||^2 + lambda*TV(u)
%
% using the Chambolle-Pock Algorithm. 
%
%
% Input description: 
% ------------------
%
% + u0 : input (noisy) image
% + lambda : regularity parameter (TV weight)
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
% v1.0 (12/2016): initial version
%
% Versions: 
% ---------
%
%   + (2017/01/13) v 1.0 : initial version (Rémy Abergel)
%
%% Usage displayer
if(nargin < 3) 
    fprintf('\nUsage: [u,E] = tvdenoise(u0,lambda,niter[,''init'',value,''gain'',value,''verbose'',value,''video'',value,''displayFcn'',value])\n\n'); 
    fprintf('  + u0 : input (noisy) image\n');
    fprintf('  + lambda : regularity parameter (TV weight)\n');
    fprintf('  + niter : number of iterations for the Chambolle-Pock Algorithm\n'); 
    fprintf('  + init (optional, default = u0): initialize algorithm with u = init\n');      
    fprintf('  + gain (optional, default = 1) : multiply default stepsize parameters (tau and sigma) by gain and 1/gain respectively\n');      
    fprintf('  + verbose (optional, default = false) : display energy at each iteration\n'); 
    fprintf('  + video (optional, default = false) : imshow ''displayFcn(u)'' at each iteration\n'); 
    fprintf('  + displayFcn : (optional, default = @(u)u) : display function for the video mode\n\n'); 
    fprintf('  + u : ouptut image\n');
    fprintf('  + E : computed energy (E(u)) at each iteration\n\n');    
    fprintf('Description: minimize E(u) := ||u-u0||_2^2 + lambda*TV(u) using the Chambolle-Pock Algorithm.\n\n'); 
    error('Incorrect number of input arguments.');     
end

%% Input parser
p = inputParser;
p.addRequired('u0',@(u0) (~isempty(u0)) && isreal(u0));
p.addRequired('lambda',@(lambda) lambda>0);
p.addRequired('niter',@(niter) niter>=0);
p.addParameter('init',u0,@(init) prod(size(init)==size(u0)) && isreal(init));
p.addParameter('gain',1,@(gain) gain > 0);
p.addParameter('verbose',false,@(verbose) islogical(verbose));
p.addParameter('video',false,@(video) islogical(video));
p.addParameter('displayFcn',@(u)u,@(A)isa(A,'function_handle'));
parse(p,u0,lambda,niter,varargin{:});
init = p.Results.init; 
verbose = p.Results.verbose; 
gain = p.Results.gain; 
video = p.Results.video; 
displayFcn = p.Results.displayFcn; 


%% Initialization 
[ny,nx] = size(u0); % size of Omega := domain of u 
u = init;
ubar = u;
px = zeros(ny,nx);
py = zeros(ny,nx);
needE = verbose || (nargout == 2); 
if(needE); E = zeros(1,niter); else E = []; end

% set primal & dual steps for Chambolle-Pock Algorithm
L = lambda*sqrt(8); 
coef = gain * (norm(u(:),2)/sqrt(nx*ny))/128; 
tau = coef * 0.99/L;
sigma = (1/coef) * 0.99/L;

% deal with video mode
if(video); 
    fig_hdl = figure('Name',sprintf('iteration 0.')); 
    imshow2(displayFcn(u)); 
    im_hdl = fig_hdl.Children.Children; 
end

%% main loop 
for k = 1:niter

    % update (px,py) from ubar 
    [gx,gy] = grad(ubar);
    px = px + sigma*lambda * gx;
    py = py + sigma*lambda * gy;
    nrm = px.^2 + py.^2;
    id = find(nrm > 1);
    nrm = sqrt(nrm(id));
    px(id) = px(id)./nrm;
    py(id) = py(id)./nrm;

    % update u and ubar from (px,py)
    ubar = -u;
    u = (u + tau*lambda*div(px,py) + 2*tau*u0) / (1+2*tau);
    ubar = ubar + 2*u;

    % deal with verbose mode
    if(needE)
        [gy,gx] = grad(u);
        E(k) = sum(sum((u-u0).^2 + lambda*sqrt(gx.^2+gy.^2)));
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

