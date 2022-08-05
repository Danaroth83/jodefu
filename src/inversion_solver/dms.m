function [res] = dms(varargin)
%DMS Solves the following Discrete Mumford-Shah model using the Semi-Linearized 
%    Proximal Alternating Minimization (SL-PAM) for joint image restoration 
%    and segmentation:
%
%      minimize  J(u,e) := L(Au,z) + beta.||(1-e)*Du||^2 + lambda.R(e)
%        u,e
%   
%   where   - z in R^{LM} is the initial degraded image (created from IM)
%           - u in R^{NM} is a piecewise smooth approximation of z
%           - e in R^|E| is equal to 1 when a contour change is detected, 
%                and 0 otherwise.
%           - A in R^{LMxNM} models a blurring operator
%           - D in R^{|E|xN} models a finite difference operator
%           - L(A.,z) : R^{NM} -> (-inf,+inf] is a data fidelity term
%           - R : R^|E| -> (-inf,+inf] favors sparse solutions
%           - beta > 0 is the smoothing parameter
%           - lambda > 0 controls the length of the contours
%
%   RES = DMS(IM, BETA, LAMBDA) returns a structure containing the input
%   data, and the solution of the Discrete Mumford-Shah problem for a given 
%   set of parameters. IM can be 2D grayscale or multivariate image.
%   The fieldnames of RES include:
%       RES.ground_truth : the given image IM
%       RES.data         : the degraded image z
%       RES.param        : the set of parameters
%       RES.u            : the resulting image u
%       RES.e            : the resulting contour e
%       RES.energy       : the value of the energy at the end
%       RES.last_var     : the last variation of the energy |J_[k+1]-J_[k]|
%       RES.time         : the computational time
%
%   RES = DMS(IM, BETA, LAMBDA, NAME1, VAL1,...) solves the Discrete 
%   Mumford-Shah problem using name-value pairs to control aspects of the 
%   computation.
%   
%   Parameters include:
%
%   'AddBlur'             - Three-element vector, [TYPE SZ STD], of two 
%                           non-negative scalars and a non-negative real
%                           number that specifies:
%                               TYPE: 0 for no blur;
%                                     1 for Gaussian blur; 
%                                     2 for Uniform blur.
%
%                               SZ:   the filter size.
%
%                               STD:  the standard deviation of the
%                                      Gaussian blur.
%                           By default, no blur is added, i.e., the vector
%                           is [0 0 0].
% 
%   'AddNoise'            - Two-element vector, [TYPE STD], of a 
%                           non-negative scalar and a non-negative real 
%                           number that specifies:
%                               TYPE: 0 for no noise;
%                                     1 for Gaussian noise; 
%                                     2 for Poisson noise.
%
%                               STD:  the standard deviation of the noise.
%                           By default, no noise is added, i.e., the vector
%                           is [0 0].
% 
%   'InitU'               - String, INIT_U, that specifies how to
%                           initialize u:
%                               'data': equal to ;
%                               'rand': random values drawn from the 
%                                       standard normal distribution;
%                               'rand_cst': identically equal to a random 
%                                           constant in (0,1).
%                           The default value is 'data'.
% 
%   'InitEdges'           - String, INIT_EDGES, that specifies how to
%                           initialize the contour e:
%                               '0': identically equal to 0 (no contour); 
%                               '1': identically equal to 1 (contours 
%                                    everywhere);
%                               'rand_cst': identically equal to a random 
%                                           constant in (0,1);
%                               'binomial': binomial law with probability
%                                           0.5.
%                           The default value is '1'.
%
%   'Regularization'      - String, REG, that specifies the regularization 
%                           type, between 'l0', 'l1', and the quadratic-l1
%                           regularization 'l1q' (see [1]). The default 
%                           value is 'l1q'.
% 
%   'Edges'               - String, EDGES, that specifies if the contours  
%                           is defined as 'similar' edges through all the 
%                           components, or as 'distinct' edges. By default, 
%                           we set 'similar' edges. 
% 
%   'maxNumberIteration'  - Positive scalar, N, that specifies the maximum
%                           number of iterations in the main loop. By 
%                           default, the maximum number of iteration is
%                           10000.
% 
%   'QuadL1Eps'           - Positive real number that specifies the value 
%                           of the epsilon parameter, EPS, used in the 
%                           definition of the quadratic-l1 penalization 
%                           (see [1]):
%                              sum_i  max {|v_i|,|v_i|^2/(4.eps)}
%                           Notice that, since v takes values in [0,1], 
%                           this parameter should be smaller than 1. The 
%                           default value is 0.5.
%
%   Example
%   ---------
%   This example shows how to add Gaussian noise with standard deviation
%   0.1 and to perform joint restoration and segmentation using DMS given
%   the original reference image.
% 
%   im = imread('cameraman.tif');
% 
%   res = dms(im,10,0.02,'AddNoise',[1 0.1]);
% 
%   figure
%   subplot(221); imshow(res.ground_truth); title('Input image');
%   subplot(222); imshow(res.data);         title('Degraded image');
%   subplot(223); imshow(res.u);            title('Restored image');
%   subplot(224); plot_contours(res.e);     title('Contours');
%
%   Class Support
%   -------------
%   Input array IM must be one of the following classes: uint8, int16,
%   uint16, single, or double.
% 
%   References:
%   -----------
%   [1] M. Foare, N. Pustelnik, L. Condat, "Semi-Linearized Proximal 
%       Alternating Minimization for a Discrete Mumford-Shah Model", 
%       IEEE Transactions on Image Processing, 2019.
%
%
%   Version 1.1, Feb. 28, 2020.
%   M. Foare -- marion.foare@ens-lyon.fr

narginchk(3,21);

%% Initialization
[z0,bet,lam,reg,init,blur,noise,edges,N,eps,alp] = parseInputs(varargin{:});

if isempty(z0), return; end

if isa(z0,'int16') % int16 is the only allowed signed-integer type for im.
    % Add offset for signed-integer types to bring values in the
    % non-negative range.
    z0 = (single(z0) - double(intmin('int16')))./255.;
elseif isinteger(z0)
    z0 = single(z0)./255.;
end

% degradation
z = z0;
[z,blur]  = addBlur(z,blur); 
[z,noise] = addNoise(z,noise);

%z = normalize(z);
    
% dimension of e
switch edges
    case 'similar',  dim_e = [2,1];  
    case 'distinct', dim_e = [2,size(z,3)];
end

% save parameters
p = {'bet','lam','alp','reg','init','blur','noise','edges','dim_e','N','eps'};
params = cell2struct({bet,lam,alp,reg,init,blur,noise,edges,dim_e,N,eps},p,2);

% load usual norms and proximity operators
norms = reg_norms(eps);
prox  = proximity_operators(eps);

% load objective functional
[J,L,S,R] = energy(params,norms,prox);
DMS = cell2struct({J,L,S,R},{'J','L','S','R'},2);

%% Main loop 
[un,en] = initialization(z,params);
Jn      = -ones(1,N,'double');

% descent parameters
norm_D     = sqrt(2.);
lipch_S_du = bet*norm_D*norm_D + 1e-16;

ck = 1.01*lipch_S_du;
dk = ck/1000;

% loop parameters
it        = 1;
err       = 1.;
stop_err  = 10^(-4);
Jn(it)    = J(un,en,z);

fprintf(1,'*** bet=%3g \t lam=%1.5g ... ', bet,lam);
display_progress(it,N,'percentage');
tic;
while (it < N) && (err > stop_err)
    % update of u
    un = L.prox(un - (bet/ck)*S.du(un,en),alp/ck,z);
     
    % update of e
    num = bet*S.D_sq(un) + dk/2.*en; 
    den = bet*S.D_sq(un) + dk/2.; 
    en  = R.prox(num./den,lam./(2*den));
    
    % update of the energy
    Jn(it+1) = J(un,en,z);    
          
    % loop parameters
    if mod(it,N/100) == 0, display_progress(it,N,'percentage'); end
    err = abs(Jn(it+1) - Jn(it));
    it  = it + 1;
end
t = toc;
fprintf(1,' done. (last var: %1.5g, %5.3g seconds)*** \n\n',err,t);

% output
iend = find(Jn==-1,1)-1;
if isempty(iend), iend = length(Jn);    end
r = {'ground_truth','data','param','DMS','u','e','energy','last_var','time'};
res = cell2struct({z0,z,params,DMS,un,en,Jn(1:iend),err,t},r,2);

end


function [J,L,S,R] = energy(params,norms,prox)
%% J(u,v) = L(u) + bet*S(u,e) + lam*R(e)
% L: data fidelity
% S: smoothing term, preserving edges
% R: contour length penalization, enforcing sparsity

mat_sum = @(M) sum(M(:));

% parameters
alp   = params.alp;
bet   = params.bet;
lam   = params.lam;
sig   = params.noise.std;
A     = params.blur.op;
%matA  = params.blur.mat;
otfA  = params.blur.otf;
dim_e = params.dim_e;

% L: data fidelity
switch params.noise.type
    case {'none','Gaussian'}
        if strcmp(params.blur.type,'none')
            L.expr = @(u,z)     0.5*norms.L2(u-z).^2;
            L.prox = @(x,tau,z) prox.L2_dist_z(x,tau,z);
        else
            L.expr = @(u,z)     0.5*norms.L2(A(u)-z).^2;
            L.prox = @(x,tau,z) prox.L2_dist_z_resto(x,tau,z,otfA);
        end
    case 'Poisson'
        L.expr = @(u,z)     KLD(u,z,sig);
        L.prox = @(x,tau,z) prox.KLD(x,z,sig,tau);
end

% R: contour length penalization, enforcing sparsity
switch params.reg
    case 'l0'
        R.expr = @(e) norms.L0(e);
        R.prox = prox.L0;
    case 'l1'
        R.expr = @(e) norms.L1(e);
        R.prox = prox.L1;
    case 'l1q'
        R.expr = @(e) norms.quadL1(e);
        R.prox = prox.quadL1;
end

% S: smoothing term, preserving edges
S.expr = @(u,e) mat_sum((repmat(1-e,[1 1 2/dim_e(1) size(u,3)/dim_e(2)]).*D(u)).^2);

switch params.edges
    case 'similar',  S.D_sq = @(u) sum(D(u).^2,4);
    case 'distinct', S.D_sq = @(u) D(u).^2;
end
S.du = @(u,e)  2*Dadj(repmat((1-e).^2,[1 1 2/dim_e(1) size(u,3)/dim_e(2)]).*D(u)); 
S.de = @(u,e) -2*(1-e).*S.D_sq(u);

% energy
J = @(u,e,z) alp*L.expr(u,z) + bet*S.expr(u,e) + lam*R.expr(e);

end

function [u0,e0] = initialization(z,params)

fprintf(1,'Initializing u=%s and e=%s ...\n',params.init.u,params.init.edges);

switch params.init.u
    case 'rand',     u0 = randn(size(z));
    case 'rand_cst', u0 = randn*ones(size(g));
    case 'data',     u0 = z;
end

switch params.init.edges
    case '1',        e0 = ones([size(z,1),size(z,2),params.dim_e]);
    case '0',        e0 = zeros([size(z,1),size(z,2),params.dim_e]);
    case 'binomial', e0 = binornd(1,0.5,[size(z,1),size(z,2),params.dim_e]);
    case 'rand_cst'
        tmp = abs(0.25*randn);
        while tmp > 1., tmp = abs(0.25*randn); end
        e0 = tmp*ones([size(z,1),size(z,2),params.dim_e]);
end

end

function norms = reg_norms(eps)

mat_sum = @(M) sum(M(:));

% primal norms
norms.L0      = @(x) sum(double(x(:) ~= 0));

norms.L1      = @(x) sum(abs(x(:))); 
norms.quadL1  = @(x) mat_sum(max(abs(x(:)),x(:).*x(:)./eps));   

norms.L2      = @(x) sqrt(sum(x(:).^2)); 
norms.L2sq    = @(x) sum(x(:).^2); 


% dual norms
norms.L2_dual = @(y) sqrt(sum(y.^2,3));

end

function prox = proximity_operators(eps)

% primal usual prox
prox.L0     = @(x,tau) x.*double(abs(x)>sqrt(2*tau)); % [Attouch, Bolte, Svaiter 2013]

prox.L1     = @(x,tau) x - max(min(x,tau),-tau);
prox.quadL1 = @(x,tau) max(0,min(abs(x)-tau,max(eps,abs(x)./(tau/(eps/2)+1)))).*sign(x);

prox.L2sq   = @(x,tau) x./(1+2*tau);
                            
% Kullback-Leibler divergence 
prox.KLD    = @(x,g,sig,tau) (x - tau*sig + sqrt(abs(x-tau*sig).^2 + 4*tau*g))/2; 

% squared euclidean distance to z ||A.-g||_2^2
prox.L2_dist_z       = @(x,tau,z)      (x+tau.*z)./(1+tau);    
prox.L2_dist_z_resto = @(x,tau,z,otfA) ...
   real(ifft2((fft2(x)+tau*conj(otfA).*fft2(z))./(tau*conj(otfA).*otfA + 1)));
end

function y = D(x)
% Finite difference operator 
 
y = zeros(size(x,1),size(x,2),2,size(x,3));

% horizontal differences
y(:,:,1,:) = [x(:,2:end,:)-x(:,1:end-1,:),zeros(size(x,1),1,size(x,3))]./2.;
% vertical differences
y(:,:,2,:) = [x(2:end,:,:)-x(1:end-1,:,:);zeros(1,size(x,2),size(x,3))]./2.;  

y = squeeze(y); % remove singleton dimensions
  
end

function y = Dadj(x)
% Adjoint operator to the finite difference operator 
y1 = [x(:,1,1,:), x(:,2:end-1,1,:)-x(:,1:end-2,1,:), -x(:,end-1,1,:)]./2.;
y2 = [x(1,:,2,:); x(2:end-1,:,2,:)-x(1:end-2,:,2,:); -x(end-1,:,2,:)]./2.;

y = - y1 - y2;
y = squeeze(y); % remove singleton dimensions
end

function kld = KLD(u,z,sig)
% Computes the Kullback-Leibler divergence of a fonction u with source z

nn = 200.;

log2u = log2(u); 
log2u(isnan(log2u)) = 0.;
log2u(isinf(log2u)) = 0.;

kld = double(z>0  & u>0) .*(-z.*log2u + sig*u) ...
    + double(z>0  & u<=0).*nn ...
    + double(z==0 & u>=0).*(sig*u) ...
    + double(z==0 & u<0) .*nn;

kld = sum(kld(:));

end

function [z,noise_out] = addNoise(z0,noise)

if ~strcmp(noise.type,'none')
    fprintf(1,'Adding %s noise with parameter %g ...\n',...
                                      noise.type,noise.std);
end

switch noise.type
    case 'none',     z = z0;
    case 'Gaussian', z = z0 + (noise.std).*randn(size(z0));
    case 'Poisson',  z = random('Poisson',(noise.std)*z0);
end

noise_out = noise;

end

function [z,blur_out] = addBlur(z0,blur)

blur_out = blur;
[h,w,c] = size(z0);
    
switch blur.type
    case 'none'
        fs = 0.;
        F  = 1.;
    case 'Gaussian'
        fprintf(1,'Adding %s blur of size %ix%i with std %g ...\n',...
                                blur.type,blur.size,blur.size,blur.std);
        fs     = (blur.size-1)/2;
        [x,y]  = meshgrid(-fs:fs,-fs:fs);
        arg    = -(x.*x + y.*y)/(2*blur.std*blur.std);
        
        F = exp(arg);
        %F(F<eps*max(F(:))) = 0;
        F = F./sum(F(:));
    case 'Uniform'
        fprintf(1,'Adding %s blur of size %ix%i ...\n',...
                               blur.type,blur.size,blur.size);
        fs = (blur.size-1)/2;
        F  = ones(blur.size,blur.size,'double');
        F  = F/norm(F(:));
end

% blur matrix

sz = h*w;
D  = [];
for it = 1:blur.size
    C = repmat(F(:,it)',[sz,1]);
    if it ~= blur.size
        D = cat(2, D, C, zeros(sz,h-2*fs-1));
    else
        D = cat(2, D, C);
    end
        
end

% N  = (blur.size-1)/2;
% blur_out.mat = spdiags(D,-N*h-fs:N*h+fs,sz,sz);

blur_out.otf   = psf2otf(F,[h,w]);  blur_out.otf = repmat(blur_out.otf,[1 1 c]);
blur_out.op    = @(x) real(ifft2(blur_out.otf.*fft2(x)));  
blur_out.opAdj = @(x) real(ifft2(conj(blur_out.otf).*fft2(x)));  

z = blur_out.op(z0);

end

% function zn = normalize(z)
% 
% z  = z - min(z(:));
% zn = z./max(z(:));
% 
% end

function msg = display_progress(nb,nb_tot,format,varargin)
    
if nargin == 4, fprintf(repmat(sprintf('\b'), 1, length(varargin{1})));
else, fprintf('');
end
if ~strcmp(format,'ratio') && ~strcmp(format,'percentage'), format = 'ratio'; end

switch format
    case 'ratio'
        msg = sprintf('(%2i/%2i)', nb, nb_tot); % don't forget this semicolon
        if nb==0 || nb==1, fprintf(msg);
        else, fprintf([repmat(sprintf('\b'), 1, length(msg)), msg]);
        end

    case 'percentage'
        msg = sprintf('(%3i', ceil(nb/nb_tot*100)); % don't forget this semicolon
        if nb==0 || nb==1, fprintf(msg); fprintf('%%)');
        else, fprintf([repmat(sprintf('\b'), 1, length(msg)+2), msg]); 
              fprintf('%%)');
        end
        msg = [msg '%)'];
end

end

function [im,bet,lam,reg,init,blur,noise,edges,N,eps,alp] = parseInputs(varargin)

validImageTypes = {'uint8','uint16','int16','single','double'};

im = varargin{1};   % data
validateattributes(im,validImageTypes,{'nonsparse','real'},'parseInputs','IM',1);
if (ndims(im) > 3)
    error('IM cannot have more than 3 dimensions.');
end

bet = varargin{2};  % beta: smoothing parameter
validateattributes(bet,{'numeric'},{'nonnan','nonnegative'},'parseInputs','BETA',2);

lam = varargin{3};  % lambda: contours length penalization
validateattributes(lam,{'numeric'},{'nonnan','nonnegative'},'parseInputs','LAMBDA',3);


% Default values for parameters
reg = 'l1q';

init.u     = 'data';
init.edges = '1';

blur.type = 'none'; 
blur.size = 0.;

blur.std  = 0.;

noise.type = 'none';
noise.std  = 0.;

edges = 'similar';

N   = 10000;
eps = .5;
alp = 1.;

% check inputs validity
args_names = {'AddBlur','AddNoise','InitU','InitEdges','Regularization',...
              'Edges','MaxNumberIteration','QuadL1Eps','Alpha'};
           
for i = 4:2:nargin
    arg = varargin{i};
    if ischar(arg)        
        idx = find(strncmpi(arg, args_names, numel(arg)));
        if isempty(idx)
            error(['Unknown input string: ',arg]); 
            
        elseif numel(idx) > 1
            error(['Ambiguous input string: ',arg]);
            
        elseif numel(idx) == 1
            if (i+1 > nargin) 
                error('One or more parameter name or value is missing.');             
            end
            if idx == 1
                addblur = varargin{i+1};
                validateattributes(addblur,{'numeric'},...
                    {'vector','nonempty','>=',0},'parseInputs',varargin{i},i);
                if length(addblur) < 2 || length(addblur) > 3 
                    error(['Expected input number' i ', ' varargin{i} ...
                        ', to be an array with number of elements equal to 2 or 3.']);  
                elseif addblur(1)==1 && length(addblur) < 3
                    error(['Expected input number' i ', ' varargin{i} ...
                        ', to be an array with number of elements equal to 3.']);  
                elseif addblur(1)==2 && length(addblur)==3
                    error(['Expected input number' i ', ' varargin{i} ...
                        ', to be an array with number of elements equal to 2.']);  
                else
                    switch addblur(1)
                        case 0, blur.type = 'none';
                        case 1, blur.type = 'Gaussian';
                        case 2, blur.type = 'Uniform';
                        otherwise, error(['Unknown blur identifier, expected '...
                                        'values are: 1 (Gaussian), 2 (Uniform).']); 
                    end
                    
                    validateattributes(addblur(2),{'numeric'},...
                        {'scalar','odd','positive'},'parseInputs','blur.size');
                    blur.size = addblur(2);
                    
                    if strcmp(blur.type,'Gaussian')
                        blur.std = addblur(3);
                    end
                end
                
            elseif idx == 2
                addnoise = varargin{i+1};
                validateattributes(addnoise,{'numeric'},{'vector','nonempty',...
                    '>=',0,'numel',2},'parseInputs',varargin{i},i);
                
                switch addnoise(1)
                    case 0, noise.type = 'none';
                    case 1, noise.type = 'Gaussian';
                    case 2, noise.type = 'Poisson';
                    otherwise, error(['Unknown noise identifier, expected ' ...
                                    'values are: 1 (Gaussian), 2 (Poisson).']);
                end
                
                noise.std = addnoise(2);
                              
            elseif idx == 3
                init.u = varargin{i+1};
                validatestring(init.u,{'data','rand','rand_cst'},'parseInputs',...
                    varargin{i},i);          
                
            elseif idx == 4
                init.edges = varargin{i+1};
                validatestring(init.edges,{'0','1','binomial','rand_cst'},...
                    'parseInputs',varargin{i},i);               
                
            elseif idx == 5
                reg = varargin{i+1};
                validatestring(reg,{'l0','l1','l1q'},'parseInputs',varargin{i},i);         
                
            elseif idx == 6
                edges = varargin{i+1};
                validatestring(edges,{'similar','distinct'},'parseInputs',...
                    varargin{i},i);        
                
            elseif idx == 7
                N = round(varargin{i+1});
                validateattributes(N,{'numeric'},{'nonnan','positive'},...
                    'parseInputs',varargin{i},i);        
                
            elseif idx == 8
                eps = varargin{i+1};
                validateattributes(eps,{'numeric'},{'nonnan','positive'},...
                    'parseInputs',varargin{i},i); 
                
            elseif idx == 9
                alp = varargin{i+1};
                validateattributes(alp,{'numeric'},{'nonnan','positive'},...
                    'parseInputs',varargin{i},i);         
            end
        end    
    else
        error('Parameter name must be specified.'); 
    end
end

if strcmp(noise.type,'Poisson') && blur.size > 1
    error('Deblurring with Poisson noise is not implemented yet.');
end

end