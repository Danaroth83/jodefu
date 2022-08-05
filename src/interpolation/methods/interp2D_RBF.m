% RADIAL BASIS FUNCTION 2D SCATTERED INTERPOLATION
%
% Description:
% This script works as frontend for the toolbox of the Radial Basis
% Function (RBF) [1], applied for 2D interpolation, for images in a grid
% with missing samples, marked by the zeros in the mask. The script also
% supports multidimensional images, but the interpolation is just performed
% over each 
%
% Usage:
% fo=interp2D_RBF(xi,yi,fi,xo,yo,sigma,type);
% 
% Input:
% xi,yi: Column vectors for the input sample coordinates
% fi: Column vector for the input values
% xo,yo: Column vectors for the output coordinates
% Nsamples: Number of samples to be considered in the RBF surroundings
%           (default: amount of non zero samples)
% sigma: Parameter of the RBF; it's the reciprocal of the variance for a
%        Gaussian RBF (default: 1; can be inputted as vector of sizes 1 x Ns)
% type: The shape of the RBF (default: 'Gaussian')
%    -'Gaussian':              RBF(r)=exp(-(sigma*r)^2)
%    -'InverseQuadratic':      RBF(r)=1/(1+(sigma*r)^2)
%    -'MultiQuadratic':        RBF(r)=sqrt(1+(sigma*r)^2)
%    -'InverseMultiQuadratic': RBF(r)=1/sqrt(1+(sigma*r)^2)
%    -'ThinPlate':             RBF(r)=(sigma*r)^2*ln(sigma*r)
%
% Output:
% I_out: Reconstructed image (If the image is monocromatic it has 
%        sizes: Nx x Ny x Ns, otherwise Nx x Ny x Nb x Ns)
%
% Reference:
% [1] Sarra S., "The Matlab Radial Basis Function Toolbox". Journal of 
% Open Research Software, 5: 8,

function fo=interp2D_RBF(xi,yi,fi,xo,yo,Nsamples,sv,type)

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'RBF'));

if nargin<=5 || isempty(Nsamples), Nsamples=prod(size(I_in,1),size(I_in,2)); end
if nargin<=6 || isempty(sv), sv=1; end
if nargin<=7 || isempty(type), type='Gaussian'; end

%% Loading RBF

mu = 1.5e-13;

if any(strcmpi(type,{'Gaussian','gax'}))
    phi = gax();
    safe = false; 
elseif any(strcmpi(type,{'InverseQuadratic','iqx','iq'}))
    phi = iqx();
    safe = false;
elseif any(strcmpi(type,{'InverseMultiQuadratic','imqx','imq'}))
    phi = imqx();
    safe = false;
elseif any(strcmpi(type,{'MultiQuadratic','mqx','mq'}))
    phi = mqx();
    safe = true;
elseif any(strcmpi(type,{'ThinPlate','Spline','tpx','tp'}))
    phi = tpx();
    safe = true; % Matrix is not semidefinite positive
else
    error('Unknown RBF')
end

Nsamples=min(Nsamples,length(xi));
fo=zeros(length(xo),length(sv));
fprintf('Interpolation sample:     0/%5d',length(xo))
for jj=1:length(xo)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d/%5d',jj,length(xo));
    % rf_ss=phi.distanceMatrix2d(xi,yi,xo(ii),yo(ii));
    rf_ss=(sqrt((xi-xo(jj)).^2+(yi-yo(jj)).^2)).';
    [rf_ss,idx_ss]=mink(rf_ss,Nsamples,2);

    ri_ss = phi.distanceMatrix2d(xi(idx_ss),yi(idx_ss));
    fi_ss=fi(idx_ss);
    for ii=1:length(sv)
        s = sv(ii); % fprintf('s=%.1f\n',s);      
        B = phi.rbf(ri_ss,s);
        a = phi.solve(B,fi_ss,mu,safe);
        H = phi.rbf(rf_ss,s);
        fo(jj,ii) = H*a;
    end
end
fprintf('. Done!\n');
