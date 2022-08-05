%% CLASSIC DEMOSAIC
%
% Description:
% This function performs classical demosacking for images obtained by
% mosaicking an image via a Bayer mask
%
% Usage:
% I_out=demosaic_classic(I_in,mask,method,DynamicRange,central)
%
% Input:
% I_in: Mosaicked image (size: N1 x N2)
% mask: Binary indicator of the positions of the original samples
%       in the full image (size: N1 x N2 x Nb).
%       Note: if empty, it returns the input as is
% method: Demosaicking method, to choose amongst (cell of multiple methods is allowed):
%    - For Bayer masks:
%       'RI':    Residual Interpolation [1]
%       'MLRI':  Minimized-Laplacian Residual Interpolation [2]
%       'MLRI2': Variant of the above
%       'ARI':   Adaptative Residual Interpolation [3]
%       'ARI2':  Variant of the above
%       'AP':    Alternating Projections [6]
%       'MSG':   Multiscale Gradients-Based CFA Interpolation [7]
%    - For periodic masks:
%       'SD':    Spectral Difference [4]
%       'ISD':   Iterative Spectral Difference [4]
%       'ID':    Intensity Difference [4]
%       'IID':   Iterative Intensity Difference [4]
%    - Scattered Interpolation (for non overlapping mask):
%       'WB':    Deluney triangulation bilinear interpolation [5]
%       'NN':    Nearest Neighbour Interpolation
%       'RBF':   Radial Basis Filter Interpolation [5]
%                ('RBF_gaussian','RBF_iq','RBF_spline','RBF_mq','RBF_imq')
%       'IDW':   Inverse Distance Weighting (Shepard) Interpolation [5]
% DynamicRange: Maximum input value (default: 2^nextpow2(max(I_in(:))+1)-1)
% central: Central Wavelengths vector in nm (default: Visible range)
%
% References:
% [1] Kiku D., Monno Y., Tanaka M. and Okutomi M., "Beyond Color
%     Difference: Residual interpolation for color image demosaicking",
%     IEEE TIP, vol.25 no. 3, pp. 1288-1300, Mar. 2016
% [2] Kiku D., Monno Y., Tanaka M. and Okutomi M., "Minimized-Laplacian 
%     Residual Interpolation for Color Image Demosaicking" in In Digital 
%     Photography X (Vol. 9023, p. 90230L). ISOP.
% [3] Monno Y., Kiku D., Tanaka M. and Okutomi M., "Adaptive Residual 
%     Interpolation for Color and Multispectral Image Demosaicking", 
%     Sensors, Vol.17, No.12, pp.2787-1-21, December, 2017.
% [4] Mihoubi S., Losson O., Mathon B. and Macaire L.,  "Multispectral 
%     demosaicing using pseudo-panchromatic image", IEEE TCI, IEEE, 2017, 
%     3 (4), pp.982-995.
% [5] Amidror I. , “Scattered data interpolation methods for electronic 
%     imaging systems: a survey.” Journal of Electronic Imaging. Vol. 11, 
%     No. 2, April 2002, pp. 157–176.
% [6] Lu Y. M., Karzand M., and Vetterli M., "Demosaicking by Alternating 
%     Projections: Theory and Fast One-Step Implementation," IEEE TIP, vol.
%     19, no. 8, August 2010.
% [7] Pekkucuksen I. and Altunbasak Y., "Multiscale Gradients-Based Color
%     Filter Array Interpolation", IEEE TIP, vol.22, no.1, pp.157-165, Jan.
%     2013.

function [Iout,time]=demosaic_classic(Iin,mask,method,DynamicRange,central)

current_folder=pwd;
algorithm_folder=fileparts(mfilename('fullpath'));
addpath(fullfile(algorithm_folder,'check_functions'));

%% Parse inputs
if isfield(Iin,'data'), I_in=Iin.data; else, I_in=Iin; end

if nargin<=2 || isempty(method), method='WB'; end
if ~iscell(method), method={method}; end

if nargin<=1 || isempty(mask)
    for ii=1:min(1,numel(method)), method{ii}='none'; end
    mask=ones(size(I_in));
end
if isfield(mask,'data'), mask=mask.data; end

if nargin<=3 || isempty(DynamicRange)
    if isfield(Iin,'DynamicRange'), DynamicRange=Iin.DynamicRange;
    else, DynamicRange=2.^nextpow2(max(I_in(:))+1)-1; 
    end
end
if nargin<=4 || isempty(central)
    if isfield(Iin,'wavelength'), central=Iin.wavelength;
    else, step=(700-400)/size(mask,3); central=(400+step/2:step:700).';
    end
end

%% Demosaic method redirecting

I_out=zeros([size(mask,1),size(mask,2),size(mask,3),numel(method)]);
time=zeros(numel(method),1);
fprintf('Current Demosaic Algorithm:  0/%2d',numel(method));
for jj=1:numel(method)
    if any(strcmpi(method{jj},{'none'}))
        I_out(:,:,:,jj)=I_in;
    elseif any(strcmpi(method{jj},{'RI','MLRI','MLRI2','ARI','ARI2'}))
        [mask_temp,I_temp,options]=Bayer_check(mask,I_in);
        addpath(fullfile(algorithm_folder,'RI'));
        [I_out_temp,time(jj)] = demosaic_algorithm_RI(method{jj},I_temp,mask_temp,'grbg',DynamicRange);
        I_out(:,:,:,jj)=Demosaic_reshaper(mask_temp,I_out_temp,options);
    elseif strncmpi(method{jj},'RBF',3) || strncmpi(method{jj},'IDW',3) || any(strcmpi(method{jj},{'WB','NN'}))
        addpath(fullfile(algorithm_folder,'..','Interpolation'));
        I_out(:,:,:,jj)=interp2D_mosaic(I_in,'mask',mask,'method',method{jj},'nn',100);
    elseif any(strcmpi(method{jj},{'SD','ISD','ID','IID'}))
        [mask_temp,I_temp,period]=Period_check(mask,I_in);
        cd(fullfile(algorithm_folder,'ID'));
        I_out(:,:,:,jj)=demosaic_algorithm_ID(I_temp,mask_temp,method{jj},period,central);
        cd(current_folder);
    elseif strcmpi(method{jj},'MSG')
        [mask_temp,I_temp,options]=Bayer_check(mask,I_in);
        I_temp=I_temp/DynamicRange*255;
        cd(fullfile(algorithm_folder,'MSG'));
        I_out_temp=MSG_Bayer(I_temp);
        cd(current_folder);
        I_out_temp=I_out_temp/255*DynamicRange;
        I_out(:,:,:,jj)=Demosaic_reshaper(mask_temp,I_out_temp,options);
    elseif strcmpi(method{jj},'AP')
        [mask_temp,I_temp,options]=Bayer_check(mask,I_in);
        cd(fullfile(algorithm_folder,'AP'));
        I_temp=I_temp/DynamicRange*255;
        I_out_temp=demosaic_algorithm_AP(I_temp);
        I_out_temp=I_out_temp/255*DynamicRange;
        cd(current_folder);
        I_out(:,:,:,jj)=Demosaic_reshaper(mask_temp,I_out_temp,options);
    end
    fprintf('\b\b\b\b\b%2d/%2d',jj,numel(method));
end
fprintf(' Done!\n');


%% Create struct for output data
if isstruct(Iin), Iout=Iin; end
Iout.data=I_out;
Iout.label='Demosaicked';
Iout.methods_demosaic=method;
  