%% Visualization of RGB components of an Hyperspectral image
%
% Author:
% Daniele Picone
%
% Description:
% This script returns the visualization of the RGB components taken from a
% given hyperspectral acquisition. Specifically, this script introduces a
% given illuminant, clipping and gamma correction
% (default: sRGB specifications)
% 
% Usage:
% I_out=show_RGB(I_in,'Field',Value)
% Eg:
% I_out=show_RGB(I_in,'illum',[0.01,0.99],'clip',2.2,'illum','D65','wavelength',420:10:720);
%
% Input:
% I_in: A struct with the input, whose fields are
%   - data: The data containing the image cube (3D matrix whose dimensions
%           are [Spatial Pixel Vertical, Spatial Pixel Horizontal, Bands])
%   - wavelength: The central wavelengths associated to the bands
%
% Optional Fields
% 'illum': type of illumination standard to be used


function I_out=show_RGB(I_in,varargin)

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'Color_space'));

wl=[];
label_illum='D65';
clip_level=[0.01,0.99];
gamma=2.2;

if isfield(I_in,'wavelength'), wl=I_in.wavelength; end
if isfield(I_in,'data'), I_load=I_in.data; end
if isnumeric(I_in), I_load=I_in; end
clear I_in;

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if any(strcmpi(pname,{'illum','illuminant','D'}))        % Illuminant
        label_illum=pval;
    elseif any(strcmpi(pname,{'gamma','gamma_correction'}))  % (Inverse) gamma correction level
        gamma=pval;
    elseif any(strcmpi(pname,{'wavelength','wl'}))           % Central wavelength of the input
        wl=pval;
    elseif any(strcmpi(pname,{'clip','clipping','clip_level'})) % Percentage pectral coverage of MS over PAN
        clip_level=pval;
    end
end

if isempty(wl), error('Information on central wavelengths is required'); end

if nargin<=0, wl=400:10:720; end
if nargin<=1, label_illum='D65'; end

if size(wl,1)==1, wl=wl.'; end

% load('ref4_scene4.mat','reflectances');
% size(reflectances);
% It should have size 255 x 335 x 33. To inspect reflectances, try taking a slice at a middle wavelength, say 560 nm (17th in the sequence 400, 410, , 720 nm), and displaying it. The image should look like Fig. 4 (except for the added red square, left). 

% slice = reflectances(:,:,17);
% figure; imagesc(slice); colormap('gray'); brighten(0.5);

% z = max(slice(100, 39));
% slice_clip = min(slice, z)/z;
% figure; imagesc(slice_clip.^0.4); colormap('gray');

% reflectance = squeeze(reflectances(141, 75,:));
% figure; plot(wl, reflectance);
% xlabel('wavelength, nm');
% ylabel('unnormalized reflectance');


I_load=I_load./max(I_load,[],1:3);

illum=illuminant_D(wl,label_illum);

radiances=I_load.*shiftdim(illum,-2);

% radiances_pixel = squeeze(radiances(141, 75, :));
% figure; plot(400:10:720, radiances_pixel,'b'); % blue curve
% xlabel('wavelength, nm');
% ylabel('radiance, arbitrary units');
% hold on;

% CIE color matching function
% Reference: http://www.cvrl.org/ (Section: Older CIE Standards)
load ('colormatch_cie31.mat','Wavelength','xyzbar'); % CIE 1931
%load ('colormatch_cie64.mat','Wavelength','xyzbar'); % CIE 1965
%load ('colormatch_xyz2e.mat','Wavelength','xyzbar'); % CIE 2006 (2 degrees aperture)
%load ('colormatch_xyz10.mat','Wavelength','xyzbar'); % CIE 2006 (10 degrees aperture)
 
radiances = permute(radiances, [3,1,2,4:ndims(radiances)]);
size_radiances=size(radiances);
radiances = reshape(radiances,size(radiances,1),[]);
xyzbar_new=zeros(length(wl),3);
for ii =1:3
    xyzbar_new(:,ii)=interp1(Wavelength,xyzbar(:,ii),wl);
end
xyzbar_new(wl<Wavelength(1)|wl>Wavelength(end),:)=0;

XYZ = xyzbar_new.'*radiances;
XYZ = reshape(XYZ,[3,size_radiances(2:end)]);
XYZ = ipermute(XYZ,[3,1,2,4:ndims(radiances)]);
XYZ = XYZ./max(XYZ,[],1:3);


RGB = XYZ2sRGB_exgamma(XYZ);
% RGB = min(max(RGB, 0),1); %thresholding
% z = max(RGB(244,17,:));

size_RGB = size(RGB);
RGB = reshape(RGB,size(RGB,1)*size(RGB,2)*size(RGB,3),[]);
RGB_sort= sort(RGB,1);
z=RGB_sort( ceil (clip_level(2)*size(RGB_sort,1)), :);
y=RGB_sort( floor(clip_level(1)*size(RGB_sort,1)), :);
I_out = (max(min(RGB, z),y)-y)./(z-y);
I_out=reshape(I_out,size_RGB);
I_out=I_out.^(1/gamma);

% RGB_clip=viewimage_outputonly(max(min(RGB,1),0));

figure; imshow(I_out);

end