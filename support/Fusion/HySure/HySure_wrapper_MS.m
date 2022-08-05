function Zimhat = HySure_wrapper_MS( HS, MS, downsamp_factor, overlap,center_points,im_tag )
%HySure_wrapper This is a wrapper function that calls HySure [1], to be
% used for the comparisons in [2]. The original source code, together with
% other datasets, can be obtained from https://github.com/alfaiate/HySure
% See the file README for more information.
% 
% Zimhat = HySure_wrapper( HS, PAN, downsamp_factor, overlap )
% 
% Input: 
% HS: HS image, 
% PAN: PAN image, 
% downsamp_factor: downsampling factor,
% overlap: array vector with the spectral coverage of both sensors, 
%     i.e., which HS bands are covered by the PAN band - optional
% 
% Output: 
% Zimhat: estimated image with high spatial and spectral resolution
%
% 
%   [1] M. Simoes, J. Bioucas-Dias, L. Almeida, and J. Chanussot, 
%        “A convex formulation for hyperspectral image superresolution via 
%        subspace-based regularization,” IEEE Trans. Geosci. Remote Sens.,
%        to be publised.
% 
%   [2] Laetitia Loncan, Luis B. Almeida, Jose M. Bioucas-Dias, Xavier Briottet, 
%        Jocelyn Chanussot, Nicolas Dobigeon, Sophie Fabre, Wenzhi Liao, 
%        Giorgio A. Licciardi, Miguel Simoes, Jean-Yves Tourneret, 
%        Miguel A. Veganzones, Gemine Vivone, Qi Wei and Naoto Yokoya, 
%        "Introducing hyperspectral pansharpening," Geoscience and Remote Sensing
%        Magazine, 2015.

% % % % % % % % % % % % % 
% 
% Version: 1
% 
% % % % % % % % % % % % % 
% 
% Copyright (C) 2015 Miguel Simoes
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 

if nargin<=5, im_tag='Paris'; end
% % % % % % % % % % % % % 
% 
% This function follows two steps. 
% 
% I. It starts by estimating the spectral and spatial response of the sensors.
% The regularization parameters can be adjusted here:
lambda_R = 1e1;
lambda_B = 1e1;
% For the denoising with SVD, we need to specify the number of bands we
% want to keep
p = min(10,numel(cell2mat(overlap))); 
% Corresponds to variable L_s in [1]; number of endmembers in VCA /
% number of non-truncated singular vectors
%
% 
basis_type = 'VCA';
if strncmpi(im_tag,'Moffett',8)
    lambda_phi = 1e-3;
else
    lambda_phi = 5e-4;
end
lambda_m = 1;
% 
% ms_bands = 1; % Only panchromatic (HySure works for multispectral images 
%              as well, i.e., with more than one band)

% Normalize
%% Bioucas-Dias toolbox version
% max_HS = max(max(max(HS)));
% Yhim = HS./max_HS;
% Ymim = MS./max_HS;

[nl,nc,Nb_MS]=size(MS);
Ymim=zeros(nl,nc,Nb_MS);
for ii=1:Nb_MS
    MS_temp=MS(:,:,ii);
    xmax=quantile(MS_temp(:),0.999);
    Ymim(:,:,ii)=MS_temp/xmax;
end

[L1,L2,Nb]=size(HS);
Yhim=zeros(L1,L2,Nb);
max_HS=zeros(1,Nb);
for ii=1:Nb
    HS_temp=HS(:,:,ii);
    max_HS(ii)=quantile(reshape(HS_temp,[L1*L2,1]),0.999);
    Yhim(:,:,ii)=HS_temp/max_HS(ii);
end

% % Blur kernel
% middlel = round((nl+1)/2);
% middlec = round((nc+1)/2);
% % Blur matrix
% B = zeros(nl,nc);
% % Starck-Murtagh filter
% B(middlel-2:middlel+2, middlec-2:middlec+2) = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1];
% % Circularly center B
% B = ifftshift(B);
% % Normalize
% B = B/sum(sum(B));
% % Fourier transform of the filters
% FB = fft2(B);
% 
% % % % % % % % % % % % % % 
% % Simulate the HS data
% % Spatial degradation (blur)
% Z = im2mat(Zim);
% Yh = ConvC(Z, FB, nl);
% Yhim_up = mat2im(Yh, nl);
% % Add noise
% sigmah = sqrt(sum(Yhim_up(:).^2)/(10^(SNRh/10))/numel(Yhim_up));
% Yhim_up = Yhim_up + sigmah*randn(size(Yhim_up));
% % % % % % % % % % % % Downsampling (version with reduced size)
% % Downsampling
% Yhim = downsamp_HS(Yhim_up, downsamp_factor, 1);


% Define the spectral coverage of both sensors, i.e., which HS band
% corresponds to each MS band
% e.g.: MS band 1 - HS bands 1,2,3,4,5
%       MS band 2 - HS bands 6,7,8,9
%       MS band 3 - ...
% Now imagine that there are some bands that are very noisy. We remove them
% from the hyperspectral data cube and build a vector with the number of 
% the bands there were not removed (we call this vector 'non_del_bands').
% For example, if we removed bands 3 and 4, non_del_bands = [1,2,5,6,...]
% We now define a cellarray, called 'intersection',  with length(ms_bands) 
% cells. Each cell corresponds to a multispectral band and will have a 
% vector with the number of the hyperspectral bands that are covered by it.
% Since we removed some bands as well, we need to keep track of the bands
% that are contiguous in the data cube but are not in the sensor.
% We call this other cellarray 'contiguous'. If there are no
% removed bands, it can be set to be the same as 'intersection'.
% intersection = cell(1,length(ms_bands));
% if nargin == 3
%     intersection{1} = 1:size(HS, 3);
% elseif nargin == 4
%     intersection{1} = overlap; 
% else
%     disp('Please check the usage of HySure_wrapper');
% end

contiguous = overlap;
intersection= overlap;
% Blur's support: [hsize_h hsize_w]
hsize_h = 9;
hsize_w = 9;    
shift = round(mean(center_points))-1; % 'phase' parameter in MATLAB's 'upsample' function
blur_center = 0; % to center the blur kernel according to the simluated data
[V, R_est, B_est] = sen_resp_est(Yhim, Ymim, downsamp_factor, intersection, contiguous, p, lambda_R, lambda_B, hsize_h, hsize_w, shift, blur_center);

% II. The data fusion algorithm is then called using the estimated responses
% and the observed data. 

Zimhat = data_fusion(Yhim, Ymim, downsamp_factor, R_est, B_est, p, basis_type, lambda_phi, lambda_m);

% Denoise the data again with V - optional
Zhat = im2mat(Zimhat);
Zhat_denoised = (V*V')*Zhat;
% In image form
Zimhat = mat2im(Zhat_denoised, size(MS, 1));

% deNormalize
for ii=1:Nb
    Zimhat(:,:,ii) = Zimhat(:,:,ii).*max_HS(ii);
end

end
