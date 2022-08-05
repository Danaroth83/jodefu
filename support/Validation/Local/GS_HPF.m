function [PanSharpenedImage, Coeff] = GS_HPF(I_MS,I_PAN,I_LP,S,equaliz,interp)
%
%
%
% Possible ways to handle the calculation on the blocks:
%  - if interp is not given or is empty, the coefficients are computed on
%  blocks that overlaps (the coeff computed on a block will be stored in
%  only its central pixel - i.e. each pixel might have different coeff values).
%  - if interp is given (can be 0: no interp, 1: bilinear, 2: bicubic) the
%  coefficients will be computed on each block and stored on all the pixels
%  belonging to each block. There is no overlap on the blocks. The
%  coefficients can be interpolated according to what specified by interp.

% keyboard


addpath('ps_util')
overlap = false;
if nargin < 6
    interp=0;
    overlap=true;
    fprintf('Overlapping ');
    if S==1,
        fprintf('HPM method\n');
    end
end

if nargin<5
    equaliz=0;
end

if nargin<4
    S = ones(size(I_PAN));
end

I_MS = double(I_MS);
[Height,Width,Bands] = size(I_MS);
ratio=4;
%% Generate I_LP
% switch GS_mode
%     case 'GS1'
%         fprintf('GS1 [Laben00]\n');
%     case 'GS2'
%         fprintf('GS2 [Laben00, Aiazzi07]\n');
%     case 'GSA'
%         fprintf('GSA [Aiazzi07]');
% end

%% generate the details
I_PAN = repmat(double(I_PAN), [1, 1, Bands]);
I_LP = double(I_LP);
if size(I_LP, 3) == 1
    I_LP = repmat(I_LP, [1, 1, Bands]);
end
if size(I_LP, 3) ~= size(I_PAN, 3)
    error('I_LP should have the same number of bands as PAN');
end

% Equalization PAN with I_LP
% keyboard
switch equaliz
    case 1
        for ii = 1:Bands,
            imageHR(:,:,ii) = I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii))+ mean2(I_LP(:,:,ii));
        end
    case 2
        for ii = 1:Bands,
            imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                std2(I_LP(:,:,ii))/std2(I_PAN(:,:,ii)) + mean2(I_LP(:,:,ii));
        end
    case 3
        for ii = 1:Bands,
            imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                std2(I_LP(:,:,ii))/std2(imresize(I_PAN(:,:,ii),1/ratio)) +...
                mean2(I_LP(:,:,ii));
        end
    case 4
        for ii = 1:Bands,
            imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                std2(I_LP(:,:,ii))/std2(imresize(imresize(I_PAN(:,:,ii),1/ratio),ratio)) +...
                mean2(I_LP(:,:,ii));
        end
    otherwise
        imageHR = I_PAN;
end
% keyboard
DetailsHRPan = imageHR - I_LP;

%% Coefficient computation

Coeff = zeros(size(I_MS));

for ii = 1:Bands,
    MS_Band = squeeze(I_MS(:,:,ii));
    I_LP_Band = squeeze(I_LP(:,:,ii));
    Coeff_Band=ones(size(I_LP_Band));
    
    %%%%%%%%%%Fusione
    % interpolation of the coefficients
    if ~overlap && interp ~= 0  %Anche con segmentation
        if isempty(BlockSize)
            BlockSize=4;
        end
        Coeff1 = imresize(Coeff_Band, 1/BlockSize, 'nearest');
        if interp==1
            Coeff_Band = imresize(Coeff1, [size(Coeff_Band,1),size(Coeff_Band,2)], 'bilinear');
        else
            Coeff_Band = imresize(Coeff1, [size(Coeff_Band,1),size(Coeff_Band,2)], 'bicubic');
        end
        clear Coeff1
    end
    
    Coeff(:,:,ii)=Coeff_Band;
    
end
% interpolation of the coefficients
% keyboard
% if (~isscalar(S) || (~isempty(BlockSize) && ~overlap)) && interp ~= 0  %Forse anche con segmentation
%     if isempty(BlockSize)
%         BlockSize=4;
%     end
%     naddrows = (ceil(size(Coeff,1)/BlockSize)+1) * BlockSize - size(Coeff,1);
%     naddcols = (ceil(size(Coeff,2)/BlockSize)+1) * BlockSize - size(Coeff,2);
%     Coeff1 =zeros([size(Coeff,1)+naddrows, size(Coeff,2)+naddcols, Bands]);
%     Coeff1(1:size(Coeff,1),1:size(Coeff,2),:) = Coeff;
%     Coeff1(size(Coeff,1)+1:end,:,:) = repmat(Coeff1(size(Coeff,1),:,:), [naddrows, 1, 1]);
%     Coeff1(:,size(Coeff,2)+1:end,:) = repmat(Coeff1(:,size(Coeff,2),:), [1, naddcols, 1]);
%     Coeff1(size(Coeff,1)+1:end,size(Coeff,2)+1:end,:) = repmat(Coeff1(size(Coeff,1),size(Coeff,2),:), [naddrows, naddcols, 1]);
%     Coeff1 = imresize(Coeff1, 1/BlockSize, 'nearest');
%     %     Coeff2 = imresize(Coeff1, BlockSize, 'nearest');% for check
%
%     if interp==1
%         Coeff = imresize(Coeff1, BlockSize, 'bilinear');
%     else
%         Coeff = imresize(Coeff1, BlockSize, 'bicubic');
%     end
%     clear Coeff1
%     Coeff = Coeff(naddrows:naddrows+Height-1,naddrows:naddcols+Width-1,:);
% end
%% Fusion
% keyboard
PanSharpenedImage = Coeff .* DetailsHRPan + I_MS;

%
% Final Mean Equalization
for ii = 1 : size(I_MS,3)
    h = PanSharpenedImage(:,:,ii);
    PanSharpenedImage(:,:,ii) = h - mean2(h) + mean2(I_MS(:,:,ii));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  End  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

