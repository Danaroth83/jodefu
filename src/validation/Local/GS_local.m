function [PanSharpenedImage, Coeff] = GS_local(I_MS,I_PAN,I_LP,S,equaliz,interp)
%
%
% S can be:
%  - a scalar, in this case it is the block size and the conventional
%       algorithm will be taken (use the size of image to have 'global')
%  - a 1-band image, this is a segmentation map defining the regions in the
%       image and it will be considered for all the bands of MS
%  - a N-band image, this is a set of segmentation maps one associated to
%       each MS band
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
if isempty(interp)
    interp=0;
    overlap=true;
    fprintf('Overlapping ');
    if S==1,
        fprintf('HPM method\n');
    elseif S==0,
        fprintf('HPF method\n');
    end
end

% if nargin<5
%     equaliz=0;
% end
%
% if nargin<4
%     S = ones(size(I_PAN));
% end

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
        %        for ii = 1:Bands,
        %             imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
        %                 std2(I_MS(:,:,ii))/std2(I_PAN(:,:,ii)) + mean2(I_MS(:,:,ii));
        %         end
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

% Parse S
BlockSize = [];
if isscalar(S)
    BlockSize = S;
    fprintf('Block size %d, interp %d\n', S, interp);
elseif (size(S,1) == size(I_PAN,1)) && (size(S,2) == size(I_PAN,2))
    if size(S,3) == 1
        S = repmat(S, [1, 1, Bands]);
    elseif size(S,3) ~= Bands
        error('Number of bands in S is incorrect.');
    end
    fprintf('Segmentation map\n');
else
    error('S is not correctly shaped.');
end
%     ms_orig = imresize(I_MS,1/ratio);
%     ms_LP_d = MTF(ms_orig,sensor,[],ratio);

% keyboard
for ii = 1:Bands,
    MS_Band = squeeze(I_MS(:,:,ii));
    I_LP_Band = squeeze(I_LP(:,:,ii));
    Coeff_Band=zeros(size(I_LP_Band));
    if isscalar(S)   %%Blocks
        if overlap
            if S==1
                %%%% HPM method
                Coeff_Band=MS_Band./I_LP_Band;
            elseif S==0
                %%%% HPF method
                Coeff_Band=ones(size(MS_Band));
            else
                for y=1 : Height
                    for x=1 : Width
                        startx = max(x - floor(BlockSize/2), 1);
                        starty = max(y - floor(BlockSize/2), 1);
                        endy = min(y + floor(BlockSize/2), size(I_MS,1));
                        endx = min(x + floor(BlockSize/2), size(I_MS,2));
                        BlockMS = MS_Band(starty:endy,startx:endx);
                        BlockPan = I_LP_Band(starty:endy,startx:endx);
                        %   keyboard
                        c = cov(BlockPan(:),BlockMS(:));
                        Coeff_Band(y,x) = c(1,2)/var(BlockPan(:));
                        %%%% Correlation
                        %  Coeff_Band(y,x) = (BlockPan(:)'*BlockMS(:))/(BlockPan(:)'*BlockPan(:));
                    end
                end
            end
        else
            
            for y=1 : BlockSize : Height
                for x=1 : BlockSize : Width
                    %%More general, but it requires double value to have the whole area
                    %                     startx = max(x - floor(BlockSize/2), 1);
                    %                     starty = max(y - floor(BlockSize/2), 1);
                    %                     endy = min(y + floor(BlockSize/2), size(MS,1));
                    %                     endx = min(x + floor(BlockSize/2), size(MS,2));
                    %% Good for integer ratio between size of image and blocks
                    startx = max(x, 1);
                    starty = max(y, 1);
                    endy = min(y + BlockSize, size(I_MS,1));
                    endx = min(x + BlockSize, size(I_MS,2));
                    
                    BlockMS = MS_Band(starty:endy,startx:endx);
                    BlockPan = I_LP_Band(starty:endy,startx:endx);
                    
                    c = cov(BlockPan(:),BlockMS(:));
                    Coeff_Band(starty:endy,startx:endx) = c(1,2)/var(BlockPan(:));
                    %%%% Correlation
                    %                     Coeff_Band(starty:endy,startx:endx) = (BlockPan(:)'*BlockMS(:))/(BlockPan(:)'*BlockPan(:));
                    % keyboard
                end
            end
            % keyboard
            %             for y=1 : BlockSize/ratio : Height/ratio
            %                 for x=1 : BlockSize/ratio : Width/ratio
            %
            %                     startx = x;%max(x - floor(BlockSize/2/ratio), 1);
            %                     starty = y;%max(y - floor(BlockSize/2/ratio), 1);
            %                     endy = y+BlockSize/ratio-1;%min(y + floor(BlockSize/2/ratio), size(MS,1));
            %                     endx = x+BlockSize/ratio-1;%min(x + floor(BlockSize/2/ratio), size(MS,2));
            %                     MS_Band_d=MS_Band(3:ratio:end,3:ratio:end);
            %                     I_LP_Band_d=I_LP_Band(3:ratio:end,3:ratio:end);
            %                     BlockMS = MS_Band_d(starty:endy,startx:endx);
            %                     BlockPan = I_LP_Band_d(starty:endy,startx:endx);
            %
            %                     c = cov(BlockPan(:),BlockMS(:));
            %                     Coeff_Band((y-1)*ratio+1:(y-1)*ratio+S,...
            %                         (x-1)*ratio+1:(x-1)*ratio+S) = c(1,2)/var(BlockPan(:));
            %                 end
            %             end
        end
        
    else  %% Segmentation
        S_Band = squeeze(S(:,:,ii));
        labels = unique(S_Band);
        % PanSharpenedImage = zeros(size(imageLR));
        for i=1:length(labels)
            
            idx = S_Band==labels(i);
            
            %%% Coefficients
            c = cov(I_LP_Band(idx),MS_Band(idx));
            Coeff_Band(idx) = c(1,2)/var(I_LP_Band(idx));
            %%% Correlation
            %             Coeff_Band2(idx) = (I_LP_Band(idx)'*MS_Band(idx))/(I_LP_Band(idx)'*I_LP_Band(idx));
            
            %     % Final Mean Equalization  ??? To assess
            %     for ii = 1 : Bands
            %         PanSharpenedImage(idx,ii) = PanSharpenedImage(idx,ii) - mean(PanSharpenedImage(idx,ii)) + mean(imageLR(idx,ii));
            %     end
        end
    end
    
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
% for ii = 1 : size(I_MS,3)
%     h = PanSharpenedImage(:,:,ii);
%     PanSharpenedImage(:,:,ii) = h - mean2(h) + mean2(I_MS(:,:,ii));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  End  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

