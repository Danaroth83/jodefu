function [PanSharpenedImage, Coeff] = PS_local(I_MS,I_PAN,PAN_LP,mode,S,Threshold,equaliz, interp)
% mode: modality for computing the coefficients for weighting the details
% to inject
%
% S can be:
%  - a scalar, in this case it is the block size and the conventional
%       algorithm will be taken
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

switch mode
    case 'CBD'
        fprintf('CBD [Aiazzi02] Th:%3.2f\t', Threshold);
    case 'ECB'
        fprintf('ECB [Aiazzi06] Th:%3.2f\t', Threshold);
end



overlap = false;
if isempty(interp)
    interp=0;
    overlap=true;
    fprintf('Overlapping ');
end

I_MS = double(I_MS);
[Height,Width,Bands] = size(I_MS);

% generate the details
I_PAN = repmat(double(I_PAN), [1, 1, Bands]);
PAN_LP = double(PAN_LP);
if size(PAN_LP, 3) == 1
    PAN_LP = repmat(PAN_LP, [1, 1, Bands]);
elseif size(PAN_LP, 3) ~= size(I_PAN, 3)
    error('PAN_LP should have the same number of bands as PAN');
end

% Equalization PAN with I_LP
% keyboard
switch equaliz
    case 1
        for ii = 1:Bands,
            imageHR(:,:,ii) = I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii))+ mean2(PAN_LP(:,:,ii));
        end
    case 2
        for ii = 1:Bands,
            imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                std2(PAN_LP(:,:,ii))/std2(I_PAN(:,:,ii)) + mean2(PAN_LP(:,:,ii));
        end
    case 3
        for ii = 1:Bands,
            imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                std2(PAN_LP(:,:,ii))/std2(imresize(I_PAN(:,:,ii),1/ratio)) +...
                mean2(PAN_LP(:,:,ii));
        end
    case 4
        for ii = 1:Bands,
            imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                std2(PAN_LP(:,:,ii))/std2(imresize(imresize(I_PAN(:,:,ii),1/ratio),ratio)) +...
                mean2(PAN_LP(:,:,ii));
        end
    otherwise
        imageHR = I_PAN;
end
% keyboard


% DetailsHRPan = PAN - PAN_LP;
DetailsHRPan = imageHR - PAN_LP;

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

for b=1:Bands
    Global_Correlation = corr2(I_MS(:,:,b),PAN_LP(:,:,b));
    if ~isempty(BlockSize)  % regular algorithm
        if overlap
            for y=1 : Height
                for x=1 : Width
                    
                    startx = max(x - floor(BlockSize/2), 1);
                    starty = max(y - floor(BlockSize/2), 1);
                    endy = min(y + floor(BlockSize/2), size(I_MS,1));
                    endx = min(x + floor(BlockSize/2), size(I_MS,2));
                    
                    BlockMS = I_MS((starty:endy),(startx:endx),b);
                    BlockPan = PAN_LP((starty:endy),(startx:endx),b);
                    StdMS=std2(BlockMS);
                    StdPan=std2(BlockPan);
                    Correlation=corr2(BlockMS,BlockPan);
                    
                    %%%%
                    switch mode
                        case 'CBD'
                            if(Correlation>Threshold)
                                LG=StdMS/StdPan;
                            else
                                LG=0;
                            end
                        case 'ECB'
                            LG = min((Correlation./Global_Correlation) .* (StdMS/StdPan),Threshold);
                    end
                    
                    Coeff(y:endy,x:endx,b) = LG;
                end
            end
        else
            for y=1 : BlockSize : Height+BlockSize
                for x=1 : BlockSize : Width+BlockSize
                    
                    startx = max(x - floor(BlockSize/2), 1);
                    starty = max(y - floor(BlockSize/2), 1);
                    endy = min(y + floor(BlockSize/2), size(I_MS,1));
                    endx = min(x + floor(BlockSize/2), size(I_MS,2));
                    
                    %                     endy = y + BlockSize;
                    %                     endx = x + BlockSize;
                    %
                    %                     if endy > Height
                    %                         endy = Height;
                    %                     end
                    %
                    %                     if endx > Width
                    %                         endx = Width;
                    %                     end
                    
                    BlockMS = I_MS((starty:endy),(startx:endx),b);
                    BlockPan = PAN_LP((starty:endy),(startx:endx),b);
                    %                     BlockPan = PAN((starty:endy),(startx:endx),b);
                    StdMS=std2(BlockMS);
                    StdPan=std2(BlockPan);
                    Correlation=corr2(BlockMS,BlockPan);
                    
                    switch mode
                        case 'CBD'
                            if(Correlation>Threshold)
                                LG=StdMS/StdPan;
                            else
                                LG=0;
                            end
                        case 'ECB'
                            LG = min((Correlation./Global_Correlation) .* (StdMS/StdPan),Threshold);
                    end
                    
                    Coeff(starty:endy,startx:endx,b) = LG;
                end
            end
        end
    else    % local algo based on the segmentation map
        S_b = S(:,:,b);
        MS_b = I_MS(:,:,b);
        Coeff_b = zeros(size(MS_b));
        
        labels = unique(S_b);
        nregs = length(labels);
        
        %         h = waitbar(0,['Process band ', int2str(b), '/', int2str(Bands)]);
        for i=1:nregs
            
            %             if mod(i, 100) == 0
            %                 waitbar(nregs/i,h);
            %             end
            
            BlockMS = MS_b(S_b==labels(i));
            %             BlockPan = PAN(S_b==labels(i));
            BlockPan = PAN_LP(S_b==labels(i));
            StdMS=std2(BlockMS);
            StdPan=std2(BlockPan);
            Correlation=corr2(BlockMS,BlockPan);
            
            switch mode
                case 'CBD'
                    if(Correlation>Threshold)
                        LG=StdMS/StdPan;
                    else
                        LG=0;
                    end
                case 'ECB'
                    LG = min(abs((Correlation./Global_Correlation)) .* (StdMS/StdPan),Threshold);
            end
            Coeff_b(S_b==labels(i)) = LG;
        end
        %         close(h);
        Coeff(:,:,b) = Coeff_b;
    end
end

% interpolation of the coefficients
if ~isempty(BlockSize) && ~overlap && interp ~= 0
    
    naddrows = (ceil(size(Coeff,1)/BlockSize)+1) * BlockSize - size(Coeff,1);
    naddcols = (ceil(size(Coeff,2)/BlockSize)+1) * BlockSize - size(Coeff,2);
    Coeff1 =zeros([size(Coeff,1)+naddrows, size(Coeff,2)+naddcols, Bands]);
    Coeff1(1:size(Coeff,1),1:size(Coeff,2),:) = Coeff;
    Coeff1(size(Coeff,1)+1:end,:,:) = repmat(Coeff1(size(Coeff,1),:,:), [naddrows, 1, 1]);
    Coeff1(:,size(Coeff,2)+1:end,:) = repmat(Coeff1(:,size(Coeff,2),:), [1, naddcols, 1]);
    Coeff1(size(Coeff,1)+1:end,size(Coeff,2)+1:end,:) = repmat(Coeff1(size(Coeff,1),size(Coeff,2),:), [naddrows, naddcols, 1]);
    Coeff1 = imresize(Coeff1, 1/BlockSize, 'nearest');
    %     Coeff2 = imresize(Coeff1, BlockSize, 'nearest');% for check
    
    if interp==1
        Coeff = imresize(Coeff1, BlockSize, 'bilinear');
    else
        Coeff = imresize(Coeff1, BlockSize, 'bicubic');
    end
    clear Coeff1
    Coeff = Coeff(naddrows:naddrows+Height-1,naddrows:naddcols+Width-1,:);
end

PanSharpenedImage = Coeff .* DetailsHRPan + I_MS;

end