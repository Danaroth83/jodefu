function [I_LP, alpha] = genI_LP(MS, mode_ilp, PAN_LP, varargin)
% Create a low-pass version of the PAN obtained as a composition of the MS
% bands.
%
% Examples
% mode can be:
% - 'mean':  I_LP is the mean of the MS bands (GS mode 1)
% - 'MRA' :  I_LP is equal to PAN_LP obtained outside (GS mode 2)
% - 'blocks'  : I_LP is the linear combination of MS bands over blocks
% - 'adaptive': I_LP is the linear combination of MS bands over
%               segmentation regions

if nargin<2
    mode_ilp = 'mean';
end

fprintf('Estimation of I_LP by %s\n', mode_ilp);

%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------

S = [];
non_neg = false;
BlockSize = [];
overlap = [];
MS_LP = MS;
[Height,Width,Bands] = size(MS);
[Height_PAN,Width_PAN,Bands_PAN] = size(PAN_LP);
% keyboard
%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SEGM'
                S = varargin{i+1};
                fprintf(' local S(%d regs)',length(unique(S)));
            case 'NON_NEG'
                non_neg = varargin{i+1};
                fprintf(' non_neg %d',non_neg);
            case 'MS_LP'
                MS_LP = varargin{i+1};
                tmp = PAN_LP;
%                 clear PAN_LP
%                 tmp2 = S;
%                 clear S
%                 S = imresize(tmp2,[size(MS_LP,1),size(MS_LP,2)],'nearest');
                PAN_LP = imresize(tmp,[size(MS_LP,1),size(MS_LP,2)],'nearest');
                fprintf(' MS_LP');
                if size(MS_LP,1) ~= size(PAN_LP,1) || size(MS_LP,2) ~= size(PAN_LP,2)
                    error('The size of MS_LP and PAN must be the same');
                end
            case 'BLOCKSIZE'
                BlockSize = varargin{i+1};
                case 'OVERLAP'
                overlap = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end
fprintf('\n');

if ~isempty(S) && (size(MS_LP,1) ~= size(S,1) || size(MS_LP,2) ~= size(S,2))
    error('The size of MS_LP and the segmentation map must be the same');
end

if isempty(BlockSize) && strcmp(mode_ilp,'blocks')
    BlockSize=55;  %Default value
end
if isempty(overlap) && strcmp(mode_ilp,'blocks')
    overlap=0;%Default value
end
    
switch mode_ilp
    case 'mean' % For GS mode 1
        %         fprintf('GLP mode 1 [Aiazzi07]\n');
        I_LP = mean(MS,3);
    case 'MRA' % Multi Resolution Analysis for GS mode 2
        I_LP = PAN_LP;
    
    case 'blocks' % 
%         keyboard
        fprintf(' Blocks %d x %d\n',BlockSize,BlockSize);
        for y=1 : BlockSize : Height+BlockSize
            for x=1 : BlockSize : Width+BlockSize
                
                startx = max(x - floor(BlockSize/2), 1);
                starty = max(y - floor(BlockSize/2), 1);
                endy = min(y + floor(BlockSize/2), size(MS,1));
                endx = min(x + floor(BlockSize/2), size(MS,2));
                BlockMS = MS_LP(starty:endy,startx:endx,:);
                BlockPan = PAN_LP(starty:endy,startx:endx);
                BlockMS_col = [reshape(double(BlockMS), size(BlockMS,1)*size(BlockMS,2), Bands)...
                    ones(size(BlockMS,1)*size(BlockMS,2),1)];
                if non_neg
                    alpha_loc = lsqnonneg(BlockMS_col,BlockPan(:));
                else
                    alpha_loc = BlockMS_col\BlockPan(:);
                    %                 regress(reshape(imageHR0,[size(imageHR0,1)*size(imageHR0,2) 1]),reshape(cat(3,imageLR_LP0,ones(size(I_MS_LP,1),size(I_MS_LP,2))),[size(imageLR_LP0,1)*size(imageLR_LP0,2) size(imageLR_LP0,3)+1]));
                end
                alpha2(1,1,:)= alpha_loc';
                alpha(starty:endy,startx:endx,:) = repmat(alpha2, [size(BlockMS,1),size(BlockMS,2),1]);
            end
        end
%         I_LP = sum([MS, ones(size(MS,1),1)] .* alpha, 2);
%         I_LP = reshape(I_LP, Height, Width);
        I_LP = sum(cat(3,MS,ones(size(MS,1),size(MS,2))).*alpha,3);

    case 'adaptive' % GSA adaptive
        % fprintf('GLP adaptative [Aiazzi07]\n');
        
        %%% Global ([Aiazzi07]): Also achievable through 1 class segmentation
        if isempty(S)  
            S = ones(size(MS,1), size(MS,2));
        end
        
        MS = reshape(double(MS), Height*Width, Bands);
        PAN_LP = double(PAN_LP(:));
%         PAN_LP = reshape(double(PAN_LP), size(PAN_LP,1)*size(PAN_LP,2), Bands_PAN);
        MS_LP = reshape(double(MS_LP), size(MS_LP,1)*size(MS_LP,2), Bands);
        
        labels = unique(S);
        if length(labels)>1
            alpha = zeros(size(MS,1),size(MS,2)+1);
            for i=1:length(labels)
                MS_LP_loc = [MS_LP(S==labels(i),:), ones(nnz(S==labels(i)),1)];
                if non_neg
                    alpha_loc = lsqnonneg(MS_LP_loc,PAN_LP(S==labels(i)));
                else
                    alpha_loc = MS_LP_loc\PAN_LP(S==labels(i));
                    %                 regress(reshape(imageHR0,[size(imageHR0,1)*size(imageHR0,2) 1]),reshape(cat(3,imageLR_LP0,ones(size(I_MS_LP,1),size(I_MS_LP,2))),[size(imageLR_LP0,1)*size(imageLR_LP0,2) size(imageLR_LP0,3)+1]));
                end
                alpha(S==labels(i),:) = repmat(alpha_loc', [nnz(S==labels(i)),1]);
            end
        else
%             keyboard
            MS_LP = [MS_LP, ones(size(MS_LP,1),1)];
            if non_neg
                alpha = lsqnonneg(MS_LP,PAN_LP);
            else
%                 keyboard
                alpha = MS_LP\PAN_LP;
                %                 regress(reshape(imageHR0,[size(imageHR0,1)*size(imageHR0,2) 1]),reshape(cat(3,imageLR_LP0,ones(size(I_MS_LP,1),size(I_MS_LP,2))),[size(imageLR_LP0,1)*size(imageLR_LP0,2) size(imageLR_LP0,3)+1]));
            end
            squeeze(alpha)';
            alpha = repmat(alpha', [size(MS,1),1]);
        end
        
        
        I_LP = sum([MS, ones(size(MS,1),1)] .* alpha, 2);
        I_LP = reshape(I_LP, Height, Width);
end

