%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%           BDSD fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by
%           exploiting the Band-Dependent Spatial-Detail (BDSD) algorithm.
%
% Interface:
%           I_Fus_BDSD = BDSD(I_MS,I_PAN,ratio,S,sensor)
%
% Inputs:
%           I_MS:     MS image upsampled at PAN scale;
%           I_PAN:    PAN image;
%           ratio:    Scale ratio between MS and PAN. Pre-condition: Integer value;
%           sensor:   String for type of sensor (e.g. 'WV2', 'IKONOS')
%           S:        Segmentation (Image) or blocks (Scalar); see   below
%           interp:   Final spatial interpolation of coefficients
%           offset:   Use constant term in the expansion [Garzelli15]
%           ng_check: Check if expansion coefficient are positive and use global values

%
% Outputs:
%           I_Fus_BDSD:     BDSD pansharpened image.
%
% References:
%           [Garzelli08] A. Garzelli, F. Nencini, and L. Capobianco, “Optimal MMSE pan sharpening of very high resolution multispectral images,”
%                        IEEE Transactions on Geoscience and Remote Sensing, vol. 46, no. 1, pp. 228–236, January 2008.
%           [Garzelli15] A. Garzelli, “Pansharpening of Multispectral Images Based on Nonlocal Parameter Optimization,”
%                        IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 4, pp. 2096–2107, April 2015.
%           [Vivone14]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”,
%                        IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
% S can be:
%  - a scalar, in this case it is the block size and the conventional
%       algorithm will be taken (use the size of image to have 'global')
%  - a 1-band image, this is a segmentation map defining the regions in the
%       image and it will be considered for all the bands of MS
%  - a N-band image, this is a set of segmentation maps one associated to
%       each MS band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_Fus_BDSD = BDSD_local_reg(I_MS,I_PAN,ratio,sensor,S,offset,lambda,ng_check,interp)


addpath('ps_util')
%%%
% if nargin<6
%     offset=1;
%     ng_check=1;
%     lambda=1000;
% end
% if nargin<7
%     ng_check=1;
%     lambda=1000;
% end
% if nargin<8
%     ng_check=1;
% end
% 
% overlap = false;
% if nargin < 9
%     interp=0;
%     overlap=true;
%     fprintf('Overlapping ');
% end
overlap = false;
if isempty(interp)
    interp=0;
    overlap=true;
    fprintf('Overlapping ');
end

% Control of input parameters and initialization
I_MS = double(I_MS);
I_PAN = double(I_PAN);
[Height,Width,Bands] = size(I_MS);
BlockSize = [];
if isscalar(S)
    BlockSize = S;
    fprintf('Block size %d, interp %d\n', S, interp);
elseif (size(S,1) == size(I_PAN,1)) && (size(S,2) == size(I_PAN,2))
    %     if size(S,3) == 1
    %         S = repmat(S, [1, 1, Bands]);
    %     elseif size(S,3) ~= Bands
    %         error('Number of bands in S is incorrect.');
    %     end
    fprintf('Segmentation map\n');
else
    error('S is not correctly shaped.');
end
if isscalar(S)
    
    if overlap
        if S==1
            %%%% HPM method
            TBD
        else
            
            %%%% Global parameters
            % Reduced resolution
            pan_LP = MTF_PAN(I_PAN,sensor,ratio);
            %             pan_LP_d = pan_LP(3:ratio:end,3:ratio:end);
            pan_LP_d = imresize(pan_LP,1/ratio,'nearest');
            ms_orig = imresize(I_MS,1/ratio);
            %         ms_orig = MS(3:ratio:end,3:ratio:end,:);
            ms_LP_d = MTF(ms_orig,sensor,[],ratio);
            if offset
                %                 gamma_global=estimation_gains_alpha_offset(ms_LP_d,pan_LP_d,ms_orig,'global');
                gamma_global=estimation_gains_alpha_offset_reg(ms_LP_d,pan_LP_d,ms_orig,'global',lambda);
            else
                %                 gamma_global=estimation_gains_alpha(ms_LP_d,pan_LP_d,ms_orig,'global');
                gamma_global=estimation_gains_alpha_reg(ms_LP_d,pan_LP_d,ms_orig,'global',lambda);
            end
            %   To check: how to degrade within local window
            for y=1 : Height
                for x=1 : Width
                    startx = max(x - floor(BlockSize/2), 1);
                    starty = max(y - floor(BlockSize/2), 1);
                    endy = min(y + floor(BlockSize/2), size(I_MS,1));
                    endx = min(x + floor(BlockSize/2), size(I_MS,2));
                    BlockMS=I_MS(starty:endy,startx:endx,:);
                    %                     %     BlockMS_LP = imresize(BlockMS,1/ratio);
                    %                     BlockMS_LP = BlockMS(3:ratio:end,3:ratio:end,:);
                    %                     BlockMS_LP_d = MTF(BlockMS_LP,sensor,[],ratio);
                    BlockMS_LP = MTF(BlockMS_LP,sensor,[],ratio);
                    %                     BlockMS_LP_d = BlockMS(3:ratio:end,3:ratio:end,:);
                    BlockMS_LP_d = imresize(BlockMS_LP,1/ratio,'nearest');
                    
                    BlockPAN=I_PAN(starty:endy,startx:endx);
                    BlockPAN_LP = MTF_PAN(BlockPAN,sensor,ratio);
                    %                     BlockPAN_LP_d = BlockPAN_LP(3:ratio:end,3:ratio:end);
                    BlockPAN_LP_d = imresize(BlockPAN_LP,1/ratio,'nearest');
                    
                    if offset
                        %                         gamma=estimation_gains_alpha_offset(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global');
                        gamma=estimation_gains_alpha_offset_reg(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global',lambda);
                        
                        %%% Negative gamma check
                        if ng_check
                            if any(gamma(end-1,:)<0) %% Pan coefficients
                                gamma=gamma_global;
                            end
                        end
                        for kk = 1 : size(I_MS,3)
                            gamma_vec(1,1,:) = gamma(1:end,kk);
                            I_Fus_BDSD(y,x,kk) =...
                                I_MS(y,x,kk) + sum(gamma_vec.* cat(3,I_MS(y,x,:),I_PAN(y,x),1),3) ;
                        end
                    else
                        %   gamma=estimation_gains_alpha(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global');
                        gamma=estimation_gains_alpha_reg(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global',lambda);
                        %%% Negative gamma check
                        if ng_check
                            if any(gamma(end-1,:)<0) %% Pan coefficients
                                gamma=gamma_global;
                            end
                        end
                        for kk = 1 : size(I_MS,3)
                            gamma_vec(1,1,:) = gamma(1:end,kk);
                            I_Fus_BDSD(y,x,kk) =...
                                I_MS(y,x,kk) + sum(gamma_vec.* cat(3,I_MS(y,x,:),I_PAN(y,x)),3) ;
                        end
                    end
                    
                    %  I_Fus_BDSD_col(idx,:)=BlockMS+[BlockMS BlockPAN ones(size(BlockPAN,1),1)]*gamma;
                    
                    
                end
            end
        end
    else
        % Blockwise
        if (S > 1)
            if(rem(S,ratio))
                fprintf(1,'\n\n ');
                error('block size must be multiple of ratio')
            end
        end
        if(rem(Height,S)||rem(Width,S))
            fprintf(1,'\n\n ');
            keyboard
            error('x and y dims of pan must be multiple of the block size')
        end
        
        % Reduced resolution
        pan_LP = MTF_PAN(I_PAN,sensor,ratio);
        %         pan_LP_d = pan_LP(3:ratio:end,3:ratio:end);
        pan_LP_d = imresize(pan_LP,1/ratio,'nearest');
        ms_orig = imresize(I_MS,1/ratio);
        %         ms_orig = MS(3:ratio:end,3:ratio:end,:);
        ms_LP_d = MTF(ms_orig,sensor,[],ratio);
        
        %%%% Global parameters
        if offset
            %             gamma_global=estimation_gains_alpha_offset(ms_LP_d,pan_LP_d,ms_orig,'global');
            gamma_global=estimation_gains_alpha_offset_reg(ms_LP_d,pan_LP_d,ms_orig,'global',lambda);
        else
            %             gamma_global=estimation_gains_alpha(ms_LP_d,pan_LP_d,ms_orig,'global');
            gamma_global=estimation_gains_alpha_reg(ms_LP_d,pan_LP_d,ms_orig,'global',lambda);
        end
        for y=1 : S/ratio : ceil(Height/ratio)
            for x=1 : S/ratio : ceil(Width/ratio)
                startx = x ;
                starty = y;
                endy = y+S/ratio-1;
                endx = x+S/ratio-1;
                BlockMS_LP_d = ms_LP_d(starty:endy,startx:endx,:);
                BlockPAN_LP_d = pan_LP_d(starty:endy,startx:endx);
                BlockMS_LP = ms_orig(starty:endy,startx:endx,:);
                BlockMS = I_MS((y-1)*ratio+1:(y-1)*ratio+S,(x-1)*ratio+1:(x-1)*ratio+S,:);
                BlockPAN = I_PAN((y-1)*ratio+1:(y-1)*ratio+S,(x-1)*ratio+1:(x-1)*ratio+S);
                if offset
                    %                     gamma=estimation_gains_alpha_offset(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global');
                    gamma=estimation_gains_alpha_offset_reg(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global',lambda);
                else
                    %                     gamma=estimation_gains_alpha(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global');
                    gamma=estimation_gains_alpha_reg(BlockMS_LP_d,BlockPAN_LP_d,BlockMS_LP,'global',lambda);
                end
                %%% Negative gamma check
                if ng_check
                    if any(gamma(end-1,:)<0) %% Pan coefficients
                        gamma=gamma_global;
                    end
                end
                %                 keyboard
                for kk = 1 : size(I_MS,3)
                    gamma_vec(1,1,:) = gamma(1:end,kk);
                    if offset
                        I_Fus_BDSD((y-1)*ratio+1:(y-1)*ratio+S,(x-1)*ratio+1:(x-1)*ratio+S,kk) =...
                            BlockMS(:,:,kk) + sum(repmat(gamma_vec,[size(BlockMS,1) size(BlockMS,2) 1])...
                            .* cat(3,BlockMS,BlockPAN,ones(size(BlockPAN))),3) ;
                    else
                        I_Fus_BDSD((y-1)*ratio+1:(y-1)*ratio+S,(x-1)*ratio+1:(x-1)*ratio+S,kk) =...
                            BlockMS(:,:,kk) + sum(repmat(gamma_vec,[size(BlockMS,1) size(BlockMS,2) 1])...
                            .* cat(3,BlockMS,BlockPAN),3) ;
                    end
                end
            end
        end
    end
else
    %     S_R=S(3:ratio:end,3:ratio:end);
    S_R = imresize(S,1/ratio,'nearest');
    pan_LP = MTF_PAN(I_PAN,sensor,ratio);
    %         pan_LP_d = pan_LP(3:ratio:end,3:ratio:end);
    pan_LP_d = imresize(pan_LP,1/ratio,'nearest');
    %     keyboard
    ms_orig = imresize(I_MS,1/ratio);
    ms_LP_d = MTF(ms_orig,sensor,[],ratio);
    
    if offset
        gamma_global=estimation_gains_alpha_offset_reg(ms_LP_d,pan_LP_d,ms_orig,'global',lambda);
    else
        gamma_global=estimation_gains_alpha_reg(ms_LP_d,pan_LP_d,ms_orig,'global',lambda);
    end
    %%%% Local parameters
    labels = unique(S);
    I_Fus_BDSD_col=zeros(Width*Height,Bands);
    for ii=1:length(labels)
        idx = S==labels(ii);
        idx_R = S_R==labels(ii);
        for kk=1:Bands
            b = ms_LP_d(:,:,kk);
            BlockMS_LP_d(:,kk) = b(idx_R);
            c = ms_orig(:,:,kk);
            BlockMS_LP(:,kk) = c(idx_R);
            d = I_MS(:,:,kk);
            BlockMS(:,kk) = d(idx);
        end
        BlockPAN_LP_d = pan_LP_d(idx_R);
        BlockPAN=I_PAN(idx);
        Diff = BlockMS_LP-BlockMS_LP_d;
        if offset
            H= [BlockMS_LP_d, BlockPAN_LP_d, ones(size(BlockPAN_LP_d,1),1)];
%             keyboard
        else
            H= [BlockMS_LP_d, BlockPAN_LP_d];
        end
        gamma= (H+lambda*eye(size(H)))\Diff;
        %%% Negative gamma check
        if ng_check
            if any(gamma(end-1,:)<0) %% Pan coefficients
                gamma=gamma_global;
            end
        end
        %     % Final Mean Equalization  ??? To assess
        %     for ii = 1 : Bands
        %         PanSharpenedImage(idx,ii) = PanSharpenedImage(idx,ii) - mean(PanSharpenedImage(idx,ii)) + mean(imageLR(idx,ii));
        %     end
        if offset
            I_Fus_BDSD_col(idx,:)=BlockMS+[BlockMS BlockPAN ones(size(BlockPAN,1),1)]*gamma;
        else
            I_Fus_BDSD_col(idx,:)=BlockMS+[BlockMS BlockPAN]*gamma;
        end
        clear BlockMS_LP_d BlockMS_LP BlockMS BlockPAN
    end
    I_Fus_BDSD=reshape(I_Fus_BDSD_col,[Height,Width,Bands]);
end


