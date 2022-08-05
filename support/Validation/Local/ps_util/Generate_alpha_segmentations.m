
% alpha=reshape(alpha_opt_Hobart1_2_RR,[size
 
% addpath('Local\Tests_mfiles')
% addpath('Local\ps_util')
addpath('..\Datasets\Local\Segmentations')
% addpath('Datasets\Local\alpha_opt')
%% New estimation alpha coefficients
new_estimation=1;

%% Segmentation

if new_estimation || ~exist(['Segmentations_alpha_opt_',filename,'.mat'],'file')
%     cd Local
    addpath('segm')
    addpath('ps_util')
    [Alpha_opt,Alpha_opt_eq]=Alpha_opt_comput(I_MS,I_PAN,I_GT,ratio);
%     cd ..
else
    eval(['load Segmentations_',filename])
end


%% New Segmentation
new_segmentation=0;

%% Segmentation

if new_segmentation || ~exist(['Segmentations_',filename,'.mat'],'file')
%     cd Local
    addpath('segm')
    addpath('ps_util')
%     Generate_alpha_segmentations

%% Segmentation of coefficients (based on difference and not on SAM)
nregs = [1 3 5 10 20 50 100 200 500 1000 2000 5000];% 10000 20000 50000 90000 179999];

k = 1;
[Sbpt_MS_alpha(:,:,k), T_MS_alpha] = segment_BPT_sqe(Alpha_opt, nregs(k));
% [Sbpt_MS_alpha(:,:,k), T_MS_alpha] = segment_BPT_sqe(I_MS, nregs(k));
for k = 2:length(nregs)
    k
    Sbpt_MS_alpha(:,:,k) = segment_BPT_sqe(T_MS_alpha, nregs(k));
    %         B = seg2bdry(Sbpt_MS(:,:,i),'imageSize');
    %         PAN_S(:,:,:,i) = repmat(imadjust(datacast(PAN,'uint8')).*uint8(B==0),[1,1,3]);
    %         PAN_S(:,:,1,i) = PAN_S(:,:,1,i) + 255.*uint8(B==1);
end

kregs = [1 3 5 10 20 50 100 200];
F1=[];
for ibands = 1 :size(Alpha_opt,3)
    a=Alpha_opt(:,:,ibands);
    %%%  Normalization
    %     F1=[F1, a(:)];
    F1=[F1, a(:)/max(a(:))];
    %[a(:)/max(a(:)),SW(:)/max(SW(:))];
end
% F1 = [F1,SW(:)/max(SW(:))];  %% Add Local PAN Intensity as a feature
for kkk = 1:length(kregs)
    kkk
    IDX = kmeans(F1,kregs(kkk));
    Skmeans_MS_alpha(:,:,kkk) = reshape(IDX,[size(I_PAN,1) size(I_PAN,2)]);
end

% eval(['save Segmentations_alpha_opt_',filename,'.mat']
    cd ..
else
    eval(['load Segmentations_',filename])
end


break
   F1=[];
for ibands = 1 :size(I_MS,3)
    a=alpha_opt_Hobart1_2_RR(:,:,ibands);
    %%%  Normalization
        F1=[F1, a(:)];
%     F1=[F1, a(:)/max(a(:))];
end
for kkk = 1:length(kregs)
    kkk
    IDX = kmeans(F1,kregs(kkk));
   alpha_kmeans(:,:,kkk) = reshape(IDX,[size(I_PAN,1) size(I_PAN,2)]);
end