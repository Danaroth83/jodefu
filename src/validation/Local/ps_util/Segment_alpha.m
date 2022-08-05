addpath('Local\Tests_mfiles')
addpath('Local\Segmentations')

    addpath('segm')
Alpha_opt=Alpha_opt_comput(I_MS,I_PAN,I_GT,ratio);
F1=[];
for ibands = 1 :size(I_MS,3)
    a=Alpha_opt(:,:,ibands);
    %%%  Normalization
        F1=[F1, a(:)];
%     F1=[F1, a(:)/max(a(:))];
    %[a(:)/max(a(:)),SW(:)/max(SW(:))];
end
kregions=[ 1 2 5 10 20];
for kkk = 1:length(kregions)
    kkk
    IDX = kmeans(F1,kregions(kkk));
    Skmeans_alpha(:,:,kkk) = reshape(IDX,[size(I_PAN,1) size(I_PAN,2)]);
end
nregions = [1 100 200 500 1000];
k1 = 1;
[Sbpt_alpha(:,:,k1), T_alpha] = segment_BPT(Alpha_opt, nregions(1));
%     B = seg2bdry(Sbpt_MS(:,:,i),'imageSize');
%     PAN_S(:,:,:,i) = repmat(imadjust(datacast(I_PAN,'uint8')).*uint8(B==0),[1,1,3]);
%     PAN_S(:,:,1,i) = PAN_S(:,:,1,i) + 255.*uint8(B==1);

for k1 = 2:length(nregions)
    k1
    Sbpt_alpha(:,:,k1) = segment_BPT(T_alpha, nregions(k1));
    %         B = seg2bdry(Sbpt_MS(:,:,i),'imageSize');
    %         PAN_S(:,:,:,i) = repmat(imadjust(datacast(I_PAN,'uint8')).*uint8(B==0),[1,1,3]);
    %         PAN_S(:,:,1,i) = PAN_S(:,:,1,i) + 255.*uint8(B==1);
end