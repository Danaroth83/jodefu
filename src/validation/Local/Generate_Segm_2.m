
%% Generate Segmentation maps
% nregs = [1 100 200 500 1000 2000];
nregs = [1 3 5 10 20 50 100 200 500 1000 2000 5000];% 10000 20000 50000 90000 179999];
%     PAN_S = zeros([size(PAN), 3, length(nregs)],'uint8');

k = 1;
[Sbpt_MS(:,:,k), T_MS] = segment_BPT(I_MS, nregs(k));
%     B = seg2bdry(Sbpt_MS(:,:,i),'imageSize');
%     PAN_S(:,:,:,i) = repmat(imadjust(datacast(PAN,'uint8')).*uint8(B==0),[1,1,3]);
%     PAN_S(:,:,1,i) = PAN_S(:,:,1,i) + 255.*uint8(B==1);

for k = 2:length(nregs)
    k
    Sbpt_MS(:,:,k) = segment_BPT(T_MS, nregs(k));
    %         B = seg2bdry(Sbpt_MS(:,:,i),'imageSize');
    %         PAN_S(:,:,:,i) = repmat(imadjust(datacast(PAN,'uint8')).*uint8(B==0),[1,1,3]);
    %         PAN_S(:,:,1,i) = PAN_S(:,:,1,i) + 255.*uint8(B==1);
end


%% Generate Segmentation maps PAN
% nregs = [1 100 200 500 1000 2000];
%     PAN_S = zeros([size(PAN), 3, length(nregs)],'uint8');
BlockSize=51;
SW=zeros(size(I_PAN));
for y=1 : size(I_PAN,1)
    for x=1 : size(I_PAN,2)
        startx = max(x - floor(BlockSize/2), 1);
        starty = max(y - floor(BlockSize/2), 1);
        endy = min(y + floor(BlockSize/2), size(I_MS,1));
        endx = min(x + floor(BlockSize/2), size(I_MS,2));
        BlockPAN = I_PAN(starty:endy,startx:endx);
        SW(y,x) = std2(BlockPAN);
    end
end
F2=cat(3,I_PAN/max(I_PAN(:)),SW/max(SW(:)));
k1 = 1;
[Sbpt_PAN(:,:,k1), T_PAN] = segment_BPT(F2, nregs(k1));
%     B = seg2bdry(Sbpt_MS(:,:,i),'imageSize');
%     PAN_S(:,:,:,i) = repmat(imadjust(datacast(I_PAN,'uint8')).*uint8(B==0),[1,1,3]);
%     PAN_S(:,:,1,i) = PAN_S(:,:,1,i) + 255.*uint8(B==1);

for k1 = 2:length(nregs)
    k1
    Sbpt_PAN(:,:,k1) = segment_BPT(T_PAN, nregs(k1));
    %         B = seg2bdry(Sbpt_MS(:,:,i),'imageSize');
    %         PAN_S(:,:,:,i) = repmat(imadjust(datacast(I_PAN,'uint8')).*uint8(B==0),[1,1,3]);
    %         PAN_S(:,:,1,i) = PAN_S(:,:,1,i) + 255.*uint8(B==1);
end
% break
%% Generate k-means segmentation maps on I_PAN
kregs = [1 3 5 10 20 50 100 200];
BlockSize=51;
SW=zeros(size(I_PAN));
for y=1 : size(I_PAN,1)
    for x=1 : size(I_PAN,2)
        startx = max(x - floor(BlockSize/2), 1);
        starty = max(y - floor(BlockSize/2), 1);
        endy = min(y + floor(BlockSize/2), size(I_MS,1));
        endx = min(x + floor(BlockSize/2), size(I_MS,2));
        BlockPAN = I_PAN(starty:endy,startx:endx);
        SW(y,x) = std2(BlockPAN);
    end
end
F=[I_PAN(:)/max(I_PAN(:)),SW(:)/max(SW(:))];
    
    
for kkk = 1:length(kregs)
    kkk
    IDX = kmeans(F,kregs(kkk));
    Skmeans_PAN(:,:,kkk) = reshape(IDX,[size(I_PAN,1) size(I_PAN,2)]);
end
    
%% Generate k-means segmentation maps on I_MS
% BlockSize=51;
% SW=zeros(size(I_ATWT));
% for y=1 : size(I_PAN,1)
%     for x=1 : size(I_PAN,2)
%         startx = max(x - floor(BlockSize/2), 1);
%         starty = max(y - floor(BlockSize/2), 1);
%         endy = min(y + floor(BlockSize/2), size(I_ATWT,1));
%         endx = min(x + floor(BlockSize/2), size(I_ATWT,2));
%         BlockPAN = I_PAN(starty:endy,startx:endx);
%         SW(y,x) = std2(BlockPAN);
%     end
% end
F1=[];
for ibands = 1 :size(I_ATWT,3)
    a=I_ATWT(:,:,ibands);
    %%%  Normalization
    %     F1=[F1, a(:)];
    F1=[F1, a(:)/max(a(:))];
    %[a(:)/max(a(:)),SW(:)/max(SW(:))];
end
% F1 = [F1,SW(:)/max(SW(:))];  %% Add Local I_PAN Intensity as a feature
kregs = [1 3 5 10 20];%50 100 200];
for kkk = 1:length(kregs)
    kkk
    IDX = kmeans(F1,kregs(kkk),'distance','cosine');
    Skmeans_MS_cosine(:,:,kkk) = reshape(IDX,[size(I_PAN,1) size(I_PAN,2)]);
end

% %%
% tic
% Spat_Kernel_vett=[100];
% Spect_Kernel_vett=[100];
% 
% for kkk = 1:length(Spat_Kernel_vett)
%     for kk=1:length(Spect_Kernel_vett)
%         hs=Spat_Kernel_vett(kkk)
%         hr=Spect_Kernel_vett(kk)
%     Smeanshift(:,:,kkk,kk) = meanshsegm(I_ATWT,hs,hr);
%     end
% end
% elapsedTime = toc

