addpath('segm\Meanshift')
% addpath('ps_util')
% addpath('Datasets\Local\Segmentations') 

%% Generate Segmentation maps
% nregs = [1 100 200 500 1000 2000];
% keyboard
% hs_vett=[5 10 20 50];
% hr_vett=[10 20 50 100 200];
hs_vett=[10 ];
hr_vett=[ 20 ];

kkk=0;

for is=1:length(hs_vett),
    for ir=1:length(hr_vett),
        kkk=kkk+1;
        Sbpt_MS_unordered(:,:,kkk)=meanshsegm(I_MS,hs_vett(is),hr_vett(ir));
        nregs_MS_unsorted(kkk)=max(max(Sbpt_MS_unordered(:,:,kkk)));
    end
end

[nregs_MS,I_sort]=sort(nregs_MS_unsorted);
Sbpt_MS=Sbpt_MS_unordered(:,:,I_sort);
