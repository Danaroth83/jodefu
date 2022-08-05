%% View All


%  flagPAN:            Flag. If flagPAN == 1, the first image to plot is the panchromatic image otherwise it is the ground-truth.
if strcmp(filename(end-1),'R')
    flag_PAN=0;
    MatrixImage(:,:,:,1) = HRMS;
    
    flag_cut_bounds = 1;
    dim_cut = 11;
else
    flag_PAN=1;
    MatrixImage(:,:,:,1) = repmat(PAN,[1 1 size(MS,3)]);
    %% Cut Final Image
    flag_cut_bounds = 0;
    dim_cut = 0;
end
if size(MS,3) == 4
    vect_index_RGB = [3,2,1];
else
    vect_index_RGB = [5,3,1];
end

%% % Compare Segmentations
I_LP_input = genPAN_LP(PAN, 'ATWT', Resize_fact);
equaliz = 0;
interp = 0;
MatrixImage(:,:,:,2) = MS;
MatrixImage(:,:,:,3) = GS_local(MS, PAN, I_LP_input,Sbpt_MS(:,:,1),equaliz,interp);
MatrixImage(:,:,:,4) = GS_local(MS, PAN, I_LP_input,Sbpt_MS(:,:,5),equaliz,interp);
MatrixImage(:,:,:,5) = GS_local(MS, PAN, I_LP_input,Sbpt_MS(:,:,7),equaliz,interp);
MatrixImage(:,:,:,6) = GS_local(MS, PAN, I_LP_input,Sbpt_MS(:,:,10),equaliz,interp);
MatrixImage(:,:,:,7) = GS_local(MS, PAN, I_LP_input,128,equaliz,interp);
MatrixImage(:,:,:,8) = GS_local(MS, PAN, I_LP_input,64,equaliz,interp);
MatrixImage(:,:,:,9) = GS_local(MS, PAN, I_LP_input,16,equaliz,interp);

titleImages = {'PAN','MS',['BPT Segm ',num2str(nregs(1))],['BPT Segm ',num2str(nregs(5))],...
    ['BPT Segm ',num2str(nregs(7))],['BPT Segm ',num2str(nregs(10))],...
    'Blocks 128 \times 128','Blocks 64 \times 64','Blocks 16 \times 16'};




figure,
MatrixPrint = showImagesAll(MatrixImage,titleImages,vect_index_RGB,flag_cut_bounds,dim_cut,flag_PAN);


for ii = 1 : size(MatrixPrint,4)
    printImage(MatrixPrint(:,:,:,ii),sprintf('Outputs/%d.eps',ii));
end


