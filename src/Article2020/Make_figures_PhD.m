clearvars; close all;

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'visualization'));
input_folder=fullfile(current_folder,'..','..','data','output','test_compression');
output_folder=fullfile(current_folder,'..','..','data','output','figures_PhD');

% im_tag='Beijing_WV3_WV3_4';
% im_tag='Hobart';
% im_tag='SanFrancisco_QB_QB_4';
% im_tag='Janeiro_cut256_4';
% im_tag='Janeiro_cut256_all';
% im_tag='Janeiro_cut256_RGB';
% im_tag='Stockholm_cut256_4';
% im_tag='Stockholm_cut256_all';
% im_tag='Stockholm_cut256_RGB';
im_tag='Washington_cut256_4';
% im_tag='Washington_cut256_RGB';

flag_tol=0; % if flag_tol==1, all the images are shown in the same reference
cd(fullfile(input_folder,im_tag));
filelist=dir('*.mat');
idx_best='q2^n';

output_folder=fullfile(output_folder,im_tag);
mkdir(output_folder);

for ii=1:numel(filelist)
    matfilename=filelist(ii).name;
    filename=matfilename(1:end-4);
    load(matfilename);
    if contains(filename,'default')
        I_GT=I_load{4};
        [a,tol]=viewimage_outputonly(I_GT.data(:,:,I_GT.Bands_to_display));
        figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[im_tag,'_GT.png']));
        if flag_tol==0, tol=[]; end       
        I_EXP=I_load{3};
        a=viewimage_outputonly(I_EXP.data(:,:,I_EXP.Bands_to_display),tol);
        figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[im_tag,'_EXP.png']));
        I_PAN=I_load{2};
        a=viewimage_outputonly(I_PAN.data(:,:,I_PAN.Bands_to_display));
        figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[im_tag,'_PAN.png']));
        break;
    end
end

for ii=1:numel(filelist)
    matfilename=filelist(ii).name;
    filename=matfilename(1:end-4);
    load(matfilename);
    
    if contains(filename,'compression')
        for jj=1:3
            if jj==1, test=[]; end
            if jj==2, MR=MR_radres; MF=MF_radres; test='_radres'; end
            if jj==3, MR=MR_jpeg; MF=MF_jpeg; test='_jpeg'; end
            for kk=1:numel(qindex_list), if strcmpi(qindex_list,idx_best), idx_cho=kk; break; end, end

            [~,idx_lambda]=max(MR(idx_cho,:));
            a=viewimage_outputonly(MF.data(:,:,MF.Bands_to_display,idx_lambda),tol);
            figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[filename,test,'_INVBEST.png']));
            
        end
    elseif contains(filename,'fusdem')
        for jj=1:numel(MR.qindex), if strcmpi(MR.qindex{jj},idx_best), idx_cho=jj; break; end, end
        idx_lambda=find(MR.bestindex(idx_cho,:)==1);
        idx_lambda2=find(MR.bestindex(idx_cho,:)==2);
        a=viewimage_outputonly(I_out.data(:,:,I_out.Bands_to_display,idx_lambda),tol);
        figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[filename,'_INVBEST.png']));
        a=viewimage_outputonly(I_out.data(:,:,I_out.Bands_to_display,idx_lambda),tol);
        figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[filename,'_INV2ND.png']));
    else
        for jj=1:numel(MR.qindex), if strcmpi(MR.qindex{jj},idx_best), idx_cho=jj; break; end, end
        idx_lambda=find(MR.bestindex(idx_cho,:)==1);
        %idx_lambda2=find(MR.bestindex(idx_cho,:)==2);
        a=viewimage_outputonly(I_out.data(:,:,I_out.Bands_to_display,idx_lambda),tol);
        figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[filename,'_INVBEST.png']));
        if ~contains(filename,'nomask')
            a=viewimage_outputonly(I_acq.data);
            figure; imshow(a,[]); imwrite(a,fullfile(output_folder,[filename,'_COMP.png']))
        end
    end
end

