%% FUSION ALGORITHMS (PANSHARPENING)
%
% Description:
%
% This function realizes the fusion among a low resolution (LR) image with 
% high spectral diversity and a monochromatic high resolution (HR) image
%
% Interface:
%
% MatrixImage=Fusion_Algorithms_Pansharpening(MS,PAN,'Field',Value)
% eg: MI=Fusion_Algorithms_Pansharpening(MS,PAN,'methods_list',{'EXP','GSA','PCA'});
%
% Inputs:
%
% MS: Struct representing the LR image (if interpolated, its label field has to be 'EXP')
% PAN: Struct representing the HR image
%     The two structs above should contain the fields:
%         data: the image itself
%         sensor: String describing the sensor
%         KerBlu0: The degradation matrix (assumption: same for all bands)
%
% Optional Fields:
%
% 'EXP': LR image upsampled to the HR scale
% 'methods_list': cell containing the strings for the methods to test
% 'ratio_tgt_REF': scale ratio of the fused image, if different from ratio between HR and LR
% 'flag_degradation': Method to for generation of Low Pass version of the HR image 
%       (0= no degradation; 1= with bicubic filter; 2=with  HR MTF-matched
%       filter; 3=with Wavelet; 4=with Hamming windowed sinc; 5=with LR MTF-matched filter)
% 'flag_equalization': Equalization method for equalization of low pass HR image
%       (0=classic; 1=with GT; 2=regression based; 3=PCA regression based)
% 'GT': Ground Truth (Should not in general be necessary for fusion)
% 'printMAT': if 1, prints the inputs and the outputs on a .mat file
% 'output_file': name of output file mat
%
% Available Fusion Methods (accepted strings in methods_list):
%
% 'EXP': Simple interpolation of the LR image;
% 'PCA': Principal Component Analysis;
% 'IHS': Generalized Intensity Hue Saturation;
% 'Brovey': Brovey Transform;
% 'BDSD': Band-Dependent Spatial Detail;
% 'G_BDSD': Band-Dependent Spatial Detail (Global);
% 'L_BDSD': Band-Dependent Spatial Detail (Local);
% 'C_BDSD': Band-Dependent Spatial Detail (with Segmentation);
% 'GS': Gram–Schmidt decomposition;
% 'GSA': Gram–Schmidt Adaptive;
% 'PRACS': Partial Replacement Adaptive Component Substitution;
% 'HPF': High Pass Filtering;
% 'SFIM': Smoothing-Filter-based Intensity Modulation;
% 'Indusion': Indusion;
% 'ATWT': A Trous Wavelet Transform;
% 'ATWT_mod': A Trous Wavelet Transform (modified version);
% 'AWLP': Additive Wavelet Luminance Proportional;
% 'ATWT_M2': A Trous Wavelet Transform - Model 2;
% 'ATWT_M2': A Trous Wavelet Transform - Model 3;
% 'MTF_GLP': Generalized Laplacian Pyramid with MTF-matched filter;
% 'MTF_GLP_CBD': 'MTF_GLP' with Context Based Decision injection;
% 'MTF_GLP_HPM': 'MTF_GLP' with High Pass Modulation injection;
% 'MTF_GLP': Generalized Laplacian Pyramid with MTF-matched filter;
% 'MTF_GLP_HPM_PP': MTF-GLP-HPM with Post Processing;
% 'MF-HG': Morphological Half Gradient
% 'GFPCA': Guided Filter in the PCA Domain
% 'CNMF': Convolutive Nonnegative Matrix Factorization;
% 'HySure': Hyperspectral Superresolution;
% 'BayesNaive': Bayesian inversion with Naive a priori;
% 'BayesSparse': Bayesian inversion with Sparse a priori;
% 'BayesML': Bayesian inversion with ML a priori;
% 'BayesNLM: Bayesian inversion with NLM a priori;
% 'BayesNaive_FE': 'BayesNaive' with degradation filter estimation;
% 'BayesSparse_FE': 'BayesSparse' with degradation filter estimation;
%
% Output:
%
% MatrixImage: struct containing a stack of fused images whose fields are:
%     - data: the stack of fused images themselves
%     - methods_list: a descriptor of the fusion methods
%     - titleImages: list of labels for each fusion method
%     - time: vector of times used for each fusion
%     - ratio: Scale ratio between the fused stack and the LR image

function Fused=Fusion_Algorithms_Pansharpening(Image_MS,Image_PAN,varargin)

%% Definition of folders

original_folder=pwd;
current_folder=fileparts(mfilename('fullpath'));
support_folder=fullfile(current_folder,'..');
fusion_folder=fullfile(support_folder,'Fusion');
output_folder=fullfile(current_folder,'..','..','data','output');

addpath(fullfile(support_folder,'Load_info'),...
        fullfile(support_folder,'Scaling'),...
        fullfile(support_folder,'Filtering'));

%% Loading input image and its metadata

I_MS_LR=[];  % If 1, skip interpolation of MS_LR
I_MS=[];  % If 1, skip downsampling of MS
I_GT=[];

if isfield(Image_PAN,'data'), I_PAN=Image_PAN.data; else, I_PAN=Image_PAN; end
if isfield(Image_MS,'label')
    if strcmpi(Image_MS.label,'EXP')
        I_MS=Image_MS.data;
        if isempty(I_PAN), ratio=1;
        elseif isfield(Image_PAN,'ratio'), ratio=Image_PAN.ratio;
        elseif isfield(Image_MS,'ratio'), ratio=Image_MS.ratio;
        else, ratio=4; disp('Warning: Fusion will assume scale ratio of 4');
        end
        if isfield(Image_MS,'time'), time_EXP=Image_MS.time; else, time_EXP=0; end
    else
        I_MS_LR=Image_MS.data;
        if isfield(Image_MS,'edge_cut'), edge_cut=Image_MS.edge_cut; else, edge_cut=0; end
        if isempty(I_PAN), ratio=1;
        else, ratio=size(I_PAN,1)/(size(I_MS_LR,1)-2*edge_cut);
        end
    end
else
    if isfield(Image_MS,'data'), I_MS_LR=Image_MS.data; else, I_MS_LR=Image_MS; end
    if isempty(I_PAN), ratio=1;
    else, ratio=size(I_PAN,1)/size(I_MS_LR,1);
    end
end

if isfield(Image_MS,'im_tag'), im_tag=Image_MS.im_tag; else, im_tag=[]; end
if isfield(Image_MS,'sensor'), sensor=Image_MS.sensor; else, sensor='none'; end
if isfield(Image_PAN,'sensor'), sensor_PAN=Image_PAN.sensor; else, sensor_PAN='none'; end
if isfield(Image_MS,'type'), type=Image_MS.type; else, type='MS'; end
if isfield(Image_PAN,'type'), type_PAN=Image_PAN.type; else, type_PAN='PAN'; end
if isfield(Image_MS,'Bands_selected'), Bands_to_sharpen=Image_MS.Bands_selected; else, Bands_to_sharpen=1:size(I_MS_LR,3); end
if isfield(Image_MS,'GNyq'), GNyq=Image_MS.GNyq; else, GNyq=load_MTF(sensor,Bands_to_sharpen,im_tag,type); end
% if isfield(Image_MS,'KerBlu'), KerBlu=Image_MS.KerBlu; else, KerBlu=load_KerBlu(im_tag,ratio,GNyq,min(size(I_PAN,1),size(I_PAN,2)),false); end
if isfield(Image_PAN,'GNyq'), GNyq_PAN=Image_PAN.GNyq; else, GNyq_PAN=load_MTF_PAN(sensor_PAN); end
if isfield(Image_MS,'IntBits'), L=Image_MS.IntBits; else, L=load_resolution(sensor, im_tag, type); end
if isfield(Image_MS,'Band_overlap_PAN'), Band_overlap_PAN=Image_MS.Band_overlap_PAN; else, Band_overlap_PAN=load_Band_overlap(sensor,sensor_PAN,type,type_PAN); end

[KerBlu0,start_pos]=load_KerBlu(im_tag,ratio,GNyq,min(size(I_PAN,1),size(I_PAN,2)));
% KerBlu=load_KerBlu(im_tag,ratio,GNyq,min(size(I_PAN,1),size(I_PAN,2)),false);

%% Loading Parameters for the Fusion

methods_list=load_methodslist('originalToolbox');
flag_thvalues=0;
flag_FR=0;
flag_degradation=0;
flag_equalization=0;
ratio_tgt_REF=ratio;
flag_interpolation=0;
flag_imresize_RR=0;

printMAT=0;
output_file=fullfile('Pansharpening','Output.mat');

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if any(strcmpi(pname,{'I_MS','MS','EXP','MS-LR','MS_LR'}))
        if isfield(pval,'data')
            if isfield(pval,'label') && strcmpi(pval.label,'MS-LR')
                I_MS_LR=pval.data;
                if isfield(pval,'edge_cut'), edge_cut=pval.edge_cut; else, edge_cut=0; end
            else
                I_MS=pval.data;
            end
        else
            I_MS=pval;
        end
        if isfield(pval,'time'), time_EXP=pval.time; else, time_EXP=0; end
    elseif any(strcmpi(pname,{'I_GT','GT'}))          % Allows passage of Ground Truth for extra operations (this is not allowed in Full Resolution)
        if isfield(pval,'data'), I_GT=pval.data; else, I_GT=pval; end
        flag_FR=1;
    elseif strcmpi(pname,'flag_thvalues')             % Thresholding the image within the allowed range after fusion
        flag_thvalues=pval;
    elseif strcmpi(pname,'flag_degradation')          % Choice on how to generate the low resolution equivalent of the HR image (0, No degradation; 1, with imresize; 2, with GNyq) 
        flag_degradation=pval;
    elseif strcmpi(pname,'flag_equalization')         % Choice on which histogram to equalize the HR image, 0=classic; 1=with GT; 2=regression based; 3=PCA regression based
        flag_equalization=pval;
    elseif any(strcmpi(pname,{'ratio_tgt_REF','ratio_REF','ratio_GT'})) % Target ratio for the sharpened image (this program performs the fusion first, then interpolates the result)
        ratio_tgt_REF=pval;
    elseif any(strcmpi(pname,{'flag_interpolation','interpolation','upsample','interp'}))       % Flag defining the method of interpolation
        flag_interpolation=pval;
    elseif strcmpi(pname,'flag_imresize_RR')          % Flag to decide how to upscale the LR image in case of Reduced Resolution validation
        flag_imresize_RR=pval;
    elseif any(strcmpi(pname,{'methods_list','methods'}))              % List of fusion methods to test (as a cell)
        methods_list=load_methodslist(pval);
    elseif strcmpi(pname,'printMAT')                  % Name of the output file (will be saved in output folder)
        printMAT=pval; 
    elseif strcmpi(pname,'output_file')               % Name of the output file (will be saved in output folder)
        output_file=pval; 
    end
end

%% Upsampling (or downsampling) of Low Resolution image
if isempty(I_MS)
    [I_MS,~,time_EXP]=scale_image(I_MS_LR,ratio,flag_interpolation,[],'edge',edge_cut);
end
if isempty(I_MS_LR)
    I_MS_LR=I_MS(start_pos(1):ratio:end,start_pos(2):ratio:end,:);
else
    I_MS_LR=I_MS_LR(edge_cut+1:end-edge_cut,edge_cut+1:end-edge_cut,:);
end


if ratio==1 || isempty(I_PAN)
    % Condition for no fusion
    for ii=1:numel(methods_list), methods_list{ii}='none'; end
    I_PAN_LR=[]; I_PAN_eq=[];
else
    % Low-resolution PAN for equalization
    if flag_FR==0 && ~isempty(I_GT)
        I_GT_LR=scale_image(I_GT,size(I_PAN,1)/size(I_GT,1),flag_interpolation,flag_imresize_RR,'GNyq',GNyq);
        [ I_PAN_LR, I_PAN_eq ] = PAN_equalization( I_PAN, I_MS, flag_degradation, flag_equalization, ratio, GNyq, GNyq_PAN, I_GT_LR );
    else
        [ I_PAN_LR, I_PAN_eq ] = PAN_equalization( I_PAN, I_MS, flag_degradation, flag_equalization, ratio, GNyq, GNyq_PAN);
    end
end

%% Data Fusion
NumOutputs=numel(methods_list);
L1=size(I_MS,1); L2=size(I_MS,2); Nb=size(I_MS,3);
MatrixImage=zeros(L1,L2,Nb,NumOutputs);
time=zeros(NumOutputs,1);
methods_list_idx=1:NumOutputs;
titleImages=strrep(methods_list,'_','-');
cd(fusion_folder);
for ii=1:NumOutputs
    method=methods_list{ii};
    switch method
        case 'none'
            MatrixImage(:,:,:,ii) = I_MS;
            time(ii) = 0;
        case 'EXP'
            MatrixImage(:,:,:,ii) = I_MS;
            time(ii) = time_EXP;
        case 'PCA'
            cd CS
            t2=tic;
            MatrixImage(:,:,:,ii) = PCA_mod(I_MS,I_PAN,I_PAN_LR);
            time(ii) = toc(t2);
            cd ..
        case 'IHS'
            cd CS
            t2=tic;
            MatrixImage(:,:,:,ii) = IHS_mod(I_MS,I_PAN,I_PAN_LR);
            time(ii) = toc(t2);
            cd ..
        case 'Brovey'
            cd CS
            t2=tic;
            MatrixImage(:,:,:,ii) = Brovey_mod(I_MS,I_PAN,I_PAN_LR);
            time(ii) = toc(t2);
            cd ..
        case 'BDSD'
            gcd_temp=gcd(L1,L2);
            if gcd_temp>=128
                BDSD_blocksize=127+find(rem(gcd_temp,128:gcd_temp)==0,1);
            elseif gcd_temp>=32
                BDSD_blocksize=gcd_temp;
            else
                BDSD_blocksize=0;
            end
            if BDSD_blocksize>=Nb+1
                cd BDSD
                t2=tic;
                MatrixImage(:,:,:,ii) = BDSD_mod(I_MS,I_PAN,ratio,BDSD_blocksize,GNyq_PAN,GNyq,start_pos);
                time(ii) = toc(t2);
                cd ..
            else
                methods_list_idx=methods_list_idx(methods_list_idx~=ii);
                fprintf('BDSD was suppressed (try using G-BDSD instead)\n');
            end
        case 'GS'
            cd CS
            t2=tic;
            MatrixImage(:,:,:,ii) = GS_mod(I_MS,I_PAN,I_PAN_LR);
            time(ii) = toc(t2);
            cd ..
        case 'GSA'
            cd GS
            t2=tic;
            MatrixImage(:,:,:,ii) = GSA(I_MS,I_PAN,I_MS_LR,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'GSA_local'
            cd GS
            t2=tic;
            MatrixImage(:,:,:,ii) = GSA_local(I_MS,I_PAN,I_MS_LR,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'PRACS'
            cd PRACS
            t2=tic;
            MatrixImage(:,:,:,ii) = PRACS(I_MS,I_PAN,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'HPF'
            cd HPF
            t2=tic;
            MatrixImage(:,:,:,ii) = HPF(I_MS,I_PAN,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'SFIM'
            cd CS
            t2=tic;
            MatrixImage(:,:,:,ii) = SFIM_mod(I_MS,I_PAN,ratio,I_PAN_eq);
            time(ii) = toc(t2);
            cd ..
        case 'Indusion'
            cd Indusion
            t2=tic;
            MatrixImage(:,:,:,ii) = Indusion(I_PAN,I_MS_LR,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'ATWT'
            cd Wavelet
            t2=tic;
            MatrixImage(:,:,:,ii) = ATWT(I_MS,I_PAN,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'AWLP'
            cd Wavelet
            t2=tic;
            MatrixImage(:,:,:,ii) = AWLP(I_MS,I_PAN,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'ATWT_M2'
            cd Wavelet
            t2=tic;
            MatrixImage(:,:,:,ii) = ATWT_M2(I_MS,I_PAN,ratio);
            time(ii) = toc(t2);
            cd ..
        case 'ATWT_M3'
            cd Wavelet
            t2=tic;
            % MatrixImage(:,:,:,ii) = ATWT_M3(I_MS,I_PAN,ratio); 
            %%% There was no warning for complex numbers in RR, but I changed to real regardless
            MatrixImage(:,:,:,ii) = real(ATWT_M3(I_MS,I_PAN,ratio));
            time(ii) = toc(t2);
            cd ..
        case 'ATWT_mod'
            cd Wavelet
            t2=tic;
            MatrixImage(:,:,:,ii) = ATWT_mod(I_MS,I_PAN,ratio,I_PAN_eq);
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP'
            cd GLP
            t2=tic;
            MatrixImage(:,:,:,ii) = MTF_GLP(I_PAN,I_MS,GNyq,ratio);
            % MatrixImage(:,:,:,ii) = MTF_GLP_ECB(I_MS,I_PAN,ratio,[9 9],2.5,GNyq);
            % MatrixImage(:,:,:,ii) = MTF_GLP_CBD(I_MS,I_PAN,ratio,[55 55],-Inf,GNyq);
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP_HPM_PP'
            % if 2^nextpow2(ratio)==ratio
            if mod(ratio,2)==0 && ratio>2 
                cd GLP
                t2=tic;
                MatrixImage(:,:,:,ii) = MTF_GLP_HPM_PP_(I_PAN,I_MS_LR,GNyq,ratio);
                time(ii) = toc(t2);
                cd ..
            else
                methods_list_idx=methods_list_idx(methods_list_idx~=ii);
                fprintf('MTF-GLP-HPM With PostProcessing was suppressed\n');
            end
        case 'MF_HG'
            individuals=[1 0 1 0 0 0]; %unused
            % ReDefines also operator after MF_Adapt_Settings 1st and 2nd bits:
            operators = [0 0 0];
            % 000: gradient; 010: top-hat; 100: toggle; 110: Bai (toggle+top-hat)
            % 001: LOCO; 011:median; 101: reconstuction
            textse=repmat([ 0 1 0 1 1 1 0 1 0 0 1],1,3);
            cd MF_Adapt
            save chromosomes_file textse individuals operators
            t2=tic;
            if flag_FR==0
                % [MatrixImage(:,:,:,ii),dum,D_MF_HG,P_LP_MF_HG]=MF_Adapt_gradLoad(I_PAN,I_MS,L,ratio,sensor,I_GT_LR);%whole
                [MatrixImage(:,:,:,ii),~,~,~]=MF_Adapt_gradLoad(I_PAN,I_MS,L,ratio,sensor,I_GT_LR);%whole
            else
                %[MatrixImage(:,:,:,ii),dum,D_MF_HG,P_LP_MF_HG]=MF_Adapt_gradLoad(I_PAN,I_MS,L,ratio,sensor,[]);%whole
                [MatrixImage(:,:,:,ii),~,~,~]=MF_Adapt_gradLoad(I_PAN,I_MS,L,ratio,sensor,[]);%whole
            end
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP_HPM'
            cd GLP
            t2=tic;
            % MatrixImage(:,:,:,ii) = MTF_GLP_HPM(I_PAN,I_MS,GNyq,ratio);
            % MatrixImage(:,:,:,ii) = MTF_GLP_HPM2(I_PAN,I_MS,GNyq,ratio);  %Tryouts
            MatrixImage(:,:,:,ii) = MTF_GLP_HPM_mod(I_PAN,I_MS,GNyq,ratio,I_PAN_eq);
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP_CBD'
            cd GS
            t2=tic;
            MatrixImage(:,:,:,ii) = GS2_GLP_mod(I_MS,I_PAN,ratio,GNyq);
            time(ii) = toc(t2);
            cd ..
        case 'GFPCA'
            No_PCs=4;    % number of principal components used by guided filtering with pan/RGB image, e.g. No_PCs=5
            r=8;         % the size of local sliding window, e.g. r=8
            eps=0.001^2; % regularization parameter determining the degree of blurring for the guided filter, e.g. eps = 0.01^2        
            cd GFPCA
            t2=tic;
            MatrixImage(:,:,:,ii) = GFPCA(I_MS_LR, I_PAN, min(No_PCs,size(I_MS_LR,3)), r, eps);
            time(ii) = toc(t2);
            cd ..
        case 'CNMF'
            cd CNMF
            t2=tic;
            MatrixImage(:,:,:,ii) = CNMF_fusion(I_MS_LR,I_PAN);
            time(ii) = toc(t2);
            cd ..
        case {'BayesNaive','BayesSparse'}
            if strcmpi(method,'BayesNaive'), prior='Gaussian'; end
            if strcmpi(method,'BayesSparse'), prior='Sparse'; end
            cd bayesfusion
            t2=tic;
            setup;
            MatrixImage(:,:,:,ii)= BayesianFusion_mod2(I_MS_LR,I_PAN,{Band_overlap_PAN},KerBlu0,ratio,prior,start_pos);
            time(ii) = toc(t2);
            cd ..
        case {'BayesNaive_FE','BayesSparse_FE','BayesML','BayesNLM'}
            if strcmpi(method,'BayesNaive_FE'), prior='Gaussian'; end
            if strcmpi(method,'BayesSparse_FE'), prior='Sparse'; end
            if strcmpi(method,'BayesML'), prior='ML'; end
            if strcmpi(method,'BayesNLM'), prior='NLM'; end
            cd BlindFuse-master
            t2=tic;
            MatrixImage(:,:,:,ii)=main_BlindFuse_mod(I_PAN,I_MS_LR,ratio,{Band_overlap_PAN},start_pos,prior,'PCA');
            time(ii) = toc(t2);
            cd ..
        case 'HySure'
            cd HySure
            t2=tic;
            % MatrixImage(:,:,:,ii) = HySure_wrapper_mod(I_MS_LR, I_PAN, ratio, Band_overlap_PAN);
            MatrixImage(:,:,:,ii) = HySure_wrapper_MS(I_MS_LR,I_PAN,ratio,{Band_overlap_PAN},start_pos,im_tag);
            time(ii) = toc(t2);
            cd ..
        case 'G_BDSD'
            cd BDSD
            t2=tic;
            MatrixImage(:,:,:,ii) = GBDSD_Fusion(I_PAN,I_MS,ratio,GNyq_PAN,GNyq);
            time(ii) = toc(t2);
            cd ..
        case 'L_BDSD'
            cd BDSD
            t2=tic;
            MatrixImage(:,:,:,ii) = IBDSD_Fusion(I_PAN,I_MS,ratio,GNyq_PAN,GNyq);
            % MatrixImage(:,:,:,ii) = LBDSD_Fusion(I_PAN,I_MS,ratio,GNyq_PAN,GNyq);
            time(ii) = toc(t2);
            cd ..
        case 'C_BDSD'
            cd BDSD
            t2=tic;
            MatrixImage(:,:,:,ii) = CBDSD_Fusion(I_PAN,I_MS,ratio,GNyq_PAN,GNyq);
            % MatrixImage(:,:,:,ii) = SCBDSD_Fusion(I_PAN,I_MS,ratio,GNyq_PAN,GNyq);
            time(ii) = toc(t2);
            cd ..
    end
    if ismember(ii,methods_list_idx)
        fprintf('Elaboration time %s: %.2f [sec]\n',titleImages{ii},time(ii));
    end
end
cd(original_folder);
        
%% Removal of failed fusions
methods_list=methods_list(methods_list_idx);
time=time(methods_list_idx);
MatrixImage=MatrixImage(:,:,:,methods_list_idx);
titleImages(methods_list_idx);

%% Upsampling of fused images
MatrixImage_in=MatrixImage;
MatrixImage=[];
for ii=1:size(MatrixImage_in,4)
    MatrixImage_temp=scale_image(MatrixImage_in(:,:,:,1),ratio_tgt_REF/ratio,...
    flag_interpolation,flag_imresize_RR,'GNyq',GNyq);
    MatrixImage_in(:,:,:,1)=[];
    MatrixImage=cat(4,MatrixImage,MatrixImage_temp);
end

%% Cuts values outside of radiometric range
if flag_thvalues==1
    MatrixImage(MatrixImage > 2^L) = 2^L;
    MatrixImage(MatrixImage < 0) = 0;
end

%% Generation of the Output

Fused.data=MatrixImage; clear MatrixImage;
Fused.methods_fusion=methods_list;
Fused.titleImages=titleImages;
Fused.fusion_type='P';
Fused.flag_degradation=flag_degradation;
Fused.flag_equalization=flag_equalization;
Fused.titleImages=titleImages;
Fused.ratio=ratio_tgt_REF;
Fused.flag_interpolation=flag_interpolation;
Fused.ratio_upscale=ratio_tgt_REF/ratio;
Fused.time=time;
Fused.flag_thvalues=flag_thvalues;
Fused.IntBits=L;
Fused.sensor=sensor;
Fused.edge_cut=0;
Fused.label='Fused';
if isfield(Image_MS,'type'), Fused.type=Image_MS.type; else, Fused.type='MS'; end
if isfield(Image_MS,'Bands_to_display'), Fused.Bands_to_display=Image_MS.Bands_to_display; else, Fused.Bands_to_display=1:3; end
if isfield(Image_MS,'Bands_selected'), Fused.Bands_selected=Image_MS.Bands_selected; else, Fused.Bands_selected=1:size(MatrixImage.data,3); end
if isfield(Image_MS,'GSD'), Fused.MS=Image_PAN.GSD*ratio_tgt_REF; end


if printMAT==1
   save( fullfile(output_folder,output_file), 'Fused','-append'); 
end
