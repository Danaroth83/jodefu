% HYPERSPECTRAL PANSHARPENING
%
% This is a demo code to run several hyperspectral pansharpening methods
% presented in the following paper.
%
% Laetitia Loncan, Luis B. Almeida, Jose M. Bioucas-Dias, Xavier Briottet, 
% Jocelyn Chanussot, Nicolas Dobigeon, Sophie Fabre, Wenzhi Liao, 
% Giorgio A. Licciardi, Miguel Simoes, Jean-Yves Tourneret, 
% Miguel A. Veganzones, Gemine Vivone, Qi Wei and Naoto Yokoya, 
% "Introducing hyperspectral pansharpening," Geoscience and Remote Sensing
% Magazine, 2015.

% clc
clearvars
close all
addpath(genpath('../Methods'));
addpath(genpath('../Quality_Indices'));
addpath(genpath('../Outputs'));
%% Load dataset

p = fileparts(mfilename('fullpath'));

PRECISION    = 'double';
OFFSET       = 0 ;
INTERLEAVE   = 'bsq';
BYTEORDER    = 'ieee-le';

FILENAME_REF = [p '\REF']; % where you put the data
% SIZE_REF     = [396,184,176]; 
SIZE_REF     = [500,500,113];
I_REF        = multibandread(FILENAME_REF, SIZE_REF, PRECISION, OFFSET, INTERLEAVE, BYTEORDER);

% start_r=5;
% start_c=5;
% I_REF=padarray(I_REF,[start_r-1 start_c-1],'circular','post');
% I_REF=I_REF(start_r:end,start_c:end,:);

% FILENAME_HS = [p '\HS']; % where you put the data
% % SIZE_HS     = [99,46,176]; 
% SIZE_HS     = [125,125,113];
% I_HS        = multibandread(FILENAME_HS, SIZE_HS, PRECISION, OFFSET, INTERLEAVE, BYTEORDER);

% FILENAME_PAN = [p '\PAN']; % where you put the data
% % SIZE_PAN     = [396,184,1];
% SIZE_PAN     = [500,500,1];
% I_PAN        = multibandread(FILENAME_PAN, SIZE_PAN, PRECISION, OFFSET, INTERLEAVE, BYTEORDER);

% cd ..
% cd Methods
%% Generating the HS and PAN image from the reference image
% ratio = size(I_PAN,1)/size(I_HS,1);
ratio = 5;
% overlap = 1:41; % commun bands (or spectral domain) between I_PAN and I_HS
overlap = 1:25; 
size_kernel=[9 9];
sig = (1/(2*(2.7725887)/ratio^2))^0.5;
start_pos(1)=1; % The starting point of downsampling
start_pos(2)=1; % The starting point of downsampling
% start_pos(1)=(ratio+1)/2;
% start_pos(2)=(ratio+1)/2;
[I_HS,KerBlu]=conv_downsample(I_REF,ratio,size_kernel,sig,start_pos);
I_PAN = mean(I_REF(:,:,overlap),3);
%% Hyperspectral pansharpening

%%%%%% NOTE: Put the value of L for the different datasets.
L = 15; %Sensor Radiometric Resolution, obtained by finding the minimal value n which
% verified max(max(max(I_HS)) < 2^n
% dataset1 => L = 13
% dataset2 => L = 15


% SFIM
cd 'SFIM'
tic
I_SFIM = SFIM(I_HS,I_PAN,ratio,1,L);
disp('Comp. time (SFIM): ');
toc
cd ..
QI_SFIM = QualityIndices(I_SFIM,I_REF,ratio);
cd ..
multibandwrite(I_SFIM,'Outputs/I_SFIM',INTERLEAVE);
cd Methods


% MTF GLP
cd 'GLP'
tic
I_MTF_GLP = MTF_GLP(I_HS,I_PAN,'GaussKernel',ratio,1,L);
disp('Comp. time (MTF GLP): ');
toc
cd ..
QI_MTF_GLP = QualityIndices(I_MTF_GLP,I_REF,ratio);
cd ..
multibandwrite(I_MTF_GLP,'Outputs/I_MTF_GLP',INTERLEAVE);
cd Methods

% MTF GLP HPM
cd 'GLP'
tic
I_MTF_GLP_HPM = MTF_GLP_HPM(I_HS,I_PAN,'GaussKernel',ratio,1,L);
disp('Comp. time (MTF GLP HPM): ');
toc
cd ..
QI_MTF_GLP_HPM = QualityIndices(I_MTF_GLP_HPM,I_REF,ratio);
cd ..
multibandwrite(I_MTF_GLP_HPM,'Outputs/I_MTF_GLP_HPM',INTERLEAVE);
cd Methods


% GS
cd 'GS'
tic
I_GS = GS(I_HS,I_PAN);
disp('Comp. time (GS): ');
toc
cd ..
QI_GS = QualityIndices(I_GS,I_REF,ratio);
cd ..
multibandwrite(I_GS,'Outputs/I_GS',INTERLEAVE);
cd Methods


% GSA
cd 'GS'
tic
I_GSA = GSA(I_HS,I_PAN);
disp('Comp. time (GSA): ');
toc
cd ..
QI_GSA = QualityIndices(I_GSA,I_REF,ratio);
cd ..
multibandwrite(I_GSA,'Outputs/I_GSA',INTERLEAVE);
cd Methods


% PCA
cd 'PCA'
tic
I_PCA = PCA(I_HS,I_PAN);
disp('Comp. time (PCA): ');
toc
cd ..
QI_PCA = QualityIndices(I_PCA,I_REF,ratio);
cd ..
multibandwrite(I_PCA,'Outputs/I_PCA',INTERLEAVE);
cd Methods


% GFPCA

% CNMF
cd 'CNMF'
tic
I_CNMF = CNMF_fusion(I_HS,I_PAN);
disp('Comp. time (CNMF): ');
toc
cd ..
QI_CNMF = QualityIndices(I_CNMF,I_REF,ratio); % measure cc, sam, rmse, ergas
cd ..
multibandwrite(I_CNMF,'Outputs/I_CNMF',INTERLEAVE);
cd Methods


% Bayesian Naive
cd 'BayesFusion'
setup;
tic
[I_BayesNaive]= BayesianFusion(I_HS,I_PAN,overlap,KerBlu,ratio,'Gaussian',start_pos);
toc
cd ..
QI_BayesNaive = QualityIndices(I_BayesNaive,I_REF,ratio);
cd ..
multibandwrite(I_BayesNaive,'Outputs/I_BayesNaive',INTERLEAVE);
cd Methods

% Bayesian Sparse
cd 'BayesFusion'
setup;
tic;
[I_BayesSparse]= BayesianFusion(I_HS,I_PAN,overlap,KerBlu,ratio,'Sparse',start_pos);
toc
cd ..
QI_BayesSparse = QualityIndices(I_BayesSparse,I_REF,ratio);
cd ..
multibandwrite(I_BayesSparse,'Outputs/I_BayesSparse',INTERLEAVE);
cd Methods

% HySure
cd 'HySure'
tic
I_HySure = HySure_wrapper(I_HS, I_PAN, ratio, overlap);
disp('Comp. time (HySure): ');
toc
cd ..
QI_HySure = QualityIndices(I_HySure, I_REF, ratio); % measure cc, sam, rmse, ergas
cd ..
multibandwrite(I_HySure,'Outputs/I_HySure',INTERLEAVE);
cd Methods


% [rmse_total, ergas, sam, uiqi] = quality_assessment(I_REF, I_HySure, 4, 1/ratio);
% fprintf('Quality indices:\n RMSE = %2.3f\n ERGAS = %2.3f\n SAM = %2.3f\n UIQI = %2.3f\n', ...
%     rmse_total,ergas,sam,uiqi)
