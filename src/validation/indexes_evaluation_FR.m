%% REDUCED RESOLUTION QUALITY INDICES EVALUATION
% Description: 
%     Reduced resolution quality indexes calculation. 
% 
% Interface:
%     [MR,qi_label] = indexes_evaluation_FR(F,PAN,I_Field,Value)
%     eg: [MR,qi_label]= indexes_evaluation_RR(F,PAN,'qindex_list',{'ssim','psnr'})
%
% Inputs:
%     F:                  Fused Image; (Its first three sizes must be the same as the ground truth; the function manages multiple tests on the 4th dimensions)
%     GT:                 A struct whose fields are:
%         image:          The Ground Truth image itself
%         ratio:          Scale ratio between MS and PAN.
%         IntBits:        Image radiometric resolution; 
%         Qblocks_size:   Block size of the Q-index locally applied;
%
% Description of Field, Value options
%     'qindex_list':      List of quality indices to perform the tests (default: {'Q2^n','SAM','ERGAS','SCC','UIQI'});
%                         (some default options for validation can be selected with numerical vaules 0-4)
%     'flag_cut_bounds':  Cut the boundaries of the viewed Panchromatic image;
%     'dim_cut':          Define the dimension of the boundary cut;
%     'th_values':        Flag. If th_values == 1, apply an hard threshold to the dynamic range.
%
% Outputs:
%     MR:        Matrix of the computed quality indices
%                (sizes: Ni x Nt, where Ni is the length of qindex_list, Nt the 4th dimension of F)
%     qi_label:  Labels for Quality Indices
% 
% Available choices for index computation in qindex_list:
%     'SAM':     Spectral Angle Mapper (SAM) index;
%     'ERGAS':   Erreur Relative Globale Adimensionnelle de Synthèse (ERGAS) index;
%     'sCC':     Spatial Correlation Coefficient between fused and ground-truth images;
%     'Q2n':     Q2n index;
%     'SSIM':    Structural Similarity (mean over bands);
%     'UIQI':    Universal Image Quality Index (mean over bands);
%     'PSNR':    Peak Signal to Noise Ratio;
%     'MSE':     Mean Square Error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MatrixResults,qindex_list] = indexes_evaluation_FR(F,PAN,EXP,varargin)

% current_folder=pwd;
% q2nnew_folder=[fileparts(mfilename('fullpath')),'\Q2n_new\'];

%% Parameters implementation

flag_L=0;
flag_ratio=0;

if isfield(F,'image'), F=F.image; end
if isfield(EXP,'image'), I_EXP=EXP.image; else, I_EXP=EXP; end
if isfield(PAN,'ratio'), ratio=PAN.ratio; flag_ratio=1; else, ratio=4; end
if isfield(EXP,'IntBits'), L=GT.IntBits; flag_L=1; else, L=11; end
if isfield(EXP,'DynamicRange'), DynamicRange=EXP.DynamicRange; flag_L=1; else, DynamicRange=2^L-1; end
if isfield(EXP,'Qblocks_size'), Qblocks_size=GT.Qblocks_size; else, Qblocks_size=32; end

qindex_list={'D_lambda','D_S','QNR','SAM','sCC'};
flag_cutbounds=0;
dim_cut=0;
th_values=0;

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};   
    if any(strcmpi(pname,{'ratio','ratio_tgt_REF'}))     % Label for a mask to generate mosaicked observation
       ratio=pval;
       flag_ratio=1;
    elseif any(strcmpi(pname,{'L','bpp'}))     % Shift to mask
       L=pval;
       flag_L=1;
    elseif any(strcmpi(pname,{'flag_cutbounds','flag_cutbounds_val','flag_dimcut'}))
        flag_cutbounds=pval;
    elseif strcmpi(pname,'dim_cut')
        dim_cut=pval;
    elseif any(strcmpi(pname,{'flag_thvalues','flag_thvalues_val','th_values'}))
        th_values=pval;    
    elseif any(strcmpi(pname,{'Q_blocks_size','Q_blocks'}))
        Qblocks_size=pval;
    elseif any(strcmpi(pname,{'qindex_list','qindex','index','index_list'}))
        qindex_list=pval;
        if qindex_list==4
            qindex_list={'SCCo_laplacian','SCCo_sobel','SCCo_sobel2',...
                'SCCo_prewitt','SCCo_prewitt2','SCC_laplacian',...
                'SCC_sobel','SCC_sobel2','SCC_prewitt','SCC_prewitt2'};
        elseif qindex_list==3
            qindex_list={'Q2n_b16s16','Q2^n_b32s32','Q2n_b32s8','Q2n_b36s36',...
                'Q2n_b32s32c1','Q2n_b32s32c10','Q2n_alt_b36s36','Q2n_alt_b36s36'};
        elseif qindex_list==2
            qindex_list={'Q2^n','SAM','ERGAS','SCC','UIQI','Q2n_alt'};
        elseif qindex_list==1
            qindex_list={'Q2^n','SAM','ERGAS','SCC','UIQI'};
        end
    end
end

%% Cutting bounds

if flag_cutbounds
    I_EXP = EXP(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    F = F(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
end

%% Thresholding of radiometric intensity
    
if th_values
    F(F > DynamicRange) = DynamicRange;
    F(F < 0) = 0;
    if flag_L==0, warning('Bit per pixel not setup for thresolding'); end
end


%% Definition of consistent labels for Quality indices

N_index=numel(qindex_list);

blocknum=zeros(1,N_index);
shiftnum=zeros(1,N_index);
cut=zeros(1,N_index);
edge=cell(1,N_index);

for jj=1:N_index
    if any(strcmpi(qindex_list{jj},{'SCC','CC'}))
        qindex_list{jj}='sCC';
    elseif any(strcmpi(qindex_list{jj},{'SAM'}))
        qindex_list{jj}='SAM';
    elseif strncmpi(qindex_list{jj},'SCC_',4)
        idx_unds=find(qindex_list{jj},'_');
        qindex_list_temp=qindex_list{jj};
        edge{jj}=qindex_list_temp(idx_unds+1:end);
        if any(strcmpi(edge{jj},{'laplacian','Lap','{Lap}','{Laplacian}'}))
            edge{jj}='laplacian'; edge2='{Lap}';
        elseif any(strcmpi(edge{jj},{'prewitt','pre','{pre}','{prewitt}'}))
            edge{jj}='prewitt'; edge2='{Pre}';
        elseif any(strcmpi(edge{jj},{'prewitt2','pre2','{pre2}','{pr2}','{prewitt2}'}))
            edge{jj}='prewitt2'; edge2='{Pr2}';
        elseif any(strcmpi(edge{jj},{'sobel','sob','{sobel}','{sob}'}))
            edge{jj}='sobel'; edge2='{Sob}';
        elseif any(strcmpi(edge{jj},{'sobel2','sob2','{sobel2}','{sob2}','{so2}'}))
            edge{jj}='sobel2'; edge2='{So2}';
        end
        qindex_list{jj}=['sCC_',edge2];
    elseif strncmpi(qindex_list{jj},'SCCold_',7) || strncmpi(qindex_list{jj},'SCCo_',5)
        idx_unds=find(qindex_list{jj},'_');
        qindex_list_temp=qindex_list{jj};
        edge{jj}=qindex_list_temp(idx_unds+1:end);
        if any(strcmpi(edge{jj},{'laplacian','Lap','{Lap}','{Laplacian}'}))
            edge{jj}='laplacian'; edge2='{Lap}';
        elseif any(strcmpi(edge{jj},{'prewitt','pre','{pre}','{prewitt}'}))
            edge{jj}='prewitt'; edge2='{Pre}';
        elseif any(strcmpi(edge{jj},{'prewitt2','pre2','{pre2}','{pr2}','{prewitt2}'}))
            edge{jj}='prewitt2'; edge2='{Pr2}';
        elseif any(strcmpi(edge{jj},{'sobel','sob','{sobel}','{sob}'}))
            edge{jj}='sobel'; edge2='{Sob}';
        elseif any(strcmpi(edge{jj},{'sobel2','sob2','{sobel2}','{sob2}','{so2}'}))
            edge{jj}='sobel2'; edge2='{So2}';
        end
        qindex_list{jj}=['sCCo_',edge2];
    end
end

%% Calculation of Quality Index

%% Fix for I_F being more than 4 dimensions
sizes_F=size(F);
if length(sizes_F)<=2, sizes_F(3)=1; end
if length(sizes_F)<=3, sizes_F(4)=1; end
N_tests=prod(sizes_F(4:end));
F=reshape(F,[sizes_F(1:3),N_tests]);

MatrixResults=zeros(N_index,N_tests);
I_PANr=repmat(I_PAN,[1,1,size(I_EXP,3)]);

fprintf('Calculating quality indices:   0/%3d',N_tests);

for kk=1:N_tests
    I_F=F(:,:,:,kk);
    for jj=1:N_index
        if any(strcmpi(qindex_list{jj},{'SCC','CC'}))
            MatrixResults(jj,kk)=SCC(I_PANr,I_F,'laplacian','global1');
        elseif any(strcmpi(qindex_list{jj},{'SAM'}))
            MatrixResults(jj,kk)=SAM(I_EXP,I_F);
        elseif strncmpi(qindex_list{jj},'SCC_',4)
            MatrixResults(jj,kk)=SCC(I_PANr,I_F,edge{jj},'global1');
        elseif strncmpi(qindex_list{jj},'SCCold_',7) || strncmpi(qindex_list{jj},'SCCo_',5)
            MatrixResults(jj,kk)=SCC(I_PANr,I_F,edge{jj});
        elseif any(strcmpi,{'D_s','Ds'})
            
        end
    end
    fprintf([repmat('\b',[1,7]),'%3d/%3d'],kk,N_tests);
end
fprintf('. Done!\n');

MatrixResults=reshape(MatrixResults,[N_index,sizes_F(4:end)]);

end