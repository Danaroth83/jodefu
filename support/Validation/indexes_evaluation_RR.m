%% REDUCED RESOLUTION QUALITY INDICES EVALUATION
% Description: 
%     Reduced resolution quality indexes calculation. 
% 
% Interface:
%     [MR,qi_label,order] = indexes_evaluation_RR(F,GT,Field,Value)
%     eg: [MR,qi_label,order]= indexes_evaluation_RR(F,GT,'qindex_list',{'ssim','psnr'})
%
% Inputs:
%     F:    Fused Image; (Its first three sizes must be the same as the ground truth; the function manages multiple tests on the 4th dimensions)
%     GT:   A struct whose fields are:
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
%     order:     A matrix of sizes Ni x Nt, for which each row defines the
%                index, in increasing order, of the best results
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

function [MatrixResults,qindex_list,MatrixResults_order] = indexes_evaluation_RR(F,GT,varargin)

current_folder=pwd;
q2nnew_folder=fullfile(fileparts(mfilename('fullpath')),'Q2n_new');

%% Parameters implementation

flag_L=0;
flag_ratio=0;

if isfield(F,'data'), F=F.data; end
if isfield(GT,'data'), I_GT=GT.data; else, I_GT=GT; end
if isfield(GT,'ratio'), ratio=GT.ratio; flag_ratio=1; else, ratio=4; end
if isfield(GT,'IntBits'), L=GT.IntBits; flag_L=1; else, L=11; end
if isfield(GT,'DynamicRange'), DynamicRange=GT.DynamicRange; flag_L=1; else, DynamicRange=2^L-1; end
if isfield(GT,'Qblocks_size'), Qblocks_size=GT.Qblocks_size; else, Qblocks_size=32; end

qindex_list={'Q2^n','SAM','ERGAS','SCC','UIQI'};
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
        if iscell(qindex_list)
            qindex_list=pval;
        elseif ischar(qindex_list)
            qindex_list={pval};
        elseif qindex_list==4
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
    I_GT = I_GT(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
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
    if any(strcmpi(qindex_list{jj},{'Q2n','Q2^n','Q'}))
        qindex_list{jj}='Q2^n';
    elseif any(strcmpi(qindex_list{jj},{'SCC','CC'}))
        qindex_list{jj}='sCC';
    elseif any(strcmpi(qindex_list{jj},{'PSNR'}))
        qindex_list{jj}='PSNR';
    elseif any(strcmpi(qindex_list{jj},{'PSNRmax','PSNRm'}))
        qindex_list{jj}='PSNRm';
    elseif any(strcmpi(qindex_list{jj},{'MSE'}))
        qindex_list{jj}='MSE';
    elseif any(strcmpi(qindex_list{jj},{'ERGAS'}))
        qindex_list{jj}='ERGAS';
        if flag_ratio==0, warning('Scale ratio for ERGAS non inputted, automatically set as 4'); end
     elseif any(strcmpi(qindex_list{jj},{'Q_avg','UIQI','Q_{avg}'}))
        qindex_list{jj}='Q_{avg}';
    elseif any(strcmpi(qindex_list{jj},{'SAM'}))
        qindex_list{jj}='SAM';
    elseif any(strcmpi(qindex_list{jj},{'Q2n_alt','Q2^n_alt','Q2n_{alt}','Q2^n_{alt}'}))
        qindex_list{jj}='Q2^n_{alt}';
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
    elseif strncmpi(qindex_list{jj},'Q2n_b',5) || strncmpi(qindex_list{jj},'Q2^n_b',6)
        qindex_list{jj}=strrep(qindex_list{jj},'^','');
        qindex_list_temp=qindex_list{jj};
        idx_c=find(qindex_list_temp=='c');
        if isempty(idx_c), cut(jj)=0; else, cut(jj)=qindex_list_temp(idx_c+1:end); idx_c=length(qindex_list_temp)+1; end
        idx_s=find(qindex_list_temp=='s');
        blocknum(jj)=str2double(qindex_list_temp(6:idx_s-1));
        shiftnum(jj)=str2double(qindex_list_temp(idx_s+1:idx_c-1));
        if cut(jj)==0
            qindex_list{jj}=sprintf('Q2^n_{b%d}^{s%d}',blocknum(jj),shiftnum(jj));
        else
            qindex_list{jj}=sprintf('Q2^n_{b%d}^{s%dc%d}',blocknum(jj),shiftnum(jj),cut(jj));
        end
    elseif strncmpi(qindex_list{jj},'Q2n_alt_b',9) || strncmpi(qindex_list{jj},'Q2^n_alt_b',10) || ...
            strncmpi(qindex_list{jj},'Q2n_{alt}_b',11) || strncmpi(qindex_list{jj},'Q2^n_{alt}_b',12)
        qindex_list{jj}=strrep(strrep(strrep(qindex_list{jj},'^',''),'{',''),'}','');
        qindex_list_temp=qindex_list{jj};
        idx_s=find(qindex_list_temp=='s');
        idx_c=find(qindex_list_temp=='c');
        if isempty(idx_c), cut(jj)=0; else, cut(jj)=qindex_list_temp(idx_c+1:end); idx_c=length(qindex_list_temp)+1; end
        blocknum(jj)=str2double(qindex_list_temp(9:idx_s-1));
        shiftnum(jj)=str2double(qindex_list_temp(idx_s+1:idx_c-1));
        if cut(jj)==0
            qindex_list{jj}=sprintf('Q2^n_{altb%d}^{s%d}',blocknum(jj),shiftnum(jj));
        else
            qindex_list{jj}=sprintf('Q2^n_{altb%d}^{s%dc%d}',blocknum(jj),shiftnum(jj),cut(jj));
        end
    elseif strncmpi(qindex_list{jj},'ssim',4)
        qindex_list{jj}='SSIM';
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
MatrixResults_order=zeros(N_index,N_tests);

fprintf('Calculating quality indices:   0/%3d',N_tests);

for kk=1:N_tests
    I_F=F(:,:,:,kk);
    for jj=1:N_index
        if any(strcmpi(qindex_list{jj},{'Q2n','Q2^n','Q'}))
            MatrixResults(jj,kk)=q2n(I_GT,I_F,Qblocks_size,Qblocks_size);
        elseif any(strcmpi(qindex_list{jj},{'SCC','CC'}))
            MatrixResults(jj,kk)=SCC(I_GT,I_F,'laplacian','global1');
        elseif any(strcmpi(qindex_list{jj},{'PSNR'}))
            MatrixResults(jj,kk)=PSNR(I_GT,I_F,L);
        elseif any(strcmpi(qindex_list{jj},{'PSNRmax','PSNRm'}))
            MatrixResults(jj,kk)=PSNR(I_GT,I_F,[]);
        elseif any(strcmpi(qindex_list{jj},{'MSE'}))
            MatrixResults(jj,kk)=MSE(I_GT,I_F);
        elseif any(strcmpi(qindex_list{jj},{'ERGAS'}))
            MatrixResults(jj,kk)=ERGAS(I_GT,I_F,ratio);
         elseif any(strcmpi(qindex_list{jj},{'Q_avg','UIQI','Q_{avg}'}))
            MatrixResults(jj,kk)=Q(I_GT,I_F,DynamicRange);
        elseif any(strcmpi(qindex_list{jj},{'SAM'}))
            MatrixResults(jj,kk)=SAM(I_GT,I_F);
        elseif any(strcmpi(qindex_list{jj},{'Q2n_alt','Q2^n_alt','Q2n_{alt}','Q2^n_{alt}'}))
            cd(q2nnew_folder);
            MatrixResults(jj,kk)=Q2n_picone(I_F,I_GT,Qblocks_size,1,0);
            cd(current_folder);
        elseif strncmpi(qindex_list{jj},'SCC_',4)
            MatrixResults(jj,kk)=SCC(I_GT,I_F,edge{jj},'global1');
        elseif strncmpi(qindex_list{jj},'SCCold_',7) || strncmpi(qindex_list{jj},'SCCo_',5)
            MatrixResults(jj,kk)=SCC(I_GT,I_F,edge{jj});
        elseif strncmpi(qindex_list{jj},'Q2n_b',5) || strncmpi(qindex_list{jj},'Q2^n_b',6)
            I_GT_temp=I_GT(cut(jj)+1:end-cut(jj),cut(jj)+1:end-cut(jj),:);
            I_F_temp=I_F(cut(jj)+1:end-cut(jj),cut(jj)+1:end-cut(jj),:);
            MatrixResults(jj,kk)=q2n(I_GT_temp,I_F_temp,blocknum(jj),shiftnum(jj));
        elseif strncmpi(qindex_list{jj},'Q2n_alt_b',9) || strncmpi(qindex_list{jj},'Q2^n_alt_b',10) || ...
                strncmpi(qindex_list{jj},'Q2n_{alt}_b',11) || strncmpi(qindex_list{jj},'Q2^n_{alt}_b',12)
            cd(q2nnew_folder);
            I_GT_temp=I_GT(cut(jj)+1:end-cut(jj),cut(jj)+1:end-cut(jj),:);
            I_F_temp=I_F(cut(jj)+1:end-cut(jj),cut(jj)+1:end-cut(jj),:);
            MatrixResults(jj,kk)=q2n_new(I_GT_temp,I_F_temp,blocknum(jj),shiftnum(jj));
            cd(current_folder);
        elseif strncmpi(qindex_list{jj},'ssim',4)
            MatrixResults_temp=0;
            for ii=1:sizes_F(3)
                %MatrixResults_temp=MatrixResults_temp+ssim(I_F(:,:,ii),I_GT(:,:,ii),'DynamicRange',DynamicRange);
                MatrixResults_temp=MatrixResults_temp+ssim(I_F(:,:,ii),I_GT(:,:,ii),[0.01,0.03],fspecial('gaussian', 11, 1.5),DynamicRange);
            end
            MatrixResults(jj,kk)=MatrixResults_temp/sizes_F(3);
        end
    end
    fprintf([repmat('\b',[1,7]),'%3d/%3d'],kk,N_tests);
end
fprintf('. Done!\n');


for jj=1:N_index
    if strncmpi(qindex_list{jj},'D_',2) || strncmpi(qindex_list{jj},'SAM',3) || strncmpi(qindex_list{jj},'ERGAS',5) || strncmpi(qindex_list{jj},'MSE',3)
        MatrixResults_current=MatrixResults(jj,:);
        MatrixResults_current(isnan(MatrixResults_current))=Inf;        
        [~,MatrixResults_current]=sort(MatrixResults_current,'ascend');
        [~,MatrixResults_order(jj,:)]=sort(MatrixResults_current,'ascend');
    else
        MatrixResults_current=MatrixResults(jj,:);
        MatrixResults_current(isnan(MatrixResults_current))=0;
        [~,MatrixResults_current]=sort(MatrixResults_current,'descend');
        [~,MatrixResults_order(jj,:)]=sort(MatrixResults_current,'ascend');
    end
end

MatrixResults=reshape(MatrixResults,[N_index,sizes_F(4:end)]);
MatrixResults_order=reshape(MatrixResults_order,[N_index,sizes_F(4:end)]);

end