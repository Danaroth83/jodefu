%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Reduced resolution quality indexes. 
% 
% Interface:
%           [Q_index, SAM_index, ERGAS_index, sCC, Q2n_index] = indexes_evaluation(I_F,I_GT,ratio,L,Q_blocks_size,flag_cut_bounds,dim_cut,th_values)
%
% Inputs:
%           I_F:                Fused Image;
%           I_GT:               Ground-Truth image;
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value;
%           L:                  Image radiometric resolution; 
%           Q_blocks_size:      Block size of the Q-index locally applied;
%           flag_cut_bounds:    Cut the boundaries of the viewed Panchromatic image;
%           dim_cut:            Define the dimension of the boundary cut;
%           th_values:          Flag. If th_values == 1, apply an hard threshold to the dynamic range.
%
% Outputs:
%           Q_index:            Q index;
%           SAM_index:          Spectral Angle Mapper (SAM) index;
%           ERGAS_index:        Erreur Relative Globale Adimensionnelle de Synth�se (ERGAS) index;
%           sCC:                spatial Correlation Coefficient between fused and ground-truth images;
%           Q2n_index:          Q2n index.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatrixResults,qindex_list] = indexes_evaluation_mod3(I_F,GT,varargin)


current_folder=pwd;
q2nnew_folder=fullfile(fileparts(mfilename('fullpath')),'Q2n_new');


%% Handling input variables in a single cell
if iscell(varargin{1})
    varargin1=varargin{1};
else
    varargin1=varargin;
end


%% Handling multiple tests in I_F
N_tests=size(I_F,4);
if N_tests>1
    MatrixResults=[];
    for ii=1:N_tests
        [MR_temp,qindex_list]=indexes_evaluation_mod3(I_F(:,:,:,1),GT,varargin1);
        I_F(:,:,:,1)=[];
        MatrixResults=cat(1,MatrixResults,MR_temp);
    end
    return;
end

%% Parameters implementation

if isfield(GT,'data'), I_GT=GT.data; else, I_GT=GT; end
if isfield(GT,'ratio'), ratio=GT.ratio; else, ratio=4; end
if isfield(GT,'IntBits'), L=GT.IntBits; else, L=11; end
if isfield(GT,'DynamicRange'), DynamicRange=GT.DynamicRange; else, DynamicRange=2^L-1; end
if isfield(GT,'Qblocks_size'), Qblocks_size=GT.Qblocks_size; else, Qblocks_size=32; end

qindex_list={'Q2^n','SAM','ERGAS','SCC','UIQI'};
flag_cutbounds=0;
dim_cut=0;
th_values=0;

flag_L=0;
flag_ratio=0;

for ii=1:2:numel(varargin1)
    pname=varargin1{ii};
    pval=varargin1{ii+1};   
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
    end
end

if flag_cutbounds
    I_GT = I_GT(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    I_F = I_F(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
end
    
if th_values
    I_F(I_F > DynamicRange) = DynamicRange;
    I_F(I_F < 0) = 0;
    if flag_L==0, warning('Bit per pixel not setup for thresolding'); end
end


MatrixResults=zeros(1,numel(qindex_list));
for jj=1:numel(qindex_list)
    if any(strcmpi(qindex_list{jj},{'Q2n','Q2^n','Q'}))
        qindex_list{jj}='Q2^n';
        MatrixResults(jj)=q2n(I_GT,I_F,Qblocks_size,Qblocks_size);
    elseif any(strcmpi(qindex_list{jj},{'SCC','CC'}))
        qindex_list{jj}='sCC';
        MatrixResults(jj)=SCC(I_GT,I_F,'laplacian','global1');
    elseif any(strcmpi(qindex_list{jj},{'PSNR'}))
        qindex_list{jj}='PSNR';
        MatrixResults(jj)=PSNR(I_GT,I_F,L);
    elseif any(strcmpi(qindex_list{jj},{'PSNRmax','PSNRm'}))
        qindex_list{jj}='PSNRm';
        MatrixResults(jj)=PSNR(I_GT,I_F,[]);
    elseif any(strcmpi(qindex_list{jj},{'MSE'}))
        qindex_list{jj}='MSE';
        MatrixResults(jj)=MSE(I_GT,I_F);
    elseif any(strcmpi(qindex_list{jj},{'ERGAS'}))
        qindex_list{jj}='ERGAS';
        MatrixResults(jj)=ERGAS(I_GT,I_F,ratio);
        if flag_ratio==0, warning('Scale ratio for ERGAS non inputted, automatically set as 4'); end
     elseif any(strcmpi(qindex_list{jj},{'Q_avg','UIQI','Q_{avg}'}))
        qindex_list{jj}='Q_{avg}';
        MatrixResults(jj)=Q(I_GT,I_F,DynamicRange);
    elseif any(strcmpi(qindex_list{jj},{'SAM'}))
        qindex_list{jj}='SAM';
        MatrixResults(jj)=SAM(I_GT,I_F);
    elseif any(strcmpi(qindex_list{jj},{'Q2n_alt','Q2^n_alt','Q2n_{alt}','Q2^n_{alt}'}))
        qindex_list{jj}='Q2^n_{alt}';
        cd(q2nnew_folder);
        MatrixResults(jj)=Q2n_picone(I_F,I_GT,Qblocks_size,1,0);
        cd(current_folder);
    elseif strncmpi(qindex_list{jj},'SCC_',4)
        idx_unds=find(qindex_list{jj},'_');
        qindex_list_temp=qindex_list{jj};
        edge=qindex_list_temp(idx_unds+1:end);
        if any(strcmpi(edge,{'laplacian','Lap','{Lap}','{Laplacian}'}))
            edge='laplacian'; edge2='{Lap}';
        elseif any(strcmpi(edge,{'prewitt','pre','{pre}','{prewitt}'}))
            edge='prewitt'; edge2='{Pre}';
        elseif any(strcmpi(edge,{'prewitt2','pre2','{pre2}','{pr2}','{prewitt2}'}))
            edge='prewitt2'; edge2='{Pr2}';
        elseif any(strcmpi(edge,{'sobel','sob','{sobel}','{sob}'}))
            edge='sobel'; edge2='{Sob}';
        elseif any(strcmpi(edge,{'sobel2','sob2','{sobel2}','{sob2}','{so2}'}))
            edge='sobel2'; edge2='{So2}';
        end
        qindex_list{jj}=['sCC_',edge2];
        MatrixResults(jj)=SCC(I_GT,I_F,edge,'global1');
    elseif strncmpi(qindex_list{jj},'SCCold_',7) || strncmpi(qindex_list{jj},'SCCo_',5)
        idx_unds=find(qindex_list{jj},'_');
        qindex_list_temp=qindex_list{jj};
        edge=qindex_list_temp(idx_unds+1:end);
        if any(strcmpi(edge,{'laplacian','Lap','{Lap}','{Laplacian}'}))
            edge='laplacian'; edge2='{Lap}';
        elseif any(strcmpi(edge,{'prewitt','pre','{pre}','{prewitt}'}))
            edge='prewitt'; edge2='{Pre}';
        elseif any(strcmpi(edge,{'prewitt2','pre2','{pre2}','{pr2}','{prewitt2}'}))
            edge='prewitt2'; edge2='{Pr2}';
        elseif any(strcmpi(edge,{'sobel','sob','{sobel}','{sob}'}))
            edge='sobel'; edge2='{Sob}';
        elseif any(strcmpi(edge,{'sobel2','sob2','{sobel2}','{sob2}','{so2}'}))
            edge='sobel2'; edge2='{So2}';
        end
        qindex_list{jj}=['sCCo_',edge2];
        MatrixResults(jj)=SCC(I_GT,I_F,edge);
    elseif strncmpi(qindex_list{jj},'Q2n_b',5) || strncmpi(qindex_list{jj},'Q2^n_b',6)
        qindex_list{jj}=strrep(qindex_list{jj},'^','');
        qindex_list_temp=qindex_list{jj};
        idx_c=find(qindex_list_temp=='c');
        if isempty(idx_c), cut=0; else, cut=qindex_list_temp(idx_c+1:end); idx_c=length(qindex_list_temp)+1; end
        idx_s=find(qindex_list_temp=='s');
        blocknum=str2double(qindex_list_temp(6:idx_s-1));
        shiftnum=str2double(qindex_list_temp(idx_s+1:idx_c-1));
        I_GT_temp=I_GT(cut+1:end-cut,cut+1:end-cut,:);
        I_F_temp=I_F(cut+1:end-cut,cut+1:end-cut,:);
        MatrixResults(jj)=q2n(I_GT_temp,I_F_temp,blocknum,shiftnum);
        if cut==0
            qindex_list{jj}=sprintf('Q2^n_{b%d}^{s%d}',blocknum,shiftnum);
        else
            qindex_list{jj}=sprintf('Q2^n_{b%d}^{s%dc%d}',blocknum,shiftnum,cut);
        end
    elseif strncmpi(qindex_list{jj},'Q2n_alt_b',9) || strncmpi(qindex_list{jj},'Q2^n_alt_b',10) || ...
            strncmpi(qindex_list{jj},'Q2n_{alt}_b',11) || strncmpi(qindex_list{jj},'Q2^n_{alt}_b',12)
        qindex_list{jj}=strrep(strrep(strrep(qindex_list{jj},'^',''),'{',''),'}','');
        qindex_list_temp=qindex_list{jj};
        idx_s=find(qindex_list_temp=='s');
        idx_c=find(qindex_list_temp=='c');
        if isempty(idx_c), cut=0; else, cut=qindex_list_temp(idx_c+1:end); idx_c=length(qindex_list_temp)+1; end
        blocknum=str2double(qindex_list_temp(9:idx_s-1));
        shiftnum=str2double(qindex_list_temp(idx_s+1:idx_c-1));
        cd(q2nnew_folder);
        I_GT_temp=I_GT(cut+1:end-cut,cut+1:end-cut,:);
        I_F_temp=I_F(cut+1:end-cut,cut+1:end-cut,:);
        MatrixResults(jj)=q2n_new(I_GT_temp,I_F_temp,blocknum,shiftnum);
        cd(current_folder);
        if cut==0
            qindex_list{jj}=sprintf('Q2^n_{altb%d}^{s%d}',blocknum,shiftnum);
        else
            qindex_list{jj}=sprintf('Q2^n_{altb%d}^{s%dc%d}',blocknum,shiftnum,cut);
        end
    elseif strncmpi(qindex_list{jj},'ssim',4)
        for ii=1:size(I_GT,3)
            %MatrixResults_temp=MatrixResults_temp+ssim(I_F(:,:,ii),I_GT(:,:,ii),'DynamicRange',DynamicRange);
            MatrixResults_temp=MatrixResults_temp+ssim(I_F(:,:,ii),I_GT(:,:,ii),[0.01,0.03],fspecial('gaussian', 11, 1.5),DynamicRange);
        end
        MatrixResults(jj)=MatrixResults_temp/size(I_GT,3);
        qindex_list{jj}='SSIM';
    end
end

end