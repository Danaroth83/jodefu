
%%Morphological pyramid options
diff_oper=0;
diff_strel=0;
dim_strel=3;


% Sampling kind
% sampl_kind='mean';
sampl_kind='downs';

% int_meth='nearest';
int_meth='bilinear';
% int_meth='bicubic';
% int_meth='lanczos2';
% int_meth='lanczos3';
% int_meth='interp23tap';

% Random Hill Climbing restarts
num_restart = 1; 
% num_restart = 2^(dim_strel^2+2); %=all

% Random Hill Climbing iterations
num_max_iter = 1; % 1 per MC

% Chromosome length
if diff_oper==1,
    len_crom = 11; 
else
    len_crom=5;
end
if diff_strel==0
    len_crom_strel=dim_strel^2+2;
else
    len_crom_strel=3*(dim_strel^2+2);
end
len_tot=len_crom+len_crom_strel;

% Monte Carlo trials
% num_popul=2;
%num_popul=2^len_tot;%all
num_popul=2^len_crom_strel;%all strels
%%%%%% Fit type
% fit_method='Q2n';
% fit_method='QPan';
% fit_method='Q2nPan';
% fit_method='Qavg';
% fit_method='QNR';
% fit_method='Q2nDegraded';
% fit_method='QavgDegraded';
% fit_method='QEquivPAN';

% fit_method='Q2n_All'; 
% fit_method='QPan_All';
% fit_method='Q2nPan_All';
% fit_method='Qavg_All';
% fit_method='QNR_All';
% fit_method='Q2nDegraded_All';
% fit_method='QavgDegraded_All';
fit_method='QEquivPAN_All';

%% Fusion method
% fusion_method='HPF';
% fusion_method='SDM';
fusion_method='HPM'; 
% fusion_method='GS'; 
% fusion_method='Contrast';

%% 
operators = [0 0 0];
% 000: gradient; 010: top-hat; 100: toggle; 110: Bai (toggle+top-hat)
% 001: LOCO; 011:median

