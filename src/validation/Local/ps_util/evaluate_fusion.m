function [sam_v, ergas_v, q2n_v] = evaluate_fusion(I_ref, I_fus)


%% Quality Index Blocks
Qblocks_size = 32;

%% Cut Final Image
flag_cut_bounds = 1;
dim_cut = 11;

%% Threshold values out of dynamic range
th_values = 0;

%% Resize Factor
Resize_fact = 4;

%% Radiometric Resolution
L = 11;


if flag_cut_bounds
    I_ref = I_ref(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    I_fus = I_fus(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
end

if th_values
    I_fus(I_fus > 2^L) = 2^L;
    I_fus(I_fus < 0) = 0;
end

sam_v = SAM(I_ref, I_fus);

ergas_v = ERGAS(I_ref, I_fus, Resize_fact);


q2n_v = q2n(I_ref, I_fus, Qblocks_size, Qblocks_size);

fprintf('SAM:%.3f\tERGAS:%.3f\tQ2n:%.3f\n', sam_v, ergas_v, q2n_v);