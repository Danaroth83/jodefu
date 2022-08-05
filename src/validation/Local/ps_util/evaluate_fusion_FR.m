function [D_lambda,D_S,QNR_index] = evaluate_fusion_FR(I_MS,I_PAN, I_fus, sensor, Resize_fact)
% [D_lambda,D_S,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,sensor,tag,ratio)


[QNR_index,D_lambda,D_S]= QNR(I_fus,I_MS,I_PAN,sensor,Resize_fact);

fprintf('QNR:%.3f\tDlambda:%.3f\tDS:%.3f\n', QNR_index,D_lambda,D_S);