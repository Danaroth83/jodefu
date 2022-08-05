
%% GSlocal Comparison
%%% GS1:
% I_LP_input = genPAN_LP(PAN, 'mean', Resize_fact);
%%% GS2:
I_LP_input = genPAN_LP(PAN, 'ATWT', Resize_fact);
% I_LP_input = genPAN_LP(PAN, 'GLP', Resize_fact);
% h = genPANMTF(Resize_fact, sensor); I_LP_input = genPAN_LP(PAN, 'GLP_MTF', Resize_fact, h);
% h = genMTF(4, sensor); I_LP_input = genPAN_LP(PAN, 'GLP_MTF', Resize_fact, h);
%%% GSA: Use GS_local with S=ones(size(PAN));
% PAN_LP = genPAN_LP(PAN, 'ATWT', Resize_fact);
% % PAN_LP = genPAN_LP(PAN, 'GLP', Resize_fact);
% % h = genPANMTF(Resize_fact, sensor); PAN_LP = genPAN_LP(PAN, 'GLP_MTF', Resize_fact, h);
% % h = genMTF(4, sensor); PAN_LP = genPAN_LP(PAN, 'GLP_MTF', Resize_fact, h);
% I_LP_input = genI_LP(MS, 'adaptive', PAN_LP, 'segm', ones(size(PAN)));

%%%%%%%%  Segmentation choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% S = Sbpt_MS(:,:,10);
% S = Sbpt_PAN(:,:,11);

% S = Skmeans_MS(:,:,3);
% S = Skmeans_PAN(:,:,3);
% S=100;
% S=100;
clear sam_GS_manual
clear ergas_GS_manual 
clear q2n_GS_manual
clear am_GS_manual_block
clear ergas_GS_manual_block
clear q2n_GS_manual_block

for kkk=1:size(Sbpt_MS,3)
    S = Sbpt_MS(:,:,kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local(MS, PAN, I_LP_input,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_PS_BPT_MS(kkk), ergas_GS_PS_BPT_MS(kkk), q2n_GS_PS_BPT_MS(kkk)] = evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_PS_BPT_MS(kkk), D_S_GS_PS_BPT_MS(kkk), QNR_GS_PS_BPT_MS(kkk)] = evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

for kkk=1:size(Sbpt_PAN,3)
    S = Sbpt_PAN(:,:,kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local(MS, PAN, I_LP_input,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_PS_BPT_PAN(kkk), ergas_GS_PS_BPT_PAN(kkk), q2n_GS_PS_BPT_PAN(kkk)] = evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_PS_BPT_PAN(kkk), D_S_GS_PS_BPT_PAN(kkk), QNR_GS_PS_BPT_PAN(kkk)] = evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

for kkk=1:size(Skmeans_MS,3)
    S = Skmeans_MS(:,:,kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local(MS, PAN, I_LP_input,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_PS_kmeans_MS(kkk), ergas_GS_PS_kmeans_MS(kkk), q2n_GS_PS_kmeans_MS(kkk)] = evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_PS_kmeans_MS(kkk), D_S_GS_PS_kmeans_MS(kkk), QNR_GS_PS_kmeans_MS(kkk)] = evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

for kkk=1:size(Skmeans_PAN,3)
    S = Skmeans_PAN(:,:,kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local(MS, PAN, I_LP_input,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_PS_kmeans_PAN(kkk), ergas_GS_PS_kmeans_PAN(kkk), q2n_GS_PS_kmeans_PAN(kkk)] = evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_PS_kmeans_PAN(kkk), D_S_GS_PS_kmeans_PAN(kkk), QNR_GS_PS_kmeans_PAN(kkk)] = evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

Blockvett=[300 150 100 75 50 25 15 10 5];
for kkk=1:length(Blockvett)
    S = Blockvett(kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local(MS, PAN, I_LP_input,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_block(kkk), ergas_GS_block(kkk), q2n_GS_block(kkk)] = ...
            evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_block(kkk), D_S_GS_block(kkk), QNR_GS_block(kkk)] = ...
            evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

%%
if strcmp(filename(end-1),'R')

h2=figure;
 subplot(221)
plot(nregs(1:end),q2n_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
plot(nregs(1:end),q2n_GS_PS_BPT_PAN,'b--','Linewidth',2)
plot(kregs(1:end),q2n_GS_PS_kmeans_MS,'r-','Linewidth',2)
plot(kregs(1:end),q2n_GS_PS_kmeans_PAN,'r--','Linewidth',2)
plot(90000./(Blockvett.^2),q2n_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
% legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)


 subplot(222)
plot(nregs(1:end),ergas_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
plot(nregs(1:end),ergas_GS_PS_BPT_PAN,'b--','Linewidth',2)
plot(kregs(1:end),ergas_GS_PS_kmeans_MS,'r-','Linewidth',2)
plot(kregs(1:end),ergas_GS_PS_kmeans_PAN,'r--','Linewidth',2)
plot(90000./(Blockvett.^2),ergas_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
% legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('ERGAS','FontSize', 12)
set(gca, 'FontSize', 12)

 subplot(223)
plot(nregs(1:end),sam_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
plot(nregs(1:end),sam_GS_PS_BPT_PAN,'b--','Linewidth',2)
plot(kregs(1:end),sam_GS_PS_kmeans_MS,'r-','Linewidth',2)
plot(kregs(1:end),sam_GS_PS_kmeans_PAN,'r--','Linewidth',2)
plot(90000./(Blockvett.^2),sam_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
hleg1=legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks')
set(hleg1, 'Position', [.67,.24,.1,.1]);
xlabel('Number of parts','FontSize', 12)
ylabel('SAM','FontSize', 12)
set(gca, 'FontSize', 12)

if strcmp(filename,'China_RR')
    saveas(h2,'China_Comparison_Segmentations.fig','fig')
    saveas(h2,'China_Comparison_Segmentations.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_Comparison_Segmentations.fig','fig')
    saveas(h2,'Rome_Comparison_Segmentations.eps','psc2')
end

%%
h2=figure;
 subplot(221)
semilogy(nregs(1:end),q2n_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
semilogy(nregs(1:end),q2n_GS_PS_BPT_PAN,'b--','Linewidth',2)
semilogy(kregs(1:end),q2n_GS_PS_kmeans_MS,'r-','Linewidth',2)
semilogy(kregs(1:end),q2n_GS_PS_kmeans_PAN,'r--','Linewidth',2)
semilogy(90000./(Blockvett.^2),q2n_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
% legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)


 subplot(222)
semilogy(nregs(1:end),ergas_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
semilogy(nregs(1:end),ergas_GS_PS_BPT_PAN,'b--','Linewidth',2)
semilogy(kregs(1:end),ergas_GS_PS_kmeans_MS,'r-','Linewidth',2)
semilogy(kregs(1:end),ergas_GS_PS_kmeans_PAN,'r--','Linewidth',2)
semilogy(90000./(Blockvett.^2),ergas_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
% legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('ERGAS','FontSize', 12)
set(gca, 'FontSize', 12)

subplot(223)
semilogy(nregs(1:end),sam_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
semilogy(nregs(1:end),sam_GS_PS_BPT_PAN,'b--','Linewidth',2)
semilogy(kregs(1:end),sam_GS_PS_kmeans_MS,'r-','Linewidth',2)
semilogy(kregs(1:end),sam_GS_PS_kmeans_PAN,'r--','Linewidth',2)
semilogy(90000./(Blockvett.^2),sam_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
hleg1=legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks')
set(hleg1, 'Position', [.67,.24,.1,.1]);
xlabel('Number of parts','FontSize', 12)
ylabel('SAM','FontSize', 12)
set(gca, 'FontSize', 12)

if strcmp(filename,'China_RR')
    saveas(h2,'China_Comparison_Segmentations.fig','fig')
    saveas(h2,'China_Comparison_Segmentations.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_Comparison_Segmentations.fig','fig')
    saveas(h2,'Rome_Comparison_Segmentations.eps','psc2')
end
else
    h2=figure;
 subplot(221)
plot(nregs(1:end),QNR_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
plot(nregs(1:end),QNR_GS_PS_BPT_PAN,'b--','Linewidth',2)
plot(kregs(1:end),QNR_GS_PS_kmeans_MS,'r-','Linewidth',2)
plot(kregs(1:end),QNR_GS_PS_kmeans_PAN,'r--','Linewidth',2)
plot(90000./(Blockvett.^2),QNR_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
% legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('QNR','FontSize', 12)
set(gca, 'FontSize', 12)


 subplot(222)
plot(nregs(1:end),D_lambda_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
plot(nregs(1:end),D_lambda_GS_PS_BPT_PAN,'b--','Linewidth',2)
plot(kregs(1:end),D_lambda_GS_PS_kmeans_MS,'r-','Linewidth',2)
plot(kregs(1:end),D_lambda_GS_PS_kmeans_PAN,'r--','Linewidth',2)
plot(90000./(Blockvett.^2),D_lambda_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
% legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('D_{\lambda}','FontSize', 12)
set(gca, 'FontSize', 12)

 subplot(223)
plot(nregs(1:end),D_S_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
plot(nregs(1:end),D_S_GS_PS_BPT_PAN,'b--','Linewidth',2)
plot(kregs(1:end),D_S_GS_PS_kmeans_MS,'r-','Linewidth',2)
plot(kregs(1:end),D_S_GS_PS_kmeans_PAN,'r--','Linewidth',2)
plot(90000./(Blockvett.^2),D_S_GS_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
xlim([1 200])
hleg1=legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks')
set(hleg1, 'Position', [.67,.24,.1,.1]);
xlabel('Number of parts','FontSize', 12)
ylabel('D_S','FontSize', 12)
set(gca, 'FontSize', 12)

if strcmp(filename,'Toulouse_RR')
    saveas(h2,'Toulouse_Comparison_Segmentations.fig','fig')
    saveas(h2,'Toulouse_Comparison_Segmentations.eps','psc2')
elseif strcmp(filename,'Rio_FR')
    saveas(h2,'Rio_Comparison_Segmentations.fig','fig')
    saveas(h2,'Rio_Comparison_Segmentations.eps','psc2')
end
end