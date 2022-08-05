
%% GSlocal ideal
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

for kkk=1:11
    S = Sbpt_MS(:,:,kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_manual(kkk), ergas_GS_manual(kkk), q2n_GS_manual(kkk)] = evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_manual(kkk), D_S_GS_manual(kkk), QNR_GS_manual(kkk)] = evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

% for kkk=1:11
%     S = Sbpt_MS(:,:,kkk);
%     equaliz = 0;
%     interp = 0;
%     tic
%     [PS_GS,Coeff] = GS_local_ideal_mean(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
%     time_GS = toc
%     if strcmp(filename(end-1),'R')
%         [sam_GS_manual_mean(kkk), ergas_GS_manual_mean(kkk), q2n_GS_manual_mean(kkk)] = evaluate_fusion(HRMS, PS_GS);
%     else
%         [D_lambda_GS_manual_mean(kkk), D_S_GS_manual_mean(kkk), QNR_GS_manual_mean(kkk)] = evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
%     end
% end


Blockvett=[300 150 100 75 50 25 15 10 5];
for kkk=1:length(Blockvett)
    S = Blockvett(kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,alpha_opt,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_manual_block(kkk), ergas_GS_manual_block(kkk), q2n_GS_manual_block(kkk)] = ...
            evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_manual_block(kkk), D_S_GS_manual_block(kkk), QNR_GS_manual_block(kkk)] = ...
            evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

for kkk=1:length(Blockvett)
    S = Blockvett(kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local_ideal_mean(MS, PAN, I_LP_input,alpha_opt,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_manual_block_mean(kkk), ergas_GS_manual_block_mean(kkk), q2n_GS_manual_block_mean(kkk)] = ...
            evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_manual_block_mean(kkk), D_S_GS_manual_block_mean(kkk), QNR_GS_manual_block_mean(kkk)] = ...
            evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

%%
h2=figure;
 subplot(221)
plot(nregs,q2n_GS_manual,'b-','Linewidth',2)
hold on
plot(nregs,q2n_GS_manual_mean,'b--','Linewidth',2)
plot(90000./(Blockvett.^2),q2n_GS_manual_block,'r-','Linewidth',2)
plot(90000./(Blockvett.^2),q2n_GS_manual_block_mean,'r--','Linewidth',2)
xlim([1 max(nregs)])
legend('Segmentation','Segmentation Mean','Blocks','Blocks Mean',4)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)

% title('Q^{2n}')

%  figure,
 subplot(222)
plot(nregs,ergas_GS_manual,'Linewidth',2)
hold on
plot(nregs,ergas_GS_manual_mean,'b--','Linewidth',2)
plot(90000./(Blockvett.^2),ergas_GS_manual_block,'r','Linewidth',2)
plot(90000./(Blockvett.^2),ergas_GS_manual_block_mean,'r--','Linewidth',2)
xlim([1 max(nregs)])
legend('Segmentation','Segmentation Mean','Blocks','Blocks Mean')
ylabel('ERGAS','FontSize', 12)

xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
% title('ERGAS')

%  figure,
 subplot(223)
plot(nregs,sam_GS_manual,'Linewidth',2)
hold on
plot(nregs,sam_GS_manual_mean,'b--','Linewidth',2)
plot(90000./(Blockvett.^2),sam_GS_manual_block,'r','Linewidth',2)
plot(90000./(Blockvett.^2),sam_GS_manual_block_mean,'r--','Linewidth',2)
xlim([1 max(nregs)])
legend('Segmentation','Segmentation Mean','Blocks','Blocks Mean')
% title('SAM')

ylabel('SAM','FontSize', 12)
xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
if strcmp(filename,'China_RR')
    saveas(h2,'China_IndicesVsNoParts.fig','fig')
    saveas(h2,'China_IndicesVsNoParts.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_IndicesVsNoParts.fig','fig')
    saveas(h2,'Rome_IndicesVsNoParts.eps','psc2')
end




