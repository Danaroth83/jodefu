%% Figures Segmentation
ind_regs=5;
nparts=nregs(ind_regs);
titfig=[filename,'_Segmentations_',num2str(nparts)]

h2=figure;
subplot(221)
imagesc(Sbpt_MS(:,:,ind_regs))
% colorbar
title('BPT MS')
axis off

subplot(222)
imagesc(Sbpt_PAN(:,:,ind_regs))
% colorbar
title('BPT PAN')
axis off
subplot(223)
imagesc(Skmeans_MS(:,:,ind_regs))
% colorbar
title('K-means MS')
axis off
subplot(224)
imagesc(Skmeans_PAN(:,:,ind_regs))
% colorbar
title('K-means PAN')
axis off

% 
% 
% hleg1=legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks')
% set(hleg1, 'Position', [.67,.24,.1,.1]);
% xlabel('Number of parts','FontSize', 12)
% ylabel('SAM','FontSize', 12)
% set(gca, 'FontSize', 12)
 saveas(h2,titfig,'fig')
  saveas(h2,[titfig,'.eps'],'psc2')
% if strcmp(filename,'China_RR')
%     saveas(h2,'China_Segmentations_5.fig','fig')
%     saveas(h2,'China_Segmentations_5.eps','psc2')
% elseif strcmp(filename,'Rome_RR')
%     saveas(h2,'Rome_Segmentations_5.fig','fig')
%     saveas(h2,'Rome_Segmentations__5.eps','psc2')
% end

% if strcmp(filename,'China_RR')
%     saveas(h2,'China_Segmentations_5.fig','fig')
%     saveas(h2,'China_Segmentations_5.eps','psc2')
% elseif strcmp(filename,'Rome_RR')
%     saveas(h2,'Rome_Segmentations_5.fig','fig')
%     saveas(h2,'Rome_Segmentations__5.eps','psc2')
% end
