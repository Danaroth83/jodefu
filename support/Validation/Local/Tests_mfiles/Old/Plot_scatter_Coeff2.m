% clear MSE
% clear ERGAS_vett
% clear MSE_imag
clear std_vett
clear std_cumul
clear var_vett
clear var_cumul
% clear alpha_mean
% clear std_vett_Coeff
% S_Band = S;
S = Sbpt_MS(:,:,3);
% S = Sbpt_PAN(:,:,3);
labels = unique(S);
std_cumul=0;
var_cumul=0;
for kk=1:length(labels)
    idx = S==labels(kk);
    a1=[];
    for iband=1:size(MS,3)
        a=alpha_opt(:,:,iband);
        a1=[a1 a(idx)];
    end
    % amean=mean(a1,1);
    std_vett(kk)=mean(std(a1,0,1));
    std_cumul=std_cumul+std_vett(kk)*numel(a1)/numel(alpha_opt);
    % var_vett(kk)=mean(var(a1,0,1));
    % var_cumul=var_cumul+var_vett(kk)*numel(a1)/numel(alpha_opt);
    alpha_mean(kk,:)=mean(a1,1);
    % l=(a1-repmat(alpha_mean,size(a1,1),1));
    
    l=(a1-repmat(alpha_mean(kk,:),size(a1,1),1));
    clear norml
    for il=1:size(l,1)
        norml(il,:)=norm(l(il,:));
    end
    mean_norm_2=mean(norml.^2);
    var_vett(kk)=mean_norm_2;
    var_cumul=var_cumul+var_vett(kk)*numel(a1)/numel(alpha_opt);
    
end
std_vett;
alpha_mean
% std_cumul
var_vett
var_cumul
Figu=zeros(size(PAN));
Figu(S==length(labels))=PAN(S==length(labels));
figure,imagesc(Figu)
% mean(std_vett)
%%

clear std_vett_blocks
clear std_cumul_blocks
clear var_vett_blocks
clear var_cumul_blocks
% clear std_vett_Coeff
nblocks=[1 2 5 10 20 50 100 300];
BlockSize=100;
% for kk=1:length(nblocks)
%     BlockSize=nblocks(kk);
kk=0;
std_cumul_blocks=0;
var_cumul_blocks=0;
clear var_vett_blocks
S=zeros(size(alpha_opt,1),size(alpha_opt,2));
for y=1 : BlockSize : size(alpha_opt,1)
    for x=1 : BlockSize : size(alpha_opt,2)
        kk=kk+1;
        %         startx = max(x - floor(BlockSize/2), 1);
        %         starty = max(y - floor(BlockSize/2), 1);
        %         endy = min(y + floor(BlockSize/2), size(MS,1));
        %         endx = min(x + floor(BlockSize/2), size(MS,2));
        startx = max(x, 1);
        starty = max(y, 1);
        endy = min(y + BlockSize-1, size(MS,1));
        endx = min(x + BlockSize-1, size(MS,2));
        S(starty:endy,startx:endx)=kk;
        a1=[];
        for iband=1:size(MS,3)
            a=alpha_opt(starty:endy,startx:endx,iband);
            a1=[a1 a(:)];
        end
        %         amean=mean(a1,1);
        %         keyboard
        %         std_vett_blocks(kk)=mean(std(a1,0,1));
        %         std_cumul_blocks=std_cumul_blocks+std_vett_blocks(kk)*numel(a1)/numel(alpha_opt);
        %             var_vett_blocks(kk)=mean(var(a1,0,1));
        %         var_cumul_blocks=var_cumul_blocks+var_vett_blocks(kk)*numel(a1)/numel(alpha_opt);
        l=(a1-repmat(amean,size(a1,1),1));
        clear norml
        for il=1:size(l,1)
            norml(il,:)=norm(l(il,:));
        end
        mean_norm_2=mean(norml.^2);
        var_vett_blocks(kk)=mean_norm_2;
        var_cumul_blocks=var_cumul_blocks+var_vett_blocks(kk)*(numel(a1))/numel(alpha_opt);
        
    end
end

Figu=zeros(size(PAN));
Figu(S==2)=PAN(S==2);
figure,imagesc(Figu)



% std_vett_blocks
% std_cumul_blocks
var_vett_blocks
var_cumul_blocks


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% MSE Coeff sull'immagine %%%%%%%%

clear err_Coeff
clear err_Coeff_cumul
clear err_Image
clear err_Image_cumul
% S = Sbpt_MS(:,:,5);
% S = Sbpt_PAN(:,:,3);
labels = unique(S);
err_Coeff_cumul=zeros(1,size(alpha_opt,3));
err_Image_cumul=zeros(1,size(alpha_opt,3));
for kk=1:length(labels)
    idx = S==labels(kk);
    a1=[];
    b1=[];
    aa1=[];
    bb1=[];
    for iband=1:size(MS,3)
        a=Coeff(:,:,iband);
        aa=alpha_opt(:,:,iband);
        b=PS_GS(:,:,iband);
        bb=HRMS(:,:,iband);
        a1=[a1 a(idx)];
        aa1=[aa1 aa(idx)];
        b1=[b1 b(idx)];
        bb1=[bb1 bb(idx)];
    end
    err_Coeff(kk,:)=sum((a1-aa1).^2,1)/size(a1,1);
    err_Coeff_cumul=err_Coeff_cumul+err_Coeff(kk,:)*size(a1,1)/size(alpha_opt,1)/size(alpha_opt,2);
    err_Image(kk,:)=sum((b1-bb1).^2,1)/size(b1,1);
    err_Image_cumul=err_Image_cumul+err_Image(kk,:)*size(b1,1)/size(alpha_opt,1)/size(alpha_opt,2);
    
end
% sqrt(err_Coeff_cumul)
% sqrt(err_Image_cumul)
% sqrt(err_Coeff)
% sqrt(err_Image)
% mean(err_Coeff,2)
mean(err_Image,2)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% MSE Coeff sull'immagine %%%%%%%%

clear err_Coeff_blocks
clear err_Coeff_cumul_blocks
clear err_Image_blocks
clear err_Image_cumul_blocks
% % clear std_vett_Coeff
% nblocks=[1 2 5 10 20 50 100 300];
BlockSize=S;
% for kk=1:length(nblocks)
%     BlockSize=nblocks(kk);

kk=0;
err_Coeff_cumul_blocks=zeros(1,size(alpha_opt,3));
err_Image_cumul_blocks=zeros(1,size(alpha_opt,3));
for y=1 : BlockSize : size(alpha_opt,1)
    for x=1 : BlockSize : size(alpha_opt,2)
        kk=kk+1;
        %         startx = max(x - floor(BlockSize/2), 1);
        %         starty = max(y - floor(BlockSize/2), 1);
        %         endy = min(y + floor(BlockSize/2), size(MS,1));
        %         endx = min(x + floor(BlockSize/2), size(MS,2));
        startx = max(x, 1);
        starty = max(y, 1);
        endy = min(y + BlockSize-1, size(MS,1));
        endx = min(x + BlockSize-1, size(MS,2));
        
        a1=[];
        b1=[];
        aa1=[];
        bb1=[];
        for iband=1:size(MS,3)
            a=Coeff(starty:endy,startx:endx,iband);
            aa=alpha_opt(starty:endy,startx:endx,iband);
            b=PS_GS(starty:endy,startx:endx,iband);
            bb=HRMS(starty:endy,startx:endx,iband);
            a1=[a1 a(:)];
            aa1=[aa1 aa(:)];
            b1=[b1 b(:)];
            bb1=[bb1 bb(:)];
        end
        err_Coeff_blocks(kk,:)=sum((a1-aa1).^2,1)/size(a1,1);
        err_Coeff_cumul_blocks=err_Coeff_cumul_blocks+err_Coeff_blocks(kk,:)...
            *size(a1,1)/size(alpha_opt,1)/size(alpha_opt,2);
        err_Image_blocks(kk,:)=sum((b1-bb1).^2,1)/size(b1,1);
        err_Image_cumul_blocks=err_Image_cumul_blocks+err_Image_blocks(kk,:)...
            *size(b1,1)/size(alpha_opt,1)/size(alpha_opt,2);
        
    end
end
% sqrt(err_Coeff_blocks)
% sqrt(err_Coeff_cumul_blocks)
% sqrt(err_Image_blocks)
% sqrt(err_Image_cumul_blocks)
% mean(err_Coeff_blocks,2)
mean(err_Image_blocks,2)
% Figu=zeros(size(PAN));
% Figu(S==2)=PAN(S==2);
% figure,imagesc(Figu)


%% Disegno scatter coeff
S = Sbpt_MS(:,:,3);
flag_plot=0;
for ii = 1:size(MS,3),
    MS_Band = squeeze(MS(:,:,ii));
    I_LP_Band = squeeze(I_LP_input(:,:));
    [b,bint,r,rint,stats] = regress(I_LP_Band(:),[MS_Band(:),ones(size(MS_Band(:),1),1)]);
    R2_tot(ii)=stats(1);
    if flag_plot
        figure,plot(MS_Band(:),I_LP_Band(:),'o')
        xlim([0 500])
        ylim([0 1000])
        hold on
        xax=linspace(0,500,1000);
        yax=polyval(b,xax);
        plot(xax,yax)
        
        title(['R^2 = ',num2str(R2_tot)])
    end
    labels = unique(S);
    for kk=1:length(labels)
        idx = S==labels(kk);
        MS_Band_label = MS_Band(idx);
        I_LP_Band_label = I_LP_Band(idx);
        [b,bint,r,rint,stats] = regress(I_LP_Band_label,[MS_Band_label,ones(size(MS_Band_label,1),1)]);
        R2_segm(kk,ii)=stats(1);
        if flag_plot
            figure,plot(MS_Band_label(:),I_LP_Band_label(:),'o')
            xlim([0 500])
            ylim([0 1000])
            hold on
            xax=linspace(0,500,1000);
            yax=polyval(b,xax);
            plot(xax,yax)
            title(['R^2 = ',num2str(stats(1))])
        end
    end
    
end

R2_tot
R2_segm


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

for kkk=1:16
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

Blockvett=[300 150 100 75 50 25 15 10 5 2 1];
for kkk=1:length(Blockvett)
    S = Blockvett(kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_manual_block(kkk), ergas_GS_manual_block(kkk), q2n_GS_manual_block(kkk)] = ...
            evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_manual_block(kkk), D_S_GS_manual_block(kkk), QNR_GS_manual_block(kkk)] = ...
            evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

%%
h2=figure;
 subplot(221)
plot(nregs(1:end-1),q2n_GS_manual(1:end-1),'Linewidth',2)
hold on
plot(90000./(Blockvett.^2),q2n_GS_manual_block,'r','Linewidth',2)
% xlim([1 max(nregs)])
legend('Segmentation','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)

% title('Q^{2n}')

%  figure,
 subplot(222)
plot(nregs(1:end-1),ergas_GS_manual(1:end-1),'Linewidth',2)
hold on
plot(90000./(Blockvett.^2),ergas_GS_manual_block,'r','Linewidth',2)
% xlim([1 max(nregs)])
legend('Segmentation','Blocks')
ylabel('ERGAS','FontSize', 12)

xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
% title('ERGAS')

%  figure,
 subplot(223)
plot(nregs(1:end-1),sam_GS_manual(1:end-1),'Linewidth',2)
hold on
plot(90000./(Blockvett.^2),sam_GS_manual_block,'r','Linewidth',2)
xlim([1 max(nregs)])
legend('Segmentation','Blocks')
% title('SAM')

ylabel('SAM','FontSize', 12)
xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
if strcmp(filename,'China_RR')
    saveas(h2,'China_IndicesVsNoParts_PixelWise.fig','fig')
    saveas(h2,'China_IndicesVsNoParts_PixelWise.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_IndicesVsNoParts_PixelWise.fig','fig')
    saveas(h2,'Rome_IndicesVsNoParts_PixelWise.eps','psc2')
end

%%
h2=figure;
 subplot(221)
semilogx(nregs(1:end-1),q2n_GS_manual(1:end-1),'Linewidth',2)
hold on
semilogx(90000./(Blockvett.^2),q2n_GS_manual_block,'r','Linewidth',2)
% xlim([1 max(nregs)])
legend('Segmentation','Blocks',2)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)

% title('Q^{2n}')

%  figure,
 subplot(222)
semilogx(nregs(1:end-1),ergas_GS_manual(1:end-1),'Linewidth',2)
hold on
semilogx(90000./(Blockvett.^2),ergas_GS_manual_block,'r','Linewidth',2)
% xlim([1 max(nregs)])
legend('Segmentation','Blocks',3)
ylabel('ERGAS','FontSize', 12)

xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
% title('ERGAS')

%  figure,
 subplot(223)
semilogx(nregs(1:end-1),sam_GS_manual(1:end-1),'Linewidth',2)
hold on
semilogx(90000./(Blockvett.^2),sam_GS_manual_block,'r','Linewidth',2)
xlim([1 max(nregs)])
legend('Segmentation','Blocks',3)
% title('SAM')

ylabel('SAM','FontSize', 12)
xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
if strcmp(filename,'China_RR')
    saveas(h2,'China_IndicesVsNoPartslog_PixelWise.fig','fig')
    saveas(h2,'China_IndicesVsNoPartslog_PixelWise.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_IndicesVsNoPartslog_PixelWise.fig','fig')
    saveas(h2,'Rome_IndicesVsNoPartslog_PixelWise.eps','psc2')
end

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

for kkk=1:size(Sbpt_MS,3)
    S = Sbpt_MS(:,:,kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
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
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
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
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
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
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_PS_kmeans_PAN(kkk), ergas_GS_PS_kmeans_PAN(kkk), q2n_GS_PS_kmeans_PAN(kkk)] = evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_PS_kmeans_PAN(kkk), D_S_GS_PS_kmeans_PAN(kkk), QNR_GS_PS_kmeans_PAN(kkk)] = evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

Blockvett=[300 150 100 75 50 25 15 10 5 2 1];
for kkk=1:length(Blockvett)
    S = Blockvett(kkk);
    equaliz = 0;
    interp = 0;
    tic
    [PS_GS,Coeff] = GS_local_ideal(MS, PAN, I_LP_input,HRMS,S,equaliz,interp);
    time_GS = toc
    if strcmp(filename(end-1),'R')
        [sam_GS_manual_block(kkk), ergas_GS_manual_block(kkk), q2n_GS_manual_block(kkk)] = ...
            evaluate_fusion(HRMS, PS_GS);
    else
        [D_lambda_GS_manual_block(kkk), D_S_GS_manual_block(kkk), QNR_GS_manual_block(kkk)] = ...
            evaluate_fusion_FR(MS,PAN, PS_GS,sensor,Resize_fact);
    end
end

%%
h2=figure;
%  subplot(221)
plot(nregs(1:end),q2n_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
plot(nregs(1:end),q2n_GS_PS_BPT_PAN,'b--','Linewidth',2)
plot(kregs(1:end),q2n_GS_PS_kmeans_MS,'r-','Linewidth',2)
plot(kregs(1:end),q2n_GS_PS_kmeans_PAN,'r--','Linewidth',2)
plot(90000./(Blockvett.^2),q2n_GS_manual_block,'g','Linewidth',2)
% xlim([1 max(nregs)])
legend('BPT_MS','BPT_PAN','Kmeans_MS','Kmeans_PAN','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)

% % title('Q^{2n}')
% 
% %  figure,
%  subplot(222)
% plot(nregs(1:end-1),ergas_GS_manual,'Linewidth',2)
% hold on
% plot(90000./(Blockvett.^2),ergas_GS_manual_block,'r','Linewidth',2)
% % xlim([1 max(nregs)])
% legend('Segmentation','Blocks')
% ylabel('ERGAS','FontSize', 12)
% 
% xlabel('Number of parts','FontSize', 12)
% set(gca, 'FontSize', 12)
% % title('ERGAS')
% 
% %  figure,
%  subplot(223)
% plot(nregs(1:end-1),sam_GS_manual,'Linewidth',2)
% hold on
% plot(90000./(Blockvett.^2),sam_GS_manual_block,'r','Linewidth',2)
% xlim([1 max(nregs)])
% legend('Segmentation','Blocks')
% % title('SAM')
% 
% ylabel('SAM','FontSize', 12)
% xlabel('Number of parts','FontSize', 12)
% set(gca, 'FontSize', 12)
if strcmp(filename,'China_RR')
    saveas(h2,'China_IndicesVsNoParts.fig','fig')
    saveas(h2,'China_IndicesVsNoParts.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_IndicesVsNoParts.fig','fig')
    saveas(h2,'Rome_IndicesVsNoParts.eps','psc2')
end

%%
h2=figure;
 subplot(221)
semilogx(nregs(1:end-1),q2n_GS_manual,'Linewidth',2)
hold on
semilogx(90000./(Blockvett.^2),q2n_GS_manual_block,'r','Linewidth',2)
% xlim([1 max(nregs)])
legend('Segmentation','Blocks',4)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)

% title('Q^{2n}')

%  figure,
 subplot(222)
semilogx(nregs(1:end-1),ergas_GS_manual,'Linewidth',2)
hold on
semilogx(90000./(Blockvett.^2),ergas_GS_manual_block,'r','Linewidth',2)
% xlim([1 max(nregs)])
legend('Segmentation','Blocks')
ylabel('ERGAS','FontSize', 12)

xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
% title('ERGAS')

%  figure,
 subplot(223)
semilogx(nregs(1:end-1),sam_GS_manual,'Linewidth',2)
hold on
semilogx(90000./(Blockvett.^2),sam_GS_manual_block,'r','Linewidth',2)
xlim([1 max(nregs)])
legend('Segmentation','Blocks')
% title('SAM')

ylabel('SAM','FontSize', 12)
xlabel('Number of parts','FontSize', 12)
set(gca, 'FontSize', 12)
if strcmp(filename,'China_RR')
    saveas(h2,'China_IndicesVsNoPartslog.fig','fig')
    saveas(h2,'China_IndicesVsNoPartslog.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_IndicesVsNoPartslog.fig','fig')
    saveas(h2,'Rome_IndicesVsNoPartslog.eps','psc2')
end
%%
h2=figure;
%  subplot(221)
semilogx(nregs(1:end),q2n_GS_PS_BPT_MS,'b-','Linewidth',2)
hold on
semilogx(nregs(1:end),q2n_GS_PS_BPT_PAN,'b--','Linewidth',2)
semilogx(kregs(1:end),q2n_GS_PS_kmeans_MS,'r-','Linewidth',2)
semilogx(kregs(1:end),q2n_GS_PS_kmeans_PAN,'r--','Linewidth',2)
semilogx(90000./(Blockvett.^2),q2n_GS_manual_block,'g','Linewidth',2)
xlim([1 200])
legend('BPT_{MS}','BPT_{PAN}','Kmeans_{MS}','Kmeans_{PAN}','Blocks',2)
xlabel('Number of parts','FontSize', 12)
ylabel('Q^{2n}','FontSize', 12)
set(gca, 'FontSize', 12)
if strcmp(filename,'China_RR')
    saveas(h2,'China_IndicesVsNoPartslog_All.fig','fig')
    saveas(h2,'China_IndicesVsNoPartslog_All.eps','psc2')
elseif strcmp(filename,'Rome_RR')
    saveas(h2,'Rome_IndicesVsNoPartslog_All.fig','fig')
    saveas(h2,'Rome_IndicesVsNoPartslog_All.eps','psc2')
end
