% clear MSE
% clear ERGAS_vett
% clear MSE_imag
clear std_vett
clear std_vett_Coeff
S_Band = S;
% S_Band = Sbpt_MS(:,:,6);

labels = unique(S_Band);
std_cumul=0;
std_cumul_Coeff=0;
for kk=1:length(labels)

S_Band_rep=repmat(S_Band,[1 1 8]);
idx = S_Band_rep==labels(kk);
a1=alpha_opt(idx);
c1=Coeff(idx);
std_vett(kk)=std(a1(:))/std(alpha_opt(:));
std_vett_Coeff(kk)=std(c1(:))/std(Coeff(:));
std_cumul=std_cumul+std_vett(kk)*numel(a1)/numel(alpha_opt);
std_cumul_Coeff=std_cumul_Coeff+std_vett_Coeff(kk)*numel(c1)/numel(Coeff);
end


% std_vett
% std_cumul
% max_std_vett=max(std_vett)
% 
% std_vett_Coeff
% std_cumul_Coeff
% max_std_vett_Coeff=max(std_vett_Coeff)


% for kk=1:length(labels)
% 
% S_Band_rep=repmat(S_Band,[1 1 8]);
% idx = S_Band_rep==labels(kk);
% c1=Coeff(idx);
% sum(idx(:))
% a1=alpha_opt(idx);
% figure,plot(a1(:),c1(:),'o')
% MSE(kk)=mean((a1-c1).^2)/mean((c1).^2);
% end
% 
% 
% MSE

% for kk=1:length(labels)
% idx = S_Band==labels(kk);
% idx_rep=repmat(idx,[1 1 8]);
% Err=HRMS(idx_rep)-I_Fus(idx_rep);
% % GT=HRMS(idx_rep);
% % ERGAS_vett(kk)=mean(Err(:).^2)/(mean(GT))^2; 
% % ERGAS_vett(kk)=ERGAS(I_GT,I_Fus,Resize_fact);
% MSE_imag(kk)=mean(Err(:).^2);
% end
% MSE_imag