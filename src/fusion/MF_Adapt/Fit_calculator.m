function [fit_individuals,fit_vett]=Fit_calculator(I_Fus,I_MS,I_PAN,PAN_LP,L,ratio,fit_method,sensor,Qblocks_size,I_GT)
% keyboard
% cd ..\Quality_Indices
for idim=1:size(I_Fus,3),
    I_Fus_LP(:,:,idim) = imresize(imresize(I_Fus(:,:,idim),1/ratio),ratio);
    %         imageHR_LP(:,:,idim) = imresize(imresize(imageHR(:,:,idim),1/Resize_fact),size(imageHR(:,:,idim)));
    
end

%%% Estimate alpha
%Con PAN_LP
PAN_LP_estim=PAN_LP;  
%Con Imresize
% PAN_LP_estim=imresize(I_PAN,1/ratio);  
% if  (2^round(log2(ratio)) ~= ratio)
%    PAN_LP_estim = imresize(PAN_LP_estim,[size(I_PAN,1),size(I_PAN,1)]);
% else
%     PAN_LP_estim = interp23tap(PAN_LP_estim,ratio);
% end

alpha(1,1,:) = estimation_alpha(cat(3,I_MS,ones(size(I_MS,1),size(I_MS,2))),PAN_LP_estim,'global');
equivPAN=sum(cat(3,I_MS,ones(size(I_MS,1),size(I_MS,2))) .* repmat(alpha,[size(I_MS,1) size(I_MS,2) 1]),3);

imageHR_LP = imresize(imresize(I_PAN,1/ratio),ratio);
[I_Fus_LR,dum]=resize_images(I_Fus,I_PAN,ratio,sensor);
if  (2^round(log2(ratio)) ~= ratio)
    I_Fus_LR = imresize(I_Fus_LR,[size(I_Fus,1),size(I_Fus,1)]);
else
    I_Fus_LR = interp23tap(I_Fus_LR,ratio);
end

if isempty(I_GT)
    I_GT=zeros(size(I_Fus));
else
    dim_cut=11;
    I_Fus=I_Fus(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    I_Fus_LP=I_Fus_LP(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    I_Fus_LR=I_Fus_LR(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    I_MS=I_MS(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    I_PAN=I_PAN(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    imageHR_LP=imageHR_LP(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    
    I_GT=I_GT(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
end
% keyboard
fit_vett=[];
% Qblocks_size=30;
% if strcmp(fit_method,'Qavg_All') || strcmp(fit_method,'QPan_All')...
%         || strcmp(fit_method,'QNR_All')|| strcmp(fit_method,'Q2n_All')...
%         || strcmp(fit_method,'Q2nPan_All') || strcmp(fit_method,'Q2nDegraded_All') ...
%         || strcmp(fit_method,'QavgDegraded_All') || strcmp(fit_method,'QavgDegraded_All')
% keyboard
 if   strcmp(fit_method(end-2:end),'All')
%         keyboard
    %calculate 'Q'
    fit_vett(1) = Q(I_GT,I_Fus,2^L);
    
    %calculate 'Q2n'
    %     I_GT = I_GT(11:end-11,11:end-11,:);
    %     I_Fus = I_Fus(11:end-11,11:end-11,:);
    %     I_Fus(I_Fus > 2^L) = 2^L;
    %     I_Fus(I_Fus < 0) = 0;
    fit_vett(2)= q2n(I_GT,I_Fus,Qblocks_size,Qblocks_size);
    
    %calculate 'QPan'
    t=0.5;
    fit_vett(3) = t*Q(I_MS,I_Fus_LP,2^L)+...
        (1-t)*Q(I_Fus-I_Fus_LP,repmat(I_PAN-imageHR_LP,[1 1 size(I_Fus,3)]),2^L);
    
    %calculate 'Q2nPan'
    t=0.5;
    fit_vett(4) = t*q2n(I_MS,I_Fus_LP,Qblocks_size,Qblocks_size)+...
        (1-t)*q2n(I_Fus-I_Fus_LP,repmat(I_PAN-imageHR_LP,[1 1 size(I_Fus,3)]),Qblocks_size,Qblocks_size);
    
    %calculate 'QNRI', 'D_lambda', 'D_S' in this order
    a1=floor(size(I_Fus,1)/Qblocks_size)*Qblocks_size;
    a2=floor(size(I_Fus,2)/Qblocks_size)*Qblocks_size;
    [fit_vett(5),fit_vett(6),fit_vett(7)] = QNR(I_Fus(1:a1,1:a2,:),I_MS(1:a1,1:a2,:),I_PAN(1:a1,1:a2,:),sensor,ratio);
    
    %calculate'Q2nDegraded'
    %     keyboard
    if sum(abs(I_MS(:)-I_Fus(:)))~=0
        fit_vett(8) = q2n(I_MS,I_Fus_LR,Qblocks_size,Qblocks_size);
    else
        fit_vett(8) = 0;
    end
    
    %calculate'QavgDegraded'
    if sum(abs(I_MS(:)-I_Fus(:)))~=0
        fit_vett(9) = Q(I_MS,I_Fus_LR,2^L);
    else
        fit_vett(9) = 0;
    end
    %calculate'QEquivPAN'
            
            if sum(abs(I_MS(:)-I_Fus(:)))~=0
                fit_vett(10) = Q(PAN_LP,equivPAN,2^L);
            else
                fit_vett(10) = 0;
            end

    switch fit_method
        case 'Qavg_All'
            fit_individuals= fit_vett(1);
        case 'Q2n_All'
            fit_individuals= fit_vett(2);
        case 'QPan_All'
            fit_individuals= fit_vett(3);
        case 'Q2nPan_All'
            fit_individuals= fit_vett(4);
        case 'QNR_All'
            fit_individuals= fit_vett(5);
        case 'Q2nDegraded_All'
            fit_individuals= fit_vett(8);
        case 'QavgDegraded_All'
            fit_individuals= fit_vett(9);
            case 'QEquivPAN_All'
              fit_individuals= fit_vett(10);  
    end
else
    switch fit_method
        case 'Qavg'
            fit_individuals = Q(I_GT,I_Fus,2^L);
            
        case 'Q2n'
            %             I_GT = I_GT(11:end-11,11:end-11,:);
            %             I_Fus = I_Fus(11:end-11,11:end-11,:);
            %             I_Fus(I_Fus > 2^L) = 2^L;
            %             I_Fus(I_Fus < 0) = 0;
            fit_individuals= q2n(I_GT,I_Fus,Qblocks_size,Qblocks_size);
        case 'QPan'
            
            fit_individuals =  t*Q(I_MS,I_Fus_LP,2^L)+...
                (1-t)*Q(I_Fus-I_Fus_LP,repmat(I_PAN-imageHR_LP,[1 1 size(I_Fus,3)]),2^L);
        case 'Q2nPan'
            t=0.5;
            
            fit_individuals =  t*q2n(I_MS,I_Fus_LP,Qblocks_size,Qblocks_size)+...
                (1-t)*q2n(I_Fus-I_Fus_LP,repmat(I_PAN-imageHR_LP,[1 1 size(I_Fus,3)]),Qblocks_size,Qblocks_size);
        case 'QNR'
            %             a1=floor(size(I_Fus,1)/Qblocks_size)*Qblocks_size;
            %             a2=floor(size(I_Fus,2)/Qblocks_size)*Qblocks_size;
            %             fit_individuals =QNR(I_Fus(1:a1,1:a2,:),imageLR(1:a1,1:a2,:),imageHR(1:a1,1:a2,:),sensor,Resize_fact);
            fit_individuals =QNR(I_Fus,I_MS,I_PAN,sensor,ratio);
        case 'Q2nDegraded'
            %             keyboard
            if sum(abs(I_MS(:)-I_Fus(:)))~=0
                try
                    fit_individuals = q2n(I_MS,I_Fus_LR,Qblocks_size,Qblocks_size);
                catch
                    keyboard
                end
            else
                fit_individuals = 0;
            end
        case 'QavgDegraded'
            if sum(abs(I_MS(:)-I_Fus(:)))~=0
                fit_individuals = Q(I_MS,I_Fus_LR,2^L);
            else
                fit_individuals = 0;
            end
        
        case 'QEquivPAN'
            if sum(abs(I_MS(:)-I_Fus(:)))~=0
                fit_individuals = Q(PAN_LP,equivPAN,2^L);
            else
                fit_individuals = 0;
            end
            
    end
end

% keyboard

% cd ..\MFAdapt
