function [I_Fus,fit,D,P_LP] = MF_Adapt_gradLoad(I_PAN,I_MS,L,ratio,sensor,I_GT)
% keyboard
fprintf('Adaptive Morphological Operator Test\n')

imageLR = double(I_MS);
imageHR = double(I_PAN);
if nargin == 6,
    I_GT = double(I_GT);
else
    I_GT = [];
end

%%%%%%%%%%%%%Multiple Band PAN%%%%%%%%%%%%%%%%%%
imageHR = repmat(imageHR,[1 1 size(imageLR,3)]);
%%Equalization
for ii = 1 : size(imageLR,3)
    imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(imageLR(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(imageLR(:,:,ii));
end
% %%%%% Cut negative values 
% imageHR=max(imageHR,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Parameters settings
MF_Adapt_Settings
%%%cromosomes (len_tot)
%[0 1 1 . .  0|0 0 . . . .  1];
%[  len_crom  |len_crom_strel]
%[ oper & lev |   text_se    ]

load chromosomes_file 
% keyboard
%% Best Result
[maxEl,ind_pos] = max(individuals(:,1));

ii = 1;

fit=[];%{individuals,best_textse,fit_vett};


%keyboard
%% Optimal pyramid setting display
% Number of levels
% if individuals(ii,end-1) == 0 && individuals(ii,end) == 0
%     lev = 1;
% elseif individuals(ii,end-1) == 0 && individuals(ii,end) == 1
%     lev = 2;
% elseif individuals(ii,end-1) == 1 && individuals(ii,end) == 0
%     lev = 3;
% else
%     lev = 4;
% end
 lev=ceil(log2(ratio))+1;
    
fprintf('Number of levels: %1.0f\n',lev)
% keyboard
% Operators To be checked
if diff_oper==1,
    for ilev=2:lev,
        fprintf(['Livello ',num2str(ilev-1),' operator choice: '])
        if (operators(1)==0) && (operators(2)==1) && (operators(3)==0)
            fprintf('Top Hat');
        elseif (operators(1)==1) && (operators(2)==0) && (operators(2)==0)
            fprintf('Toggle');
        elseif (operators(1)==1) && (operators(2)==1) && (operators(2)==0)
            fprintf('Bai (Toggle+Top Hat');
        elseif (operators(1)==0) && (operators(2)==0) && (operators(2)==1)
            fprintf('LOCO');
        elseif (operators(1)==0) && (operators(2)==1) && (operators(2)==1)
            fprintf('Median');
        elseif (operators(1)==0) && (operators(2)==0) && (operators(2)==0)
            fprintf('Gradient');
        else
            fprintf('Error');
        fprintf('\n')
        end
    end    %
else
    fprintf(['Tutti i livelli, operator : '])
    if (operators(1)==0) && (operators(2)==1) && (operators(3)==0)
        fprintf('Top Hat');
    elseif (operators(1)==1) && (operators(2)==0) && (operators(2)==0)
        fprintf('Toggle');
    elseif (operators(1)==1) && (operators(2)==1) && (operators(2)==0)
        fprintf('Bai (Toggle+Top Hat');
    elseif (operators(1)==0) && (operators(2)==0) && (operators(2)==1)
        fprintf('LOCO');
    elseif (operators(1)==0) && (operators(2)==1) && (operators(2)==1)
        fprintf('Median');
    elseif (operators(1)==0) && (operators(2)==0) && (operators(2)==0)
        fprintf('Gradient');
    else
        fprintf('Error');
    end
    fprintf('\n')
end

fprintf('\n')


% Structruring Element
% textse=best_textse(ii,:);
if diff_strel==1,
    for ilev=2:lev,
        textse_lin_lev=textse((ilev-2)*(dim_strel^2+2)+1:(ilev-1)*(dim_strel^2+2));
        [textse_lev]=Strel_Constr(textse_lin_lev(1:end-2),textse_lin_lev(end-1)*2+textse_lin_lev(end)+2);
        
        fprintf(['Livello ',num2str(ilev-1),', strel: '])
        fprintf('\n')
        %         disp(reshape(textse((ilev-2)*dim_strel^2+1:(ilev-1)*dim_strel^2),dim_strel,dim_strel))
        disp(textse_lev)
        fprintf('\n')
    end
else
    textse_lin_lev=textse(1:dim_strel^2+2);
    [textse_lev]=Strel_Constr_whole(textse_lin_lev(1:end-2));
    
    fprintf('Tutti i livelli, strel: ')
    fprintf('\n')
    disp(textse_lev)
    fprintf('\n')
    dim_strel1 = textse_lin_lev(end-1)*2+textse_lin_lev(end)+2;
    fprintf('Dimensione strel:  %1.0f\n',dim_strel1)
    fprintf('Tutti i livelli, strel resized: ')
    fprintf('\n')
    strel_eff=imresize(textse_lev,[dim_strel1,dim_strel1]);
    disp(strel_eff)
end



%% Best Image Construction

%%% Operators
% if diff_oper==1,
%     %%% Different Operators
%     operators=individuals(ii,2:end-2);
% else
%     %%Same operators
%     operators=repmat(individuals(ii,2:end-2),1,3);
% end

% P = Pyr_Dec_Whole(imageHR,textse,lev,operators);
P = Pyr_Dec_Whole_grad(imageHR,textse,lev,operators,sampl_kind,int_meth);
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Fusion    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(P,4)<lev
    keyboard
end
P_LP=P(:,:,:,lev);
D=squeeze(P(:,:,:,1) - P_LP);
switch fusion_method
    case 'HPF'
        I_Fus = imageLR+ D;
    case 'SDM'
        PLP = imresize(imresize(P(:,:,:,1),1/4),4);
        I_Fus = imageLR + (imageLR ./ (PLP+eps)) .* D;
    case 'HPM'
%         keyboard
        I_Fus = imageLR .* (P(:,:,:,1)./(P_LP+eps));
%         eps0=10^-3;
%         I_Fus = imageLR .* (P(:,:,:,1)./(P_LP+1));
        case 'GS'
              %%% Detail Extraction
              %%% Coefficients
              g = ones(1,size(I_MS,3));
              for ii = 1 : size(I_MS,3)
                  h = imageLR(:,:,ii);
                  h2 = P_LP(:,:,ii);
                  c = cov(h2(:),h(:));
                  g(ii) = c(1,2)/var(h2(:));
              end
              I_Fus = zeros(size(imageLR));
              for ii = 1 : size(imageLR,3)
                  I_Fus(:,:,ii) = imageLR(:,:,ii) + D(:,:,ii) .* g(ii);
              end
    case 'Contrast'
%         keyboard
        P_1=squeeze(P(:,:,:,1));
        textse_lev = Strel_Constr_whole(textse_lin_lev(1:end-2));
%         for ii=1:lev-1
        P_num=imdilate(P_1,textse_lev);
%         end
        P_den=imdilate(P_1,textse_lev)+imerode(P_1,textse_lev);
        I_Fus = imageLR .* (2*P_num./(P_den+eps));

end


% keyboard

%%%%%%%%%%%%%%%%%%%%%%% Final Equalization

% for ii = 1 : size(imageLR,3)
%   I_Fus(:,:,ii) = I_Fus(:,:,ii) - mean2(squeeze(I_Fus(:,:,ii))) + mean2(squeeze(imageLR(:,:,ii)));
% end


% for ii = 1 : size(imageLR,3)
%   I_Fus(:,:,ii) = (I_Fus(:,:,ii) - mean2(squeeze(I_Fus(:,:,ii)))).*(std2(squeeze(imageLR(:,:,ii)))./std2(squeeze(I_Fus(:,:,ii)))) + mean2(squeeze(imageLR(:,:,ii)));
% end

% for b = 1 : size(imageLR,3)
%
%     histogram_matrix=zeros(MbnHist,'uint16');
%     for hc_y = 1 : size(imageLR,1)
%         for hc_x = 1 : size(imageLR,2)
%             value = uint16(imageLR(hc_y,hc_x,b));
%             histogram_matrix(value) = histogram_matrix(value) + 1;
%         end
%     end
%
%     I_Fus(:,:,b) = spechist(I_Fus(:,:,b),histogram_matrix,2^16);
%
% end


