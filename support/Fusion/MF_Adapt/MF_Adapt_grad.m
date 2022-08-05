function [I_Fus,fit,D,P_LP] = MF_Adapt_grad(I_PAN,I_MS,L,ratio,sensor,Qblocks_size,I_GT)
% keyboard
fprintf('Adaptive Morphological Operator Test Gradient Sampling\n')
addpath('mf_util')

imageLR = double(I_MS);
imageHR = double(I_PAN);
if nargin == 7,
    I_GT = double(I_GT);
else
    I_GT = [];
end
% keyboard
%%%%%%%%%%%%%Multiple Band PAN%%%%%%%%%%%%%%%%%%
imageHR = repmat(imageHR,[1 1 size(imageLR,3)]);
%%Equalization
for ii = 1 : size(imageLR,3)
    imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(imageLR(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(imageLR(:,:,ii));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Parameters settings
MF_Adapt_Settings
%%%cromosomes (len_tot)
%[0 1 1 . .  0|0 0 . . . .  1];
%[  len_crom  |len_crom_strel]
%[ oper & lev |   text_se    ]

%%% Parte da cambiare
%%%individuals = operators & levels
%[0  | 1 1 . .  0 | 0 0  ];
%[fit| operators  |levels]
% num operators=3 if diff_strel=0
% num operators=9 if diff_strel=1

% text_se (len_crom_strel)
%[1 1 . . . . . 0 |    0 0     ];
%[    strel       | dim_strel  ];
%[    strel       | dim_strel  ];
% len_crom_strel=dim_strel^2+2 if diff_strel=0
% len_crom_strel=3*dim_strel^2 if diff_strel=1

%%% Chromosome generation
% Seed initialization
reset(RandStream.getGlobalStream,sum(100*clock));

% indexes=randperm(2^len_tot);
% keyboard
indexes=1:2^len_crom_strel;
% indexes=randperm(2^len_crom_strel);
% indexes=randi(2^len_crom_strel,num_popul,1);
chromosomes=dec2gc(indexes(1:num_popul)-1,len_tot);
individuals = [zeros(num_popul,1) chromosomes(:,1:len_crom)];
fit_individuals = zeros(size(individuals,1),1);
fit_vett = zeros(size(individuals,1),10);
% keyboard
if diff_strel==1,
    best_textse=zeros(size(individuals,1),len_crom_strel);
else
    best_textse=zeros(size(individuals,1),3*len_crom_strel);
end
%keyboard
for ii = 1 : size(individuals,1)
    
    tic
    %     keyboard
    % Levels
    if individuals(ii,end-1) == 0 && individuals(ii,end) == 0
        lev = 1;
    elseif individuals(ii,end-1) == 0 && individuals(ii,end) == 1
        lev = 2;
    elseif individuals(ii,end-1) == 1 && individuals(ii,end) == 0
        lev = 3;
    else
        lev = 4;
    end
    %     lev =3; %2 decomposition levels
    lev=floor(log2(ratio))+1;
    %     %%% Operators
    %         if diff_oper==1,
    %             %%% Different Operators
    %             operators=individuals(ii,2:end-2);
    %         else
    %             %%Same operators
    %             operators=repmat(individuals(ii,2:end-2),1,3);
    %
    %         end
    
    best_ind = [];
    best_fit = -Inf;
    best_fit_run = -Inf;
    best_ind_run = [];
    
    for rr = 1 : num_restart
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%% Pyramid construction    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Restart
        reset(RandStream.getGlobalStream,sum(100*clock));
        if diff_strel==1
            textse = chromosomes(ii,len_crom+1:end);
        else
            textse = repmat(chromosomes(ii,len_crom+1:end),1,3);
        end
        
        for iter = 1 : num_max_iter
            %             keyboard
            % P = Pyr_Dec_Whole(imageHR,textse,lev,operators);
            P = Pyr_Dec_Whole_grad(imageHR,textse,lev,operators,sampl_kind,int_meth);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%% Fusion    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %   P_LP= squeeze(P(:,:,:,lev));
            %   for i1=1:size(P_LP,3)
            %       P_LP(:,:,i1)=P_LP(:,:,i1)-mean2(P_LP(:,:,i1))+mean2(I_MS(:,:,i1));
            %       %     P_LP(:,:,i1)=(P_LP(:,:,i1)-mean2(P_LP(:,:,i1)))*std2(I_MS(:,:,i1))/std2(P_LP(:,:,i1))+mean2(I_MS(:,:,i1));
            %   end
            P2= Pyr_Dec_Whole_grad(I_PAN,textse,lev,operators,sampl_kind,int_meth);
            if size(P,4)<lev
                keyboard
            end
            P_LP=P(:,:,:,lev);
            PAN_LP=P2(:,:,:,lev);
            D=squeeze(P(:,:,:,1) - P_LP);
            switch fusion_method
                case 'HPF'
                    I_Fus = imageLR+ D;
                    if sum(isfinite(I_Fus(:)))~=numel(I_Fus)
                        % keyboard
                        break
                    end
                    
                case 'SDM'
                    PLP = imresize(imresize(P(:,:,:,1),1/4),4);
                    I_Fus = imageLR + (imageLR ./ (PLP+eps)) .* D;
                    if sum(isfinite(I_Fus(:)))~=numel(I_Fus)
                        break
                    end
                case 'HPM'
                    I_Fus = imageLR .* (P(:,:,:,1)./(P_LP+eps));
                    if sum(isfinite(I_Fus(:)))~=numel(I_Fus)
                        break
                    end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%% Final Equalization
            
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
            
            
            %             keyboard
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%% Quality evaluation and update %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             keyboard
            if strcmp(fit_method(end-2:end),'All')
                [fit_individuals(ii),fit_vett(ii,:)]=Fit_calculator(I_Fus,I_MS,I_PAN,PAN_LP,L,ratio,fit_method,sensor,Qblocks_size,I_GT);
            else
                [fit_individuals(ii),dum2]=Fit_calculator(I_Fus,I_MS,I_PAN,PAN_LP,L,ratio,fit_method,sensor,Qblocks_size,I_GT);
            end
            if fit_individuals(ii) >= best_fit
                best_fit = fit_individuals(ii);
                best_ind = textse;
            end
            textse = best_ind;
            %             if isempty(textse)
            %             keyboard
            %             end
            % Move in the neighbour
            if diff_strel==1,
                pos_flip = round(rand().*(len_crom_strel-1)) + 1;
                textse(pos_flip) = 1 - textse(pos_flip);
            else
                pos_flip = round(rand().*(len_crom_strel-1)) + 1;
                textse(pos_flip) = 1 - textse(pos_flip);
                textse=repmat(textse(1:len_crom_strel),1,3);
            end
            
        end
        if best_fit > best_fit_run
            best_fit_run = best_fit;
            best_ind_run = best_ind;
        end
        
    end
    if ~isempty(best_ind_run),
        best_textse(ii,:)=best_ind_run;
    else
        best_textse(ii,:)=zeros(size(textse));
    end
    fit_individuals(ii)=best_fit_run;
    
    if rem(ii,10) == 0,
        fprintf('Iterazione numero: %d, Tempo iterazione corrente: %3.2f \n',ii,toc),
    end
    if rem(ii,10000) == 0,
        save Morph_Test_Temporary
    end
end

% New matrices including fitness as first column
individuals = [fit_individuals, individuals(:,2:end)];


%% Best Result
[maxEl,ind_pos] = max(individuals(:,1));

ii = ind_pos;

fit={individuals,best_textse,fit_vett};

% keyboard
%% Optimal pyramid setting display
% Number of levels
if individuals(ii,end-1) == 0 && individuals(ii,end) == 0
    lev = 1;
elseif individuals(ii,end-1) == 0 && individuals(ii,end) == 1
    lev = 2;
elseif individuals(ii,end-1) == 1 && individuals(ii,end) == 0
    lev = 3;
else
    lev = 4;
end
% lev=3;
lev=floor(log2(ratio))+1;
fprintf('Number of levels: %1.0f\n',lev)

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
%
% Structruring Element
textse=best_textse(ii,:);
if diff_strel==1,
    for ilev=2:lev,
        textse_lin_lev=textse((ilev-2)*(dim_strel^2+2)+1:(ilev-1)*(dim_strel^2+2));
        [textse_lev]=Strel_Constr_whole(textse_lin_lev(1:end-2));
        
        fprintf(['Livello ',num2str(ilev-1),', strel: '])
        fprintf('\n')
        disp(textse_lev)
        fprintf('\n')
        dim_strel1 = textse_lin_lev(end-1)*2+textse_lin_lev(end)+2;
        fprintf('Dimensione strel:  %1.0f\n',dim_strel1)
        fprintf(['Livello ',num2str(ilev-1),', strel resized: '])
        fprintf('\n')
        strel_eff=imresize(textse_lev,[dim_strel1,dim_strel1]);
        disp(strel_eff)
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
%     %% Different Operators
%     operators=individuals(ii,2:end-2);
% else
%     %%Same operators
%     operators=repmat(individuals(ii,2:end-2),1,3);
% end
% keyboard
% P = Pyr_Dec_Whole(imageHR,textse,lev,operators);
P = Pyr_Dec_Whole_grad(imageHR,textse,lev,operators,sampl_kind,int_meth);

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
        I_Fus = imageLR .* (P(:,:,:,1)./(P_LP+eps));
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
end


% keyboard

%%%%%%%%%%%%%%%%%%%%%%% Final Equalization

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


