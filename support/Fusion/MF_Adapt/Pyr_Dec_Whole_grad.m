function  P = Pyr_Dec_Whole_grad(Im,textse,lev,operators,sampl_kind,int_meth)


%%% textse contains the strels from level 2 to level lev
%%% operators contains the operators  from level 2 to level lev

% keyboard
dim_strel=sqrt((numel(textse)-6)/3);
dim_oper=numel(operators)/3;

P(:,:,:,1) = Im;
Sizes(1,:)=[size(Im,1), size(Im,2)];
imageI_new=P(:,:,:,1);
first=1;
oper_lev=operators;
% keyboard
for ii = 2 : lev
    
    imageI_old = imageI_new;
    clear imageI_new
    textse_lin_lev = textse((ii-2)*(dim_strel^2+2)+1:(ii-1)*(dim_strel^2+2));
    
    dim_strel1 = textse_lin_lev(end-1)*2+textse_lin_lev(end)+2;
    %     textse_lev = Strel_Constr_whole(textse_lin_lev(1:end-2));
%     oper_lev = operators((ii-2)*dim_oper+1:(ii-1)*dim_oper);
    
    %         Resizing=dim_strel1;
    %     keyboard
    textse_lev = Strel_Constr_whole(textse_lin_lev(1:end-2));
    used_strel=imresize(textse_lev,[dim_strel1,dim_strel1]);
%     keyboard
 %% operators   
    if (oper_lev(1)==0) && (oper_lev(2)==1) && (oper_lev(3)==0)
        %  Top Hat
        PC = imclose(imageI_old,used_strel);
        PO = imopen(imageI_old,used_strel);
        WTH=imageI_old-PO;
        BTH=PC-imageI_old;
        D=WTH-BTH;
        PS = imageI_old -0.5*D;
        % PS=0.5*squeeze(PC+PO); %equivalently
        
        
    elseif (oper_lev(1)==1) && (oper_lev(2)==0) && (oper_lev(3)==0)
        
        % Toggle
        PD = imdilate(imageI_old,used_strel);
        PE = imerode(imageI_old,used_strel);
        
        TCO=PD.*((PD-imageI_old)<(imageI_old-PE))+PE.*((PD-imageI_old)>(imageI_old-PE))...
            +imageI_old.*((PD-imageI_old)==(imageI_old-PE));
        DTCO=max(TCO-imageI_old,0);
        ETCO=max(imageI_old-TCO,0);
        D=DTCO-ETCO;
        % PS=imageI_old-D;
        PS=imageI_old-0.5*D;
        
    elseif (oper_lev(1)==1) && (oper_lev(2)==1) && (oper_lev(3)==0)
        
        % Bai (Toggle+TopHat)
        PD = imdilate(imageI_old,used_strel);
        PE = imerode(imageI_old,used_strel);
        
        TCO=PD.*((PD-imageI_old)<(imageI_old-PE))+PE.*((PD-imageI_old)>(imageI_old-PE))...
            +imageI_old.*((PD-imageI_old)==(imageI_old-PE));
        DTCO=max(TCO-imageI_old,0);
        ETCO=max(imageI_old-TCO,0);
        PC = imclose(imageI_old,used_strel);
        PO = imopen(imageI_old,used_strel);
        
        WR=imageI_old-PO;
        BR=PC-imageI_old;
        D=WR-BR+DTCO-ETCO;
        PS=imageI_old-0.5*D; %modified for pansharpening 
%         PS=imageI_old-D; %original Bai
    elseif (oper_lev(1)==0) && (oper_lev(2)==0) && (oper_lev(3)==0)
        %  Gradient
        PD = imdilate(imageI_old,used_strel);
        PE= imerode(imageI_old,used_strel);
        rho_minus=imageI_old-PE;
        rho_plus=PD-imageI_old;
        D=rho_minus-rho_plus;
        PS = imageI_old -0.5*D;
        % PS = 0.5*squeeze(PD+PE); %equivalently
    elseif (oper_lev(1)==0) && (oper_lev(2)==0) && (oper_lev(3)==1)
        % LOCO
        POC = imopen(imclose(imageI_old,used_strel),used_strel);
        PCO= imclose(imopen(imageI_old,used_strel),used_strel);
        PS = 0.5*squeeze(POC+PCO);
%         keyboard
        D= imageI_old -PS;
    elseif (oper_lev(1)==0) && (oper_lev(2)==1) && (oper_lev(3)==1)
        % Median
%         keyboard
        clear PS
        num_elnotzero=sum(textse_lev(:));
        if num_elnotzero>0
        for iband=1:size(imageI_old,3),
        PS(:,:,iband) = ordfilt2(imageI_old(:,:,iband),ceil(num_elnotzero/2),textse_lev,'symmetric');
        end
        
%         keyboard
        else
           PS= imageI_old ;
        end
        D= imageI_old -PS;
    
    elseif (oper_lev(1)==1) && (oper_lev(2)==0) && (oper_lev(3)==1)
        PORCR = imopenbyrec(imclosebyrec(imageI_old,used_strel),used_strel);
        PCROR= imclosebyrec(imopenbyrec(imageI_old,used_strel),used_strel);
        PS = 0.5*squeeze(PORCR+PCROR);
%         keyboard
        D= imageI_old -PS;
     end   
%      keyboard
    clear PC
    clear PO
    if strcmp(sampl_kind,'mean')
        dim_strel2=2;%dim_strel1;
        for il=1:size(imageI_old,3)
            for i1=1:floor(size(imageI_old,1)/dim_strel2),
                for i2=1:floor(size(imageI_old,2)/dim_strel2),
                    imageI_new(i1,i2,il)=mean2(PS((i1-1)*dim_strel2+1:(i1)*dim_strel2,...
                        (i2-1)*dim_strel2+1:(i2)*dim_strel2,il));
                end
            end
        end
        % imageI_new=imresize(imageI_old,1/dim_strel1);
    elseif strcmp(sampl_kind,'downs')
        if first
            for il=1:size(imageI_old,3)
                imageI_new(:,:,il)=PS(2:2:end,2:2:end,il);
                %                imageI_new(:,:,il)=PS(3:2:end,3:2:end,il);
            end
            first=0;
        else
            for il=1:size(imageI_old,3)
                imageI_new(:,:,il)=PS(1:2:end,1:2:end,il);
            end
        end
    end
    Sizes(ii,:)=[size(imageI_new,1) size(imageI_new,1)];
    imageI_resized_old=imageI_new;
%     keyboard
    if strcmp(int_meth,'interp23tap')
        %Only for specific image size and ratio
        for il=1:size(Im,3)
            imageI_resized_new(:,:,il)  = interp23tap(imageI_resized_old(:,:,il),2^(ii-1));
        end
    imageI_resized_old=imageI_resized_new;
        clear imageI_resized_new
    else
        for ir=ii:-1:2,
        for il=1:size(Im,3)
            imageI_resized_new(:,:,il)  = imresize(imageI_resized_old(:,:,il),[Sizes(ir-1,1) Sizes(ir-1,2)],int_meth);
        end
        imageI_resized_old=imageI_resized_new;
        clear imageI_resized_new
        end
    end
    
    
    %%%%%%%%%%%%%%%%% Without sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     imageI_resized_old=PS;
    %     imageI_new=PS;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sum(isfinite(imageI_resized_old(:)))~=numel(imageI_resized_old)
        P(:,:,:,1:lev) =repmat(P(:,:,:,1),1,1,1,lev);
        break
    else
        P(:,:,:,ii) = imageI_resized_old;
    end
    %     P(:,:,:,ii) = imageI_resized_old;
    clear imageI_resized_old
end

