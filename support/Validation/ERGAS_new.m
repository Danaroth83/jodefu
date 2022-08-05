function [ERGAS_index,ERGAS_map,ERGAS_band]=ERGAS_new(I_GT,I_F,ratio)
    [L1,L2,Nb]=size(I_GT);
    ERGAS_band=zeros([L1,L2,Nb]);
    for ii=1:Nb
        I_GT_band=I_GT(:,:,ii);
        I_F_band=I_F(:,:,ii);
        ERGAS_band(:,:,ii) = (I_GT_band-I_F_band).^2/(mean2(I_GT_band))^2;
    end
    ERGAS_index = (100/ratio) * sqrt((1/Nb) * sum(ERGAS_band));
    ERGAS_map=mean(ERGAS_band,3);
    if nargout==3
        ERGAS_band=squeeze(mean(mean(ERGAS_map,1),2));
    end
end