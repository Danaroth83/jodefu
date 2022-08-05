function ERGAS = ERGAS(I,I_Fus,Resize_fact)

I = double(I);
I_Fus = double(I_Fus);

Err=I-I_Fus;
ERGAS=0;
for iLR=1:size(Err,3),
    ERGAS=ERGAS+mean2(Err(:,:,iLR).^2)/(mean2((I(:,:,iLR))))^2;   
end

ERGAS = (100/Resize_fact) * sqrt((1/size(Err,3)) * ERGAS);

end