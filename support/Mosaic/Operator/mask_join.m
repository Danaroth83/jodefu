function I_mosaic=mask_join(I_in,mask_in)

    if isfield(I_in,'HR'), I_HR=I_in.HR; else, I_HR=[]; end
    if isfield(I_in,'data'), I_in=I_in.data; elseif isfield(I_in,'LR'), I_in=I_in.LR; end
    if isfield(mask_in,'data'), mask=mask_in.data; else, mask=mask_in; end
    if isfield(mask_in,'shift'), shift=mask_in.shift; else, shift=zeros(size(mask,3),2); end
    if isfield(mask_in,'padding'), padding=mask_in.padding; else, padding=[0,0,0,0]; end

    ratio=[size(I_HR,1)/size(I_in,1),size(I_HR,2)/size(I_in,2)];
    if all(ratio>=1), I_in=imresize(I_in,[size(I_in,1)*ratio(1),size(I_in,2)*ratio(2)]); end
    I_both=cat(3,I_HR,I_in);
    
    if isempty(mask)
        I_mosaic=I_both; 
    else
        I_mosaic=opdir_mosaicshift(I_both,mask,shift,padding);
    end
end