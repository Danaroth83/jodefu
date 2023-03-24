function I_out=compress_image(I_in,comp_factor,comp_method,L)

if comp_factor<2^L/2^8 && strcmpi(comp_method,'JP2K')
    comp_method='radiometric';
    warning('No JPEG2k compression was performed');
end

I_out=I_in;

if comp_factor~=1
    if strcmpi(comp_method,'JP2K')
        I_out=uint8((I_in)*2^8/2^L);
        imwrite(I_out,'image_temp.j2k','CompressionRatio',comp_factor*2^8/2^L);
        I_out=imread('image_temp.j2k');
        I_out=double(I_out)*2^L/2^8;
    elseif strcmpi(comp_method,'radiometric')
        I_out=round(I_in/comp_factor)*comp_factor;
    elseif strcmpi(comp_method,'rad_extend')
        MI=max(I_in(:));
        mI=min(I_in(:));
        I_out=(I_in-mI)/(MI-mI)*2^L;
        I_out=round(I_out/comp_factor)*comp_factor;
        I_out=I_out/2^L*(MI-mI)+mI;
    end
end