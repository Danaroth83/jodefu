function D_BDSD = detile_BDSD (D_BDSD_col, det_met, Height, Width, Bands, ratio, TS, ol)

% keyboard
% spacchetta bande
% if size (I_Fus_CS_column,3) < C_MS
%     I_Fus_B = zeros (TS^2*Resize_fact^2, size(I_Fus_CS_column,2), C_MS);
M=size(D_BDSD_col,2)/(Bands);  %Number of atoms
    
%     I_Fus_CS_column = I_Fus_B;
% end


D_BDSD = zeros ([Height Width Bands]);
countpx = zeros ([Height Width Bands]);
nr = ceil ((Height/ratio - ol) / (TS - ol));
nc = ceil ((Width/ratio - ol) / (TS - ol));
for iband=1:Bands
        D_BDSD_col_B = D_BDSD_col(:,(0:Bands:(M-1)*Bands)+iband);
    icount=0;
    for irow=1:nr
        for icol=1:nc
            shiftr = 0; shiftc = 0;
            if irow == nr && mod(Height/ratio-ol, TS-ol) ~= 0
                shiftr = TS-ol - mod (Height/ratio-ol, TS-ol);
            end
            if icol == nc && mod(Width/ratio-ol, TS-ol) ~= 0
                shiftc = TS-ol - mod (Width/ratio-ol, TS-ol);
            end
            blockr = ((irow-1)*(TS-ol)*ratio+1 - shiftr*ratio) : ((irow*TS-(irow-1)*ol)*ratio - shiftr*ratio);
            blockc = ((icol-1)*(TS-ol)*ratio+1 - shiftc*ratio) : ((icol*TS-(icol-1)*ol)*ratio - shiftc*ratio);
            icount=icount+1;
            if strcmpi (det_met, 'average')
                D_BDSD(blockr,blockc,iband) = D_BDSD(blockr,blockc,iband) + reshape(D_BDSD_col_B(1:length(blockr)*length(blockc),icount),length(blockr),length(blockc));
                countpx(blockr,blockc,iband) = countpx(blockr,blockc,iband) +1;
            else
                D_BDSD(blockr,blockc,iband) = reshape(D_BDSD_col_B(1:length(blockr)*length(blockc),icount),length(blockr),length(blockc));
            end
        end
    end
end
if strcmpi (det_met, 'average')
    D_BDSD = D_BDSD ./ countpx;
end

end