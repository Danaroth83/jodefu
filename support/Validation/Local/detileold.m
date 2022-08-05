function [D_BDSD,Count_Tiles] = detile (D_BDSD,Count_Tiles,D_BDSD_col(:,(kkk-1)*Bands+1:kkk*Bands));
 
 nr = ceil ((Height/ratio - ol) / (TS - ol));
nc = ceil ((Width/ratio - ol) / (TS - ol));

for irow=1:nr
    for icol=1:nc
        =
        icount = icount + 1;
        shiftr = 0; shiftc = 0;
        if irow == nr && mod(size(I_MS_d,1)-ol, TS-ol) ~= 0
            shiftr = TS-ol - mod (size(I_MS_d,1)-ol, TS-ol);
        end
        if icol == nc && mod(size(I_MS_d,2)-ol, TS-ol) ~= 0
            shiftc = TS-ol - mod (size(I_MS_d,2)-ol, TS-ol);
        end
        blockr = ((irow-1)*(TS-ol)*ratio+1 - shiftr*ratio) : ((irow*TS-(irow-1)*ol)*ratio- shiftr*ratio);
        blockc = ((icol-1)*(TS-ol)*ratio+1 - shiftc*ratio) : ((icol*TS-(icol-1)*ol)*ratio- shiftc*ratio);
        
%         blockrl = ((irow-1)*(TS-ol)+1 - shiftr) : (irow*TS-(irow-1)*ol - shiftr);
%         blockcl = ((icol-1)*(TS-ol)+1 - shiftc) : (icol*TS-(icol-1)*ol - shiftc);
D_BDSD 