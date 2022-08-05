function   [D_PAN, D_MStilde, D_PAN_d, D_MS_d, D_MS_orig] = Dict_Learn(I_PAN, I_MS, I_PAN_d, I_MS_d, I_MS_orig,ratio, TS, ol,SVDflag);

%  keyboard
%%%Input variables
%
% TS=tile size
% ol = overlap

%%%Output variables
% D_PAN=High resolution Dictionary from  PAN
%  D_MStilde=High resolution Dictionary from  MStilde
% D_PAN_d=High resolution Dictionary from  PAN_d
%  D_MS_d=High resolution Dictionary from  MS_d
%  D_MS_orig=High resolution Dictionary from  MS_orig



%R = TS + (nr-1)*(TS-ol)
nr = ceil ((size(I_PAN,1)/ratio- ol) / (TS - ol));
nc = ceil ((size(I_PAN,2)/ratio- ol) / (TS - ol));
nBands = size (I_MS_d,3);
D_PAN = zeros (TS^2*ratio^2,nr*nc);
D_MStilde = zeros (TS^2*ratio^2,nBands*nr*nc);
D_PAN_d = zeros (TS^2, nr*nc);
D_MS_d = zeros (TS^2,nBands*nr*nc);
D_MS_orig = zeros (TS^2,nBands*nr*nc);


icount = 0;
for irow=1:nr
    for icol=1:nc
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
        
        blockrl = ((irow-1)*(TS-ol)+1 - shiftr) : (irow*TS-(irow-1)*ol - shiftr);
        blockcl = ((icol-1)*(TS-ol)+1 - shiftc) : (icol*TS-(icol-1)*ol - shiftc);
        
        
        colmn_PAN = I_PAN(blockr,blockc);
        colmn_PAN_LR = I_PAN_d(blockrl,blockcl);
        D_PAN(:,icount)=colmn_PAN(:);
        D_PAN_d(:,icount)=colmn_PAN_LR(:);
        for iband = 1:nBands
            colmn_MStilde = I_MS(blockr,blockc,iband);
            colmn_MS_d = I_MS_d(blockrl,blockcl,iband);
            colmn_MS_orig = I_MS_orig(blockrl,blockcl,iband);
            D_MStilde(:,(icount-1)*nBands+iband)=colmn_MStilde(:);
            D_MS_d(:,(icount-1)*nBands+iband)=colmn_MS_d(:);
            D_MS_orig(:,(icount-1)*nBands+iband)=colmn_MS_orig(:);
        end
    end
end
% end
% S = sum(Dh.*Dh).^0.5;
% Dh=Dh./repmat(S,size(Dh,1),1);
% S = sum(Dl.*Dl).^0.5;
% Dl=Dl./repmat(S,size(Dl,1),1);
% S = sum(ytilde_k.*ytilde_k).^0.5;
% ytilde_k=ytilde_k./repmat(S,size(ytilde_k,1),1);
% keyboard
% Z = [D_PAN; D_MS];
% if SVDflag
%     cd KSVD_Matlab_ToolBox
%     par.K = size_Dict;
%     par.numIteration = num_Iter;
%     par.errorFlag = 0;
%     par.L = tau;
%     par.preserveDCAtom = 0;
%     par.InitializationMethod = 'DataElements';
%     par.displayProgress = 1;
%     %     for iband=1:size(Z,3)
%     %         D(:,:,iband) = KSVD (Z(:,:,iband), par);%% Normalized Dictionary Construction
%     %     end
%     D = KSVD_Mod (Z, par); %%Non-Normalized Dictionary Construction
%     cd ..
%     D_PAN = D(1:TS^2*ratio^2*nBands, :,:);
%     D_MS = D(TS^2*ratio^2*nBands+1:end, :,:);
% else
%     D=Z;
% end

end






