function [D_BDSD] = sparse_coding_BDSD (D_PAN, D_MStilde, D_PAN_d, D_MS_d,D_MS_orig, ratio, TS, ol,num_atoms,Height,Width,Bands,opt_met, det_met)

%%% Coefficient calculation as gamma=pinv(DL)*(MS_LR-MS_LR^LP)
%%% Image Reconstruction  as  MShat-MStilde = DH*gamma;

M=size(D_MStilde,2)/(Bands);  %Number of atoms
% Same method used in Dictionary Construction
    DH = [D_MStilde, D_PAN, ones(size(D_PAN,1),1)];
    DL = [D_MS_d, D_PAN_d, ones(size(D_PAN_d,1),1)];

    
% keyboard

D_BDSD = zeros(Height,Width,Bands);
Count_Tiles = zeros(Height,Width);
for kkk = 1 : M,
    D_MS_orig_kkk=D_MS_orig(:,(kkk-1)*Bands+1:kkk*Bands);
    D_MS_d_kkk=D_MS_d(:,(kkk-1)*Bands+1:kkk*Bands);
    Diff_MS = D_MS_orig_kkk-D_MS_d_kkk; 
   switch opt_met
       case 'all'
    %%% All Coefficients
    gamma = DL\Diff_MS;
       case 'OMP'
           %%% Or
           delta=0;
           gamma=zeros(size(DH,2),Bands);
           for iband=1:Bands,
               [gamma_ord, indx] = OMP_BDSD(DL, Diff_MS(:,iband), delta, num_atoms);
               gamma(indx,iband)=gamma_ord;
           end
           %%% Negative gamma check (TBD)
           %     if ng_check
           %         if any(gamma(end-1,:)<0) %% Pan coefficients
           %             gamma=gamma_global;
           %         end
           %     end          
   end
   D_BDSD_col(:,(kkk-1)*Bands+1:kkk*Bands)=DH*gamma;
end

D_BDSD = detile_BDSD (D_BDSD_col,det_met, Height, Width, Bands, ratio, TS, ol);








