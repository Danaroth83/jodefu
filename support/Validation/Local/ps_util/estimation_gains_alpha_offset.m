function [g] = estimation_gains_alpha_offset(imageLR,imageHR,GT,type_estimation)

if strcmp(type_estimation,'global')
    %%%%%%%% Global estimation
    IHc = imageHR(:);
    ILRc = reshape(imageLR,[size(imageLR,1)*size(imageLR,2) size(imageLR,3)]);
    if size(IHc,1)~=size(ILRc,1)
        keyboard
    end
    Diffc = reshape(GT - imageLR,[size(GT,1)*size(GT,2) size(GT,3)]);
%     g = [ILRc,IHc,ones(numel(IHc),1)]\Diffc;
        g = pinv([ILRc,IHc,ones(numel(IHc),1)])*Diffc;
else  %%%ToDo
    %%%%%%%% Local estimation
    block_win = 55;
    gs = zeros(size(imageLR,3) + 1,size(imageLR,3));
    cont_bl = 0;
    for ii = 1 : block_win : size(imageLR,1)
        for jj = 1 : block_win : size(imageLR,2)
                imHRbl = imageHR(ii : min(size(imageLR,1),ii + block_win - 1), jj : min(size(imageLR,2),jj + block_win - 1));
                GTbl = GT(ii : min(size(imageLR,1),ii + block_win - 1), jj : min(size(imageLR,2),jj + block_win - 1),:);
                imageLRbl = imageLR(ii : min(size(imageLR,1),ii + block_win - 1), jj : min(size(imageLR,2),jj + block_win - 1),:);
                imageHRc = reshape(imHRbl,[numel(imHRbl) 1]);
                ILRc = reshape(imageLRbl,[size(imageLRbl,1).*size(imageLRbl,2) size(imageLRbl,3)]);
                Diffc = reshape(GTbl - imageLRbl,[numel(GTbl)./size(GTbl,3) size(GTbl,3)]);     
                gh = [ILRc,imageHRc]\Diffc;
                gs = gs + gh;
                cont_bl = cont_bl + 1;
        end
    end
    g = gs/cont_bl;
end

end