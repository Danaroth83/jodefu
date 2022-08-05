% GENERATE M FUNCTION
% 
% Description:
% Filtering with the support M function for the Intensity difference 
% demosaicking method 

function I_out=filterM(I_in,mask,period)
    M1=mask(1:period(1),1:period(2),:);
    sizeM=period+1-mod(period,2);
    M=zeros(sizeM);
    I_out=zeros(size(I_in));
    for ii=1:period(1)
        for jj=1:period(2)
            idx=M1(ii,jj,:)==1;
            M=imfilter(circshift(M1(:,:,idx),[1-ii,1-jj]),ones(sizeM),'circular','full');
            M=M(1:sizeM(1),1:sizeM(2));
            M=1./M;
            M=M/sum(M(:));
            I_temp=imfilter(I_in,M,'circular');
            I_out(ii:period(1):end,jj:period(2):end)=I_temp(ii:period(1):end,jj:period(2):end);
        end
    end
end