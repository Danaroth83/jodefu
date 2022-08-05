function [I_out,Mat_mseidx]=NLM_personal(I_in,Dim,S,neighb,sigma)

sigmasq=sigma^2;
[L1,L2,~]=size(I_in);
I_in=edge_extend(I_in,[floor((Dim-S)/2);ceil((Dim-S)/2)]);
[L1e,L2e,Nb]=size(I_in);
stepx=floor((L1e-Dim)/S+1);
stepy=floor((L2e-Dim)/S+1);


Dimsq=Dim^2;
kernel=reshape(repmat(make_kernel_personal(floor(Dim/2)),[1,1,Nb]),[Dimsq*Nb,1]);
kernel=kernel/sum(kernel);
Mat=zeros(Dimsq*Nb,stepx*stepy);
for x=1:stepx
    iA=((x-1)*S)+(1:Dim);
    for y=1:stepy
        iB=((y-1)*S)+(1:Dim);
        Mat(:,(x-1)*stepy+y)=reshape(I_in(iA,iB,:),[Dimsq*Nb,1]);
    end
end

Mat_mseidx=zeros(neighb,stepx*stepy);
Mat2=zeros(Dimsq*Nb,stepx*stepy);
for x=1:stepx
    iA=((x-1)*S)+(1:Dim);
    for y=1:stepy
        iB=((y-1)*S)+(1:Dim);
        currentblock=reshape(I_in(iA,iB,:),[Dimsq*Nb,1]);
        mse=sum(repmat(kernel,[1,stepx*stepy]).*((repmat(currentblock,[1,stepx*stepy])-Mat).^2),1);
        [mse,idx]=sort(mse,'ascend');
        Mat_mseidx(:,(x-1)*stepy+y)=idx(1:neighb);
        weight=exp(-mse(1:neighb)/sigmasq);
        Mat2(:,(x-1)*stepy+y)=sum(repmat(weight,[Dimsq*Nb,1]).*Mat(:,idx(1:neighb)),2)/sum(weight);
    end
end

% I_out=zeros(size(I_in));
% I_outcount=zeros(L1e,L2e);
% for x=1:stepx
%     iA=((x-1)*S)+(1:Dim);
%     for y=1:stepy
%         iB=((y-1)*S)+(1:Dim);
%         I_out(iA,iB,:)=I_out(iA,iB,:)+reshape(Mat2(:,(x-1)*stepy+y),[Dim,Dim,Nb]);
%         I_outcount(iA,iB)=I_outcount(iA,iB)+ones(Dim,Dim);
%     end
% end
% I_out=I_out./repmat(I_outcount,[1,1,Nb]);
% I_out=I_out(floor((Dim-S)/2)+1:end-(ceil((Dim-S)/2)),floor((Dim-S)/2)+1:end-(ceil((Dim-S)/2)),:);

I_out=zeros(L1,L2,Nb);
for x=1:stepx
    for y=1:stepy
        temp=reshape(Mat2(:,(x-1)*stepy+y),[Dim,Dim,Nb]);
        I_out(x,y,:)=temp(floor((Dim-S)/2)+1,floor((Dim-S)/2)+1,:);
    end
end
% I_out=I_out(floor((Dim-S)/2)+1:end-(ceil((Dim-S)/2)),floor((Dim-S)/2)+1:end-(ceil((Dim-S)/2)),:);


end

