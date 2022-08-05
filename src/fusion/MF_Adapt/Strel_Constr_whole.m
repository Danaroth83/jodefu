function text_se_lev=Strel_Constr_whole(text_se)

%%% Structruring Element construction
% keyboard
dim_orig=sqrt(numel(text_se));
if dim_orig==3,
    text_se_lev=zeros(3,3);
    text_se_lev(1:3,1:3)=reshape(text_se,3,3);
%     text_se_lev(1:3,4:5)=text_se_lev(1:3,2:-1:1);
%     text_se_lev(4:5,:)=text_se_lev(2:-1:1,:);
    text_se_lev=text_se_lev>0;
elseif dim_orig==2,
    text_se_lev=zeros(2,2);
    text_se_lev(1:2,1:2)=reshape(text_se,2,2);
%     text_se_lev(1:3,3:4)=text_se_lev(1:3,2:-1:1);
%     text_se_lev(3:4,:)=text_se_lev(2:-1:1,:);
    text_se_lev=text_se_lev>0;
end



