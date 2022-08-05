function [UIQI,UIQI_map] = img_qi_shift(I_1,I_2,Blocksize,shift)
% Custom shift mobile window version of UIQI
%
% keyboard
if nargin<=3 || isempty(shift)
    shift=Blocksize;
end
I_1 = double(I_1);
I_2 = double(I_2);

edge_extension=[floor((Blocksize-shift)/2);ceil((Blocksize-shift)/2)];
I_1=edge_extend(I_1,edge_extension);
I_2=edge_extend(I_2,edge_extension);
[L1e,L2e]=size(I_1);
stepx=floor((L1e-Blocksize)/shift+1);
stepy=floor((L2e-Blocksize)/shift+1);

UIQI_map=zeros(stepx,stepy);
for i1=1:stepx
    iA=((i1-1)*shift)+(1:Blocksize);
    for i2=1:stepy
        iB=((i2-1)*shift)+(1:Blocksize);
        UIQI_map(i1,i2)=img_qi(I_1(iA,iB),I_2(iA,iB),Blocksize);
    end;
end; 
UIQI=mean2(UIQI_map);

