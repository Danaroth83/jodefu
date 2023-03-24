function a=exchangepatches(I_out,coord,patches,Dim,S)

L1e=size(I_out,1)+floor((Dim-S)/2)+ceil((Dim-S)/2);
L2e=size(I_out,2)+floor((Dim-S)/2)+ceil((Dim-S)/2);
% stepx=floor((L1e-Dim)/S+1);
stepy=floor((L2e-Dim)/S+1);

%coord=coord-floor(S/2)+floor((Dim-S)/2);  %check
coord=round(coord/S);
coord=(coord(1)-1)*stepy+coord(2);

refer=patches(:,coord);
refer2=zeros(length(refer),2);
refer2(:,1)=floor(refer/stepy)+1;
refer2(:,2)=refer-(refer2(:,1)-1)*stepy;

% a=refer2*S+floor(S/2)-floor((Dim-S)/2);  % check
a=refer2;

end


