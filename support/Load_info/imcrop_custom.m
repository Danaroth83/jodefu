function [ I_out,edge,I_out_cropped ] = imcrop_custom( I_in,vcrop,hcrop,ratio,vshift,hshift )
%IMCROP_CUSTOM cuts the section [vcrop,hcrop] of an input image I_in
%Ratio allows to scale the cropping
%hshift and vshift are shifts in the target scale

% I_out_cropped is the cropped image
% I_out is the cropped image, expanded with as much extra edge as possible
%    which will be saved to edge

if nargin<=1 || isempty(vcrop), vcrop=1:size(I_in,1); end
if nargin<=2 || isempty(hcrop), hcrop=1:size(I_in,2); end
if nargin<=3 || isempty(ratio), ratio=1; end
if nargin<=4 || isempty(vshift), vshift=0; end
if nargin<=5 || isempty(hshift), hshift=0; end

[L1,L2,~]=size(I_in);

tol=10^(-8);
if rem(ratio,1)<tol
    hcrop_FR=reshape((repmat(hcrop(:)*ratio+hshift,[1,ratio])+repmat((-ratio+1:0),[length(hcrop),1])).',1,[]);
    vcrop_FR=reshape((repmat(vcrop(:)*ratio+vshift,[1,ratio])+repmat((-ratio+1:0),[length(vcrop),1])).',1,[]);
else
    hcrop_FR=((hcrop(1)-1)*ratio+1:hcrop(end)*ratio)+hshift;
    vcrop_FR=((vcrop(1)-1)*ratio+1:vcrop(end)*ratio)+vshift;
end

I_in=padarray(I_in,[max(0,vcrop_FR(end)-L1), max(0,hcrop_FR(end)-L2)],'circular','post');
I_in=padarray(I_in,[max(0,-vcrop_FR(1)+1), max(0,-hcrop_FR(1)+1)],'circular','pre');
vcrop_FR=vcrop_FR+max(0,-vcrop_FR(1)+1); hcrop_FR=hcrop_FR+max(0,-hcrop_FR(1)+1);
[L1,L2,~]=size(I_in);

edge=min([vcrop_FR(1)-1,L1-vcrop_FR(end),hcrop_FR(1)-1,L2-hcrop_FR(end)]);

I_out=I_in([vcrop_FR(1)+(-edge:-1),vcrop_FR,vcrop_FR(end)+(1:edge)],[hcrop_FR(1)+(-edge:-1),hcrop_FR,hcrop_FR(end)+(1:edge)],:);

if nargout>2
    I_out_cropped=I_in(vcrop_FR,hcrop_FR,:);
end

end

