%% MASK VISUALIZATION
% Input: Binary mask in the mask.image; amount of HR bands in mask.band_HR
% Assigns White to pixels given to the HR image (Shades of gray if more than one)
% Assigns Colors to pixels given to the LR image (Respectively RGB and Yellow)
% Assigns Black to unassigned pixels

function vis_mask= visualize_mask(mask,flag_vis)

if nargin<=1, flag_vis=1; end

for kk=1:size(mask.data,3)
    mask.data(:,:,kk)=circshift(mask.data(:,:,kk),mask.shift(kk,:));
end

color_code_PAN=[1,1,1].*(1:-1/mask.band_HR:1/mask.band_HR).';
color_code_MS=[1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
color_code=cat(1,color_code_PAN,color_code_MS);

pink=[255,192,203]/255;
difference=size(mask.data,3)-size(color_code,1);
if difference>0
    color_code=cat(1,color_code,pink.*(1:-1/(difference):1/(difference)).');
else
    color_code=color_code(1:size(mask.data,3),:);
end

vis_mask=reshape(sum(shiftdim(color_code,-2).*mask.data,3),[size(mask.data,1),size(mask.data,2),3]);

% vis_mask=reshape(vis_mask,[size(mask.image,1)*size(mask.image,2),3]);
% idx=find(sum(mask.image,3)==0);
% vis_mask(idx,:)=[1,1,1].*ones(length(idx),1);
% vis_mask=reshape(vis_mask,[size(mask.image,1),size(mask.image,2),3]);

if flag_vis==1
    figure; imshow(imresize(vis_mask,max(round(256/max(size(mask.data))),1),'nearest'),[]);
end

