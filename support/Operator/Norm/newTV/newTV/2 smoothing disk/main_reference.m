function main
	
	lambda = 6;
	radius = 32;
	% cell-average sampling approximated by point sampling on a 16x larger grid followed by downsampling
	S = 99;
	y=zeros(S*16,S*16);
	[n1,n2]=meshgrid(1:S*16,1:S*16);
	y=((n1-S*8-1/2).^2+(n2-S*8-1/2).^2<=(radius*16)^2);
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	
	tmp=2*pi*lambda*radius/(S^2-pi*radius^2);
	y=y*(1-2*lambda/radius-tmp)+tmp;
	
	figure(1);
	imshow(y);
	colormap gray
	imwrite(y,'x_reference.tif');
	
end
