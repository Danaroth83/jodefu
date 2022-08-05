function main 
	lambda=6/32;
	S=100/64; % taille du carré
	delta=16*((2*lambda-(S^2-1))^2-lambda^2*(pi-4)*(S^2-1));
	r=(-8*lambda+4*(S^2-1)-sqrt(delta))/(pi-4)/2/lambda;
	v=1-lambda/r; %background value = 0.253154205548092
	
	imv=zeros(1600,1600);
	for xx=1:1600  
	for yy=1:1600
		imv(yy,xx)=denoisedsquare((xx-800.5)/32/16,(yy-800.5)/32/16,lambda);
	end
	end
	imv2=max(imv,v);
	imv2=(imv2(1:2:end,1:2:end)+imv2(2:2:end,1:2:end)+imv2(1:2:end,2:2:end)+imv2(2:2:end,2:2:end))/4;
	imv2=(imv2(1:2:end,1:2:end)+imv2(2:2:end,1:2:end)+imv2(1:2:end,2:2:end)+imv2(2:2:end,2:2:end))/4;
	imv2=(imv2(1:2:end,1:2:end)+imv2(2:2:end,1:2:end)+imv2(1:2:end,2:2:end)+imv2(2:2:end,2:2:end))/4;
	imv2=(imv2(1:2:end,1:2:end)+imv2(2:2:end,1:2:end)+imv2(1:2:end,2:2:end)+imv2(2:2:end,2:2:end))/4;
	%sum(sum(imv2)) gives 4096.01 -> ideally, should be equal to sum(sum(y)) = 4096
	imwrite(imv2,'x_reference.tif');
	figure(1); 
	imshow(imv2);
	colormap gray
end

function v = denoisedsquare(x,y,lambda)
	x = abs(x);
	y = abs(y);
	if (x>=1)||(y>=1), v=0; return; end
	if lambda>=1/(1+sqrt(pi)/2), v=0; return; end
	% we looke for r such that |x| = 1-r + r.cos(theta) et y = 1-r + r.sin(theta)
	% r is solution to x^2+y^2+2-2*x-2*y+r*(2*x+2*y-4)+r^2 = 0
	r = (4-2*x-2*y+sqrt(8*(1-x)*(1-y)))/2; 
	if r<=lambda, v=0; return; 
	elseif r>=1/(1+sqrt(pi)/2), v=1-lambda*(1+sqrt(pi)/2); return; 
	else, v = 1-lambda/r; end
end