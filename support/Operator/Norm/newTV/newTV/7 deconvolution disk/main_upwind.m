function main 	
	
	Nbiter= 4000;
	lambda = 0.1; 
	mu = 1.2; 
	tau = 0.99/(16+mu);
	gamma = 0.99;
	blurfilter = exp(-(-3:0.2:3).^2);
	blurfilter = blurfilter/sum(blurfilter);
	noiselevel = 0.05; % std. dev. of the noise added to the blurred image
	
	% cell-average sampling approximated by point sampling on a 16x larger grid followed by downsampling
	S=99;
	y = zeros(S*16,S*16);
	[n1,n2]=meshgrid(1:S*16,1:S*16);
	y=((n1-S*8-1/2).^2+(n2-S*8-1/2).^2<=(32*16)^2);
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y = y*0.8+0.1;
	xref = y;
	figure(1);
	imshow(y);
	colormap gray
	imwrite(y,'reference.tif');
	
	opblur = @(x) imfilter(imfilter(x,blurfilter,'symmetric'),blurfilter','symmetric');
	% blur with symmetric boundary conditions. opblur is self-adjoint	
	rng(100);
	noise = randn(size(y));
	noise = noise/norm(noise,'fro')*sqrt(numel(noise))*noiselevel; % noise of std. dev. 0.1
	y = opblur(y)+noise;
	figure(2);
	imshow(y);
	colormap gray
	imwrite(y,'y.tif');	
	
	opD = @(x) cat(3,[-diff(x,1,1);zeros(1,size(x,2))],[-diff(x,1,2) zeros(size(x,1),1)],...
		[zeros(1,size(x,2));diff(x,1,1)],[zeros(size(x,1),1) diff(x,1,2)]);
	opDadj = @(u) [u(1,:,1);diff(u(:,:,1),1,1)]+[u(:,1,2) diff(u(:,:,2),1,2)]+...
		[-diff(u(:,:,3),1,1);u(end,:,3)]+[-diff(u(:,:,4),1,2) u(:,end,4)];	
	
	x = y;
	u = zeros([size(x) 4]);
	v = opD(x);
	fprintf('0 %f\n',sum(sum((opblur(x)-y).^2))/2+lambda*sum(sum(sqrt(sum(max(v,0).^2,3)))));
	for iter = 1:Nbiter
		x = x-tau*opDadj(opD(x)-v+mu*u)-tau*mu*opblur(opblur(x)-y);
		v = prox_mu_gamma_g(v+gamma*(opD(x)-v+mu*u),mu*gamma*lambda);
		u = u+(opD(x)-v)/mu;
		if mod(iter,40)==0
			%we display the primal cost value
			fprintf('%d %f\n',iter,sum(sum((opblur(x)-y).^2))/2+...
				lambda*sum(sum(sqrt(sum(max(opD(x),0).^2,3)))));
			figure(3);
			imshow(x);
			colormap flag
			drawnow
		end
	end
	figure(3);
	imshow(x);
	colormap gray
	imwrite(x,'x_upwind.tif');
	fprintf('restoration error: %f\n',norm(x-xref,'fro')); %=2.77
end

function unew = prox_mu_gamma_g(u,mugammalambda) 
	tmp = max(u,0);
	unew = u-bsxfun(@rdivide, tmp, max(sqrt(sum(tmp.^2,3))/mugammalambda,1));
end
