function main 	
	
	Nbiter= 3000;
	lambda = 6; 
	tau = 0.99/8;
	sigma = 0.99/3;
	mu = 0.1; 
	
	% cell-average sampling approximated by point sampling on a 16x larger grid followed by downsampling
	S=99;
	y=zeros(S*16,S*16);
	[n1,n2]=meshgrid(1:S*16,1:S*16);
	y=((n1-S*8-1/2).^2+(n2-S*8-1/2).^2<=(32*16)^2);
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	figure(1);
	imshow(y);
	colormap gray
	imwrite(y,'y.tif');
	
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_mu_tau_f = @(x) (x+mu*tau*y)/(1+mu*tau);
	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,3))/(lambda*mu*sigma),1));
	
	x = y;
	u = zeros([size(x) 2]);
	v = zeros([size(x) 2 3]);
	tmp = opD(x);
	v(:,:,1,1) = tmp(:,:,1);
	v(:,:,2,2) = tmp(:,:,2);
	fprintf('0 %f\n',lambda*sum(sum(sum(sqrt(sum(v.^2,3))))));
	for iter = 1:Nbiter
		x = prox_mu_tau_f(x+tau*opDadj(-opD(x)+opLadj(v)-mu*u));
		v = prox_mu_sigma_g(v-sigma*opL(-opD(x)+opLadj(v)-mu*u));
		u = u-(-opD(x)+opLadj(v))/mu;
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			fprintf('%d %f %f\n',iter,sum(sum((x-y).^2))/2+lambda*sum(sum(sum(sqrt(sum(v.^2,3))))),...
				-sum(sum((y-opDadj(u)).^2-y.^2))/2);
			figure(2);
			imshow(x);
			colormap flag
			drawnow
		end
	end
	figure(2);
	imshow(x);
	colormap gray
	imwrite(x,'x_proposed.tif');
	
	fig3=figure(3);
	imshow(x);
	colormap gray
	hold on
	[n1,n2]=meshgrid(1:S,1:S);
	quiver(n1,n2+0.5,v(:,:,2,1),v(:,:,1,1),0,'r');
	quiver(n1+0.5,n2,v(:,:,2,2),v(:,:,1,2),0,'b');
	quiver(n1,n2,v(:,:,2,3),v(:,:,1,3),0,'g');
	truesize(fig3,[S*4,S*4]);
end


function t = opL(u)
	[height,width,d]=size(u);
	t=zeros(height,width,2,3);
	t(:,:,1,1)=u(:,:,1); 
	t(1:end-1,2:end,2,1)=(u(2:end,1:end-1,2)+u(1:end-1,1:end-1,2)+u(2:end,2:end,2)+u(1:end-1,2:end,2))/4; 			t(1:end-1,1,2,1)=(u(1:end-1,1,2)+u(2:end,1,2))/4; 
	t(:,:,2,2)=u(:,:,2);
	t(2:end,1:end-1,1,2)=(u(2:end,1:end-1,1)+u(1:end-1,1:end-1,1)+u(2:end,2:end,1)+u(1:end-1,2:end,1))/4; 
	t(1,1:end-1,1,2)=(u(1,1:end-1,1)+u(1,2:end,1))/4; 	
	t(2:end,:,1,3) = (u(2:end,:,1)+u(1:end-1,:,1))/2; 
	t(1,:,1,3) = u(1,:,1)/2;
	t(:,2:end,2,3) = (u(:,2:end,2)+u(:,1:end-1,2))/2; 
	t(:,1,2,3) = u(:,1,2)/2;
end

function u = opLadj(t)
	[height,width,d,c]=size(t);
	u=zeros(height,width,2);
	u(1:end-1,2:end,1)=t(1:end-1,2:end,1,1)+(t(1:end-1,2:end,1,2)+t(1:end-1,1:end-1,1,2)+...
	t(2:end,2:end,1,2)+t(2:end,1:end-1,1,2))/4+(t(1:end-1,2:end,1,3)+t(2:end,2:end,1,3))/2;
	u(1:end-1,1,1)=t(1:end-1,1,1,1)+(t(1:end-1,1,1,2)+t(2:end,1,1,2))/4+...
	(t(1:end-1,1,1,3)+t(2:end,1,1,3))/2;
	u(2:end,1:end-1,2)=t(2:end,1:end-1,2,2)+(t(2:end,1:end-1,2,1)+t(1:end-1,1:end-1,2,1)+...
	t(2:end,2:end,2,1)+t(1:end-1,2:end,2,1))/4+(t(2:end,1:end-1,2,3)+t(2:end,2:end,2,3))/2;
	u(1,1:end-1,2)=t(1,1:end-1,2,2)+(t(1,1:end-1,2,1)+t(1,2:end,2,1))/4+...
	(t(1,1:end-1,2,3)+t(1,2:end,2,3))/2;
end
	
	