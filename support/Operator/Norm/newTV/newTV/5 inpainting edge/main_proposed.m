function main
	
	Nbiter= 3000;
	tau = 0.9/8;
	sigma = 0.9/3;
	mu = 1;

	S = 28;
	S2 = 5;
	[n1,n2]=meshgrid(1:S,1:S);
	y  = double(n1+n2>=S+2);
	idx=find(n1+n2==S+1);
	y(idx)=0.5;
	
	idx=find((n1>S2)&(n1<=S-S2)&(n2>S2)&(n2<=S-S2));
	idx2=find((n1<=S2)|(n1>S-S2)|(n2<=S2)|(n2>S-S2));
	y(idx)=0;
	yidx=y(idx2);
	ycolor=zeros(S,S,3);
	ycolor(:,:,1)=y;
	ycolor(:,:,2)=y;
	ycolor(:,:,3)=y;
	ycolor(S2+1:S-S2,S2+1:S-S2,3)=1;
	
	figure(1);
	imshow(y);
	figure(2);
	imshow(ycolor);
	imwrite(y,'y.tif');
	imwrite(ycolor,'ycolor.tif');
			
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,3))/(mu*sigma),1));

	x = y;
	u = zeros([size(x) 2]);
	v = zeros([size(x) 2 3]);
	tmp = opD(x);
	v(:,:,1,1) = tmp(:,:,1);
	v(:,:,2,2) = tmp(:,:,2);
	opDy = opD(y);
	opDy = opDy(:);
	fprintf('0 %f\n',sum(abs(opDy))); % == sum(sum(sum(sqrt(sum(v.^2,3)))))
	for iter = 1:Nbiter
		x = prox_mu_tau_f(x+tau*opDadj(-opD(x)+opLadj(v)-mu*u),idx2,yidx);
		v = prox_mu_sigma_g(v-sigma*opL(-opD(x)+opLadj(v)-mu*u));
		u = u-(-opD(x)+opLadj(v))/mu;
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			fprintf('%d %f %f\n',iter,sum(sum(sum(sqrt(sum(v.^2,3))))),...
				sum(u(:).*opDy)); % == sum(sum(opDadj(u).*y))
			figure(3);
			imshow(x);
			colormap flag
			drawnow
		end
	end
	figure(3);
	imshow(x);
	colormap gray
	imwrite(x,'x_proposed.tif');
	
	figure(4)
	imshow(x);
	colormap gray
	hold on
	quiver(n1,n2+0.5,v(:,:,2,1),v(:,:,1,1),0,'r');
	quiver(n1+0.5,n2,v(:,:,2,2),v(:,:,1,2),0,'b');
	quiver(n1,n2,v(:,:,2,3),v(:,:,1,3),0,'g');
	
end


function xout = prox_mu_tau_f (x,idx,yidx)
 	xout =x;
 	xout(idx)=yidx;
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
