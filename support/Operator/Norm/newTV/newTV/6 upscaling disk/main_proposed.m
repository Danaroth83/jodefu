function main
	
	Nbiter= 10000;
	tau = 0.9/8;
	sigma = 0.9/3;
	mu = 1;

	S=23;
	level=6;
	facto=2^level;
	y=zeros(S*facto,S*facto);
	[n1,n2]=meshgrid(1:S*facto,1:S*facto);
	y=((n1-S*facto/2-1/2).^2+(n2-S*facto/2-1/2).^2<=(9*facto)^2);
	for lev=1:level-2  
		y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	end
	xref=y*3/4+1/8+0.0004;
	for lev=level-1:level
		y=(y(1:2:end,1:2:end)+y(2:2:end,1:2:end)+y(1:2:end,2:2:end)+y(2:2:end,2:2:end))/4;
	end
	y=y*3/4+1/8+0.0004;
	figure(1);
	imshow(xref);
	colormap gray
	imwrite(xref,'x_reference.tif');
	figure(2);
	imshow(y);
	colormap gray
	imwrite(y,'y.tif');
			
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,3))/(mu*sigma),1));

	x = prox_mu_tau_f(zeros(size(y)*4),y);
	u = zeros([size(x) 2]);
	v = zeros([size(x) 2 3]);
	tmp = opD(x);
	v(:,:,1,1) = tmp(:,:,1);
	v(:,:,2,2) = tmp(:,:,2);
	fprintf('0 %f\n',sum(abs(tmp(:)))); % == sum(sum(sum(sqrt(sum(v.^2,3)))))
	for iter = 1:Nbiter
		x = prox_mu_tau_f(x+tau*opDadj(-opD(x)+opLadj(v)-mu*u),y);
		v = prox_mu_sigma_g(v-sigma*opL(-opD(x)+opLadj(v)-mu*u));
		u = u-(-opD(x)+opLadj(v))/mu;
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			tmp = opDadj(u);
			tmp = tmp(1:4:end,1:4:end)+tmp(2:4:end,1:4:end)+tmp(3:4:end,1:4:end)+tmp(4:4:end,1:4:end)+...
			tmp(1:4:end,2:4:end)+tmp(2:4:end,2:4:end)+tmp(3:4:end,2:4:end)+tmp(4:4:end,2:4:end)+...
			tmp(1:4:end,3:4:end)+tmp(2:4:end,3:4:end)+tmp(3:4:end,3:4:end)+tmp(4:4:end,3:4:end)+...
			tmp(1:4:end,4:4:end)+tmp(2:4:end,4:4:end)+tmp(3:4:end,4:4:end)+tmp(4:4:end,4:4:end);
			fprintf('%d %f %f\n',iter,sum(sum(sum(sqrt(sum(v.^2,3))))),...
				sum(sum(tmp.*y)));
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
	fprintf('||x-xref||=%f\n',sqrt(sum(sum((x-xref).^2))));
	
	figure(4)
	imshow(x);
	colormap gray
	hold on
	[n1,n2]=meshgrid(1:S*4,1:S*4);
	quiver(n1,n2+0.5,v(:,:,2,1),v(:,:,1,1),0,'r');
	quiver(n1+0.5,n2,v(:,:,2,2),v(:,:,1,2),0,'b');
	quiver(n1,n2,v(:,:,2,3),v(:,:,1,3),0,'g');
	
end


function xout = prox_mu_tau_f(x,y)
	z=y-(x(1:4:end,1:4:end)+x(2:4:end,1:4:end)+x(3:4:end,1:4:end)+x(4:4:end,1:4:end)+...
	x(1:4:end,2:4:end)+x(2:4:end,2:4:end)+x(3:4:end,2:4:end)+x(4:4:end,2:4:end)+...
	x(1:4:end,3:4:end)+x(2:4:end,3:4:end)+x(3:4:end,3:4:end)+x(4:4:end,3:4:end)+...
	x(1:4:end,4:4:end)+x(2:4:end,4:4:end)+x(3:4:end,4:4:end)+x(4:4:end,4:4:end))/16;
	xout = x;
	xout(1:4:end,1:4:end)=x(1:4:end,1:4:end)+z;
	xout(2:4:end,1:4:end)=x(2:4:end,1:4:end)+z;
	xout(3:4:end,1:4:end)=x(3:4:end,1:4:end)+z;
	xout(4:4:end,1:4:end)=x(4:4:end,1:4:end)+z;
	xout(1:4:end,2:4:end)=x(1:4:end,2:4:end)+z;
	xout(2:4:end,2:4:end)=x(2:4:end,2:4:end)+z;
	xout(3:4:end,2:4:end)=x(3:4:end,2:4:end)+z;
	xout(4:4:end,2:4:end)=x(4:4:end,2:4:end)+z;
	xout(1:4:end,3:4:end)=x(1:4:end,3:4:end)+z;
	xout(2:4:end,3:4:end)=x(2:4:end,3:4:end)+z;
	xout(3:4:end,3:4:end)=x(3:4:end,3:4:end)+z;
	xout(4:4:end,3:4:end)=x(4:4:end,3:4:end)+z;
	xout(1:4:end,4:4:end)=x(1:4:end,4:4:end)+z;
	xout(2:4:end,4:4:end)=x(2:4:end,4:4:end)+z;
	xout(3:4:end,4:4:end)=x(3:4:end,4:4:end)+z;
	xout(4:4:end,4:4:end)=x(4:4:end,4:4:end)+z;	
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
