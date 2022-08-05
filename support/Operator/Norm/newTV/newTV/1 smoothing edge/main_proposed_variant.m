% Variant of the proposed approach with the 4 times, instead of 3 times, 
% finer grid for the gradient field. Other said, the gradient is defined 
% at pixel corners, in addition to pixel centers and pixel edges.

function main 
	
	Nbiter= 2000;
	lambda = 2; 
	tau = 0.99/8;
	sigma = 0.99/3;
	mu = 0.05; 
	
	S=256;
	[n1,n2]=meshgrid(1:S,1:S);
	y=double((n1-S/2-1/2)+(n2-S/2-1/2)*3.2+0.001>=0); 
	
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_mu_tau_f = @(x) (x+mu*tau*y)/(1+mu*tau);
	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,3))/(lambda*mu*sigma),1));
	
	figure(1);
	imshow(y);
	colormap gray
	imwrite(y,'y.tif');
		
	x = y;
	u = zeros([size(x) 2]);
	v = zeros([size(x) 2 4]);
	tmp = opD(x);
	v(:,:,1,1) = tmp(:,:,1);
	v(:,:,2,2) = tmp(:,:,2);
	
	fig2=figure(2);
	imshow(y(S/2-6:S/2+7,S/2-19:S/2+19));
	colormap gray
	hold on
	quiver(n1(1:14,1:39),n2(1:14,1:39)+0.5,v(S/2-6:S/2+7,S/2-19:S/2+19,2,1),...
		v(S/2-6:S/2+7,S/2-19:S/2+19,1,1),0,'r');
	quiver(n1(1:14,1:39)+0.5,n2(1:14,1:39),v(S/2-6:S/2+7,S/2-19:S/2+19,2,2),...
		v(S/2-6:S/2+7,S/2-19:S/2+19,1,2),0,'b');
	truesize(fig2,[15*24,39*24]);

	fprintf('0 %f\n',lambda*sum(sum(sum(sqrt(sum(v.^2,3))))));
	for iter = 1:Nbiter
		x = prox_mu_tau_f(x+tau*opDadj(-opD(x)+opLadj(v)-mu*u));
		v = prox_mu_sigma_g(v-sigma*opL(-opD(x)+opLadj(v)-mu*u));
		u = u-(-opD(x)+opLadj(v))/mu;
		if mod(iter,40)==0
			fprintf('%d %f\n',iter,sum(sum((x-y).^2))/2+lambda*sum(sum(sum(sqrt(sum(v.^2,3))))));
			figure(3);
			imshow(x);
			colormap flag
			drawnow
		end
	end
	figure(3);
	imshow(x);
	colormap gray
	imwrite(x,'x_proposed_variant.tif');
	
	figure(41);
	imshow(x);
	colormap gray
	hold on
	quiver(n1,n2+0.5,v(:,:,2,1),v(:,:,1,1),0,'r');
	quiver(n1+0.5,n2,v(:,:,2,2),v(:,:,1,2),0,'b');
	quiver(n1,n2,v(:,:,2,3),v(:,:,1,3),0,'g');
	quiver(n1+0.5,n2+0.5,v(:,:,2,4),v(:,:,1,4),0,'y');
	
	fig5=figure(51);
	imshow(x(S/2-6:S/2+7,S/2-19:S/2+19)); 
	colormap gray
	hold on
	quiver(n1(1:14,1:39),n2(1:14,1:39)+0.5,v(S/2-6:S/2+7,S/2-19:S/2+19,2,1),...
		v(S/2-6:S/2+7,S/2-19:S/2+19,1,1),0,'r');
	quiver(n1(1:14,1:39)+0.5,n2(1:14,1:39),v(S/2-6:S/2+7,S/2-19:S/2+19,2,2),...
		v(S/2-6:S/2+7,S/2-19:S/2+19,1,2),0,'b');
	quiver(n1(1:14,1:39),n2(1:14,1:39),v(S/2-6:S/2+7,S/2-19:S/2+19,2,3),...
		v(S/2-6:S/2+7,S/2-19:S/2+19,1,3),0,'g');
	quiver(n1(1:14,1:39)+0.5,n2(1:14,1:39)+0.5,v(S/2-6:S/2+7,S/2-19:S/2+19,2,4),...
		v(S/2-6:S/2+7,S/2-19:S/2+19,1,4),0,'y');
	truesize(fig5,[15*24,39*24]);
end


function t = opL(u)
	t=zeros([size(u) 4]);
	t(:,:,1,1)=u(:,:,1); 
	t(1:end-1,2:end,2,1)=(u(2:end,1:end-1,2)+u(1:end-1,1:end-1,2)+...
		u(2:end,2:end,2)+u(1:end-1,2:end,2))/4;
	t(1:end-1,1,2,1)=(u(1:end-1,1,2)+u(2:end,1,2))/4;
	t(:,:,2,2)=u(:,:,2);
	t(2:end,1:end-1,1,2)=(u(2:end,1:end-1,1)+u(1:end-1,1:end-1,1)+...
		u(2:end,2:end,1)+u(1:end-1,2:end,1))/4;
	t(1,1:end-1,1,2)=(u(1,1:end-1,1)+u(1,2:end,1))/4;
	t(2:end,:,1,3) = (u(2:end,:,1)+u(1:end-1,:,1))/2;
	t(1,:,1,3) = u(1,:,1)/2;
	t(:,2:end,2,3) = (u(:,2:end,2)+u(:,1:end-1,2))/2;
	t(:,1,2,3) = u(:,1,2)/2;
	t(1:end-1,1:end-1,1,4) = (u(1:end-1,1:end-1,1)+u(1:end-1,2:end,1))/2;
	t(1:end-1,1:end-1,2,4) = (u(1:end-1,1:end-1,2)+u(2:end,1:end-1,2))/2;
end

function u = opLadj(t)
	[height,width,d,c]=size(t);
	u=zeros(height,width,2);
	u(1:end-1,2:end,1)=t(1:end-1,2:end,1,1)+(t(1:end-1,2:end,1,2)+...
		t(1:end-1,1:end-1,1,2)+t(2:end,2:end,1,2)+t(2:end,1:end-1,1,2))/4+...
		(t(1:end-1,2:end,1,3)+t(2:end,2:end,1,3))/2+...
		(t(1:end-1,1:end-1,1,4)+t(1:end-1,2:end,1,4))/2;
	u(1:end-1,1,1)=t(1:end-1,1,1,1)+(t(1:end-1,1,1,2)+t(2:end,1,1,2))/4+...
		(t(1:end-1,1,1,3)+t(2:end,1,1,3))/2+t(1:end-1,1,1,4)/2;
	u(2:end,1:end-1,2)=t(2:end,1:end-1,2,2)+(t(2:end,1:end-1,2,1)+...
		t(1:end-1,1:end-1,2,1)+t(2:end,2:end,2,1)+t(1:end-1,2:end,2,1))/4+...
		(t(2:end,1:end-1,2,3)+t(2:end,2:end,2,3))/2+...
		(t(1:end-1,1:end-1,2,4)+t(2:end,1:end-1,2,4))/2;
	u(1,1:end-1,2)=t(1,1:end-1,2,2)+(t(1,1:end-1,2,1)+t(1,2:end,2,1))/4+...
		(t(1,1:end-1,2,3)+t(1,2:end,2,3))/2+t(1,1:end-1,2,4)/2;
end	