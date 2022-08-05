% Variant of the proposed approach with the 4 times, instead of 3 times, 
% finer grid for the gradient field. Other said, the gradient is defined 
% at pixel corners, in addition to pixel centers and pixel edges.

function main 
	
	Nbiter= 2000;
	sigma = 0.99/3;
	mu = 0.5; 
	
	S=32;
	[n1,n2]=meshgrid(1:S,1:S);
	
	Pattern = 3;
	switch Pattern
		case 1, x=double((n1-S/2)>=1);
		case 2, x=double((n1-S/2)+(n2-S/2)>=1);
		case 3,	x=double((n1-S/2)+(n2-S/2)>=1);
			idx=find((n1-S/2)+(n2-S/2)==1);
			x(idx)= 0.5;
		case 4, x=double((n1-S/2)+(n2-S/2)>=1);
			idx=find((n1-S/2)+(n2-S/2)==1);
			x(idx)= 7/8;
			idx=find((n1-S/2)+(n2-S/2)==0);
			x(idx)= 1/8;
		case 5, x=double((n1-S/2)==1); 
		case 6,	x=double((n1-S/2)+(n2-S/2)==1);
		case 7, x=zeros(S);
			idx=find((n1-S/2)+(n2-S/2)==1);
			x(idx)= 1; 
			idx=find((n1-S/2)+(n2-S/2)==0);
			x(idx)= 0.5; 
			idx=find((n1-S/2)+(n2-S/2)==2);
			x(idx)= 0.5; 
		case 8, x=zeros(S); x(S/2,S/2)=1;
		case 9, x=double((-1).^n1==1);
		case 10, x=double((-1).^(n1+n2)==1);
		case 11, x=zeros(S); x(S/2,S/2)=1; x(S/2-1,S/2)=0.5;
			x(S/2+1,S/2)=0.5; x(S/2,S/2+1)=0.5; x(S/2,S/2-1)=0.5;
	end
	
	fig1=figure(11); 
	imshow(x);
	truesize(fig1,[512,512]);
	
	opDx = cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,3))/(mu*sigma),1));
		
	u = zeros([size(x) 2]);
	v = zeros([size(x) 2 4]);
	tmp = opDx;
	v(:,:,1,1) = tmp(:,:,1);
	v(:,:,2,2) = tmp(:,:,2);
	fprintf('0 %f\n',sum(sum(sum(sqrt(sum(v.^2,3))))));
	for iter = 1:Nbiter
		v = prox_mu_sigma_g(v-sigma*opL(-opDx+opLadj(v)-mu*u));
		u = u-(-opDx+opLadj(v))/mu;
		if mod(iter,400)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			fprintf('%d %f %f\n',iter,sum(sum(sum(sqrt(sum(v.^2,3))))),...
				sum(u(:).*opDx(:)));  % == sum(sum(opDadj(u).*x)));
		end
	end
	
	figure(11)
	hold on
	quiver(n1,n2+0.5,v(:,:,2,1),v(:,:,1,1),0,'r');
	quiver(n1+0.5,n2,v(:,:,2,2),v(:,:,1,2),0,'b');
	quiver(n1,n2,v(:,:,2,3),v(:,:,1,3),0,'g');
	quiver(n1+0.5,n2+0.5,v(:,:,2,4),v(:,:,1,4),0,'y');

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