function main
	
	Nbiter= 6000;
	tau = 0.1;
	sigma = 1/tau/8;
	rho = 1.9;

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
	prox_sigma_g_conj = @(u) bsxfun(@rdivide,u,max(sqrt(sum(u.^2,3)),1));

	x = y;
	u = zeros([size(x) 2]);
	opDy = opD(y);
	fprintf('0 %f\n',sum(sum(sqrt(sum(opDy.^2,3)))));
	opDy = opDy(:);
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u),idx2,yidx);
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x));
		x = xnew+(rho-1)*(xnew-x);
		u = unew+(rho-1)*(unew-u);
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			fprintf('%d %f %f\n',iter,sum(sum(sqrt(sum(opD(xnew).^2,3)))),...
				sum(unew(:).*opDy)); % == sum(sum(opDadj(unew).*y))
			figure(3);
			imshow(xnew);
			colormap flag
			drawnow
		end
	end
	figure(3);
	imshow(xnew);
	colormap gray
	imwrite(xnew,'x_isotropic.tif');
	
end


function xout = prox_tau_f(x,idx,yidx)
 	xout = x ;
 	xout(idx) = yidx;
end