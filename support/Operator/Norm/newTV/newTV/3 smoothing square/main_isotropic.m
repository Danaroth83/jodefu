function main
	
	Nbiter= 3000;
	lambda = 6; 
	tau = 0.01;
	sigma = 1/tau/8;
	rho = 1.9;

	S=100;
	[n1,n2]=meshgrid(1:S,1:S);
	y=double((n1>=19)&(n2>=19)&(n1<=82)&(n2<=82));
	figure(1);
	imshow(y);
	colormap gray
	imwrite(y,'y.tif');
			
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_tau_f = @(x) (x+tau*y)/(1+tau);
	prox_sigma_g_conj = @(u) bsxfun(@rdivide, u, max(sqrt(sum(u.^2,3))/lambda,1));
	
	x = y;
	u = zeros([size(x) 2]);
	fprintf('0 %f\n',lambda*sum(sum(sqrt(sum(opD(y).^2,3)))));
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u));
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x));
		x = xnew+(rho-1)*(xnew-x);
		u = unew+(rho-1)*(unew-u);
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			fprintf('%d %f %f\n',iter,sum(sum((xnew-y).^2))/2+lambda*sum(sum(sqrt(sum(opD(xnew).^2,3)))),...
				-sum(sum((y-opDadj(unew)).^2-y.^2))/2);
			figure(2);
			imshow(xnew);
			colormap flag
			drawnow
		end
	end
	figure(2);
	imshow(xnew);
	colormap gray
	imwrite(xnew,'x_isotropic.tif');
	
end
