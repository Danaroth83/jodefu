function main
	
	Nbiter= 1000;
	lambda = 0.16; 
	tau = 0.01;
	sigma = 1/tau/16;
	rho = 1.9;

	rng(0);
	y=double(imread('bike2.tif'))/255;
	y=y+randn(size(y))*0.18;
	imwrite(y,'y.tif');
		
	opD = @(x) cat(3,[-diff(x,1,1);zeros(1,size(x,2))],[-diff(x,1,2) zeros(size(x,1),1)],...
		[zeros(1,size(x,2));diff(x,1,1)],[zeros(size(x,1),1) diff(x,1,2)]);
	opDadj = @(u) [u(1,:,1);diff(u(:,:,1),1,1)]+[u(:,1,2) diff(u(:,:,2),1,2)]+...
		[-diff(u(:,:,3),1,1);u(end,:,3)]+[-diff(u(:,:,4),1,2) u(:,end,4)];	
	prox_tau_f = @(x) (x+tau*y)/(1+tau);
	
	x = y;
	u = zeros([size(x) 4]);
	fprintf('0 %f\n',lambda*sum(sum(sqrt(sum(max(opD(y),0).^2,3)))));
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u));
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x),lambda);
		x = xnew+(rho-1)*(xnew-x);
		u = unew+(rho-1)*(unew-u);
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			fprintf('%d %f %f\n',iter,sum(sum((xnew-y).^2))/2+...
				lambda*sum(sum(sqrt(sum(max(opD(xnew),0).^2,3)))),...
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
	imwrite(xnew,'x_upwind.tif');
	
end

function unew = prox_sigma_g_conj(u,lambda) 
	u = max(u,0);
	unew = bsxfun(@rdivide, u, max(sqrt(sum(u.^2,3))/lambda,1));
end