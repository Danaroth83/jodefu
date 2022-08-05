function main
	
	Nbiter= 1600;
	lambda = 0.16; 
	tau = 0.01;
	sigma = 1/tau/8;
	rho = 1.9;

	rng(0);
	y=double(imread('bike2.tif'))/255;
	y=y+randn(size(y))*0.18;
	figure(1);
	imshow(y);
	colormap gray
	imwrite(y,'y.tif');
			
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_tau_f = @(x) (x+tau*y)/(1+tau);
	prox_sigma_g_conj = @(u) u./max(abs(u)/lambda,1);

	x = y;
	u = zeros([size(x) 2]);
	fprintf('0 %f\n',lambda*sum(sum(sum(abs(opD(y))))));
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u));
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x));
		x = xnew+(rho-1)*(xnew-x);
		u = unew+(rho-1)*(unew-u);
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			fprintf('%d %f %f\n',iter,sum(sum((xnew-y).^2))/2+lambda*sum(sum(sum(abs(opD(xnew))))),...
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
	imwrite(xnew,'x_anisotropic.tif');
	
end
