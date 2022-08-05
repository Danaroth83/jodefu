function main
	
	Nbiter= 6000;
	tau = 0.02;
	sigma = 1/tau/16;
	rho = 1.9;

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
	
	opD = @(x) cat(3,[-diff(x,1,1);zeros(1,size(x,2))],[-diff(x,1,2) zeros(size(x,1),1)],...
		[zeros(1,size(x,2));diff(x,1,1)],[zeros(size(x,1),1) diff(x,1,2)]);
	opDadj = @(u) [u(1,:,1);diff(u(:,:,1),1,1)]+[u(:,1,2) diff(u(:,:,2),1,2)]+...
		[-diff(u(:,:,3),1,1);u(end,:,3)]+[-diff(u(:,:,4),1,2) u(:,end,4)];	
	
	x = prox_tau_f(zeros(size(y)*4),y);
	u = zeros([size(x) 4]);
	fprintf('0 %f\n',sum(sum(sqrt(sum(max(opD(x),0).^2,3)))));
	for iter = 1:Nbiter
		xnew = prox_tau_f(x-tau*opDadj(u),y);
		unew = prox_sigma_g_conj(u+sigma*opD(2*xnew-x));
		x = xnew+(rho-1)*(xnew-x);
		u = unew+(rho-1)*(unew-u);
		if mod(iter,40)==0
			%we display the primal and dual cost functions, which reach equal values at convergence
			tmp = opDadj(unew);
			tmp = tmp(1:4:end,1:4:end)+tmp(2:4:end,1:4:end)+tmp(3:4:end,1:4:end)+tmp(4:4:end,1:4:end)+...
			tmp(1:4:end,2:4:end)+tmp(2:4:end,2:4:end)+tmp(3:4:end,2:4:end)+tmp(4:4:end,2:4:end)+...
			tmp(1:4:end,3:4:end)+tmp(2:4:end,3:4:end)+tmp(3:4:end,3:4:end)+tmp(4:4:end,3:4:end)+...
			tmp(1:4:end,4:4:end)+tmp(2:4:end,4:4:end)+tmp(3:4:end,4:4:end)+tmp(4:4:end,4:4:end);
			fprintf('%d %f %f\n',iter,sum(sum(sqrt(sum(max(opD(xnew),0).^2,3)))),...
				sum(sum(tmp.*y)));
			figure(3);
			imshow(xnew);
			colormap flag
			drawnow
		end
	end
	figure(3);
	imshow(xnew);
	colormap gray
	imwrite(xnew,'x_upwind.tif');
	fprintf('||x-xref||=%f\n',sqrt(sum(sum((xnew-xref).^2))));
	
end

function xout = prox_tau_f(x,y)
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

function unew = prox_sigma_g_conj(u) 
	u = max(u,0);
	unew = bsxfun(@rdivide, u, max(sqrt(sum(u.^2,3)),1));
end