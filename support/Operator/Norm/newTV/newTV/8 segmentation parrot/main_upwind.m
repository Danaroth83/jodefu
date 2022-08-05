function main

	lambda = 0.09;	% regularization parameter
	nbiter = 2000;	% number of iterations
	tau = 2;
	sigma = 1/tau/16; 
	rho = 1.9;
	
	K = 6;	% number of color labels
	label = zeros(K,3);
	label(1,:) = [50 50 50];      % black
	label(2,:) = [83 120 44];     % green
	label(3,:) = [115 110 80];    % brown
	label(4,:) = [91 147 158];    % turquoise
	label(5,:) = [235 188 12];    % yellow
	label(6,:) = [220 212 202];   % white
	label = label/255;
	y  = double(imread('parrot.png'))/255;   % initial image
	[height,width,nbchannels] = size(y);
	figure(1)
	imshow(y);
	
	opD = @(z) cat(4,cat(2,-diff(z,1,2),zeros(K,1,width)),cat(3,-diff(z,1,3),zeros(K,height,1)),...
		cat(2,zeros(K,1,width),diff(z,1,2)),cat(3,zeros(K,height,1),diff(z,1,3)));
	opDadj = @(u) cat(2,u(:,1,:,1),diff(u(:,:,:,1),1,2))+cat(3,u(:,:,1,2),diff(u(:,:,:,2),1,3))+...
		cat(2,-diff(u(:,:,:,3),1,2),u(:,end,:,3))+cat(3,-diff(u(:,:,:,4),1,3),u(:,:,end,4));
	proj_simplex = @(z) max(bsxfun(@minus,z,max(bsxfun(@rdivide,...
 		cumsum(sort(z,1,'descend'),1)-1,(1:K)'),[],1)),0);
 	
 	label = reshape(label,K,1,1,3);
	p = sum(bsxfun(@minus, shiftdim(y,-1), label).^2,4);
	[tmp, idx] = min(p,[],1);
    z = double(bsxfun(@eq, (1:K)', idx));	
    % This initialization assigns the closest label to every pixel. 
    % Equivalently, this is the solution of the problem when lambda=0.
    figure(2)
	imshow(shiftdim(sum(bsxfun(@times, z, label),1),1));
 	u = zeros(K,height,width,4);  	
	fprintf('0 %f\n',(sum(z(:).*p(:)) + lambda*sum(sum(sum(sqrt(sum(max(opD(z),0).^2,4))))))/2);		
	for iter = 1:nbiter		
		znew = proj_simplex(z - tau*opDadj(u) - tau*p);
		unew = prox_sigma_g_conj(u + sigma*opD(2*znew-z),lambda);
    	z = znew + (rho-1)*(znew - z);
    	u = unew + (rho-1)*(unew - u);
		if mod(iter,10)==0
			fprintf('%d %f %f\n',iter,(sum(znew(:).*p(:)) + ...
				lambda*sum(sum(sum(sqrt(sum(max(opD(znew),0).^2,4))))))/2,...
				sum(sum(min(opDadj(unew)+p,[],1)))/2);	
			x = shiftdim(sum(bsxfun(@times,znew,label),1),1);
			figure(3)
			imshow(x);
		end
	end
	x = shiftdim(sum(bsxfun(@times,znew,label),1),1);
	figure(3)
	imshow(x);
	imwrite(x,'x_upwind.tif');
	
end    

function unew = prox_sigma_g_conj(u,lambda) 
	u = max(u,0);
	unew = bsxfun(@rdivide, u, max(sqrt(sum(u.^2,4))/lambda,1));
end
 