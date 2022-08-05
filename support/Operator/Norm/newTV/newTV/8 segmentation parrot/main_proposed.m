function main

	lambda = 0.09;	% regularization parameter
	nbiter = 1000;	% number of iterations
	tau = 0.99/8;
	sigma = 0.99/3;
	mu = 50; 

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
	
	opD = @(z) cat(4,cat(2,diff(z,1,2),zeros(K,1,width)),cat(3,diff(z,1,3),zeros(K,height,1)));
 	opDadj = @(u) -cat(2,u(:,1,:,1),diff(u(:,:,:,1),1,2))-cat(3,u(:,:,1,2),diff(u(:,:,:,2),1,3));	
 	proj_simplex = @(z) max(bsxfun(@minus,z,max(bsxfun(@rdivide,...
 		cumsum(sort(z,1,'descend'),1)-1,(1:K)'))),0);
 	prox_mu_sigma_g = @(t) t-bsxfun(@rdivide, t, max(sqrt(sum(t.^2,4))/(lambda*mu*sigma),1));

	label = reshape(label,K,1,1,3);
	p = sum(bsxfun(@minus, shiftdim(y,-1), label).^2,4);
	[tmp, idx] = min(p);
    z = double(bsxfun(@eq, (1:K)', idx));
    % This initialization assigns the closest label to every pixel. 
    % Equivalently, this is the solution of the problem when lambda=0.
    figure(2)
	imshow(shiftdim(sum(bsxfun(@times, z, label),1),1));
 	u = zeros(K,height,width,2);
 	v = zeros(K,height,width,2,3);
	opDz = opD(z);
	v(:,:,:,1,1) = opDz(:,:,:,1);
	v(:,:,:,2,2) = opDz(:,:,:,2);
 	opLadjv = opLadj(v);
 	fprintf('0 %f\n',(sum(z(:).*p(:)) + lambda*sum(abs(opDz(:))))/2);	
	for iter = 1:nbiter		
		z = proj_simplex(z - tau*opDadj(opDz-opLadjv+mu*u) - mu*tau*p);
		opDz = opD(z);
		v = prox_mu_sigma_g(v + sigma*opL(opDz-opLadjv+mu*u));
		opLadjv = opLadj(v);
		u = u + (opDz-opLadjv)/mu;
    	if mod(iter,10)==0
			fprintf('%d %f %f\n',iter,(sum(z(:).*p(:)) + ...
				lambda*sum(sum(sum(sum(sqrt(sum(v.^2,4)))))))/2,...
				sum(sum(min(opDadj(u)+p)))/2);
			x = shiftdim(sum(bsxfun(@times,z,label),1),1);
			figure(3)
			imshow(x);
		end
	end
	x = shiftdim(sum(bsxfun(@times,z,label),1),1);
	figure(3)
	imshow(x);
	imwrite(x,'x_proposed.tif');
	
end    
    

function t = opL(u)
	[K,height,width,d]=size(u);
	t=zeros(K,height,width,2,3);	
	t(:,:,:,1,1)=u(:,:,:,1); 
	t(:,1:end-1,2:end,2,1)=(u(:,2:end,1:end-1,2)+u(:,1:end-1,1:end-1,2)+u(:,2:end,2:end,2)+u(:,1:end-1,2:end,2))/4;
	t(:,1:end-1,1,2,1)=(u(:,1:end-1,1,2)+u(:,2:end,1,2))/4; 	
	t(:,:,:,2,2)=u(:,:,:,2); 
	t(:,2:end,1:end-1,1,2)=(u(:,2:end,1:end-1,1)+u(:,1:end-1,1:end-1,1)+u(:,2:end,2:end,1)+u(:,1:end-1,2:end,1))/4; 
	t(:,1,1:end-1,1,2)=(u(:,1,1:end-1,1)+u(:,1,2:end,1))/4; 	
	t(:,2:end,:,1,3) = (u(:,2:end,:,1)+u(:,1:end-1,:,1))/2; % correspond au coset aux centre des pixels
	t(:,1,:,1,3) = u(:,1,:,1)/2;
	t(:,:,2:end,2,3) = (u(:,:,2:end,2)+u(:,:,1:end-1,2))/2; 
	t(:,:,1,2,3) = u(:,:,1,2)/2;
end

function u = opLadj(t)
	[K,height,width,d,c]=size(t);
	u=zeros(K,height,width,2);
	u(:,1:end-1,2:end,1)=t(:,1:end-1,2:end,1,1)+(t(:,1:end-1,2:end,1,2)+t(:,1:end-1,1:end-1,1,2)+...
	t(:,2:end,2:end,1,2)+t(:,2:end,1:end-1,1,2))/4+(t(:,1:end-1,2:end,1,3)+t(:,2:end,2:end,1,3))/2;
	u(:,1:end-1,1,1)=t(:,1:end-1,1,1,1)+(t(:,1:end-1,1,1,2)+t(:,2:end,1,1,2))/4+...
	(t(:,1:end-1,1,1,3)+t(:,2:end,1,1,3))/2;
	u(:,2:end,1:end-1,2)=t(:,2:end,1:end-1,2,2)+(t(:,2:end,1:end-1,2,1)+t(:,1:end-1,1:end-1,2,1)+...
	t(:,2:end,2:end,2,1)+t(:,1:end-1,2:end,2,1))/4+(t(:,2:end,1:end-1,2,3)+t(:,2:end,2:end,2,3))/2;
	u(:,1,1:end-1,2)=t(:,1,1:end-1,2,2)+(t(:,1,1:end-1,2,1)+t(:,1,2:end,2,1))/4+...
	(t(:,1,1:end-1,2,3)+t(:,1,2:end,2,3))/2;
end
