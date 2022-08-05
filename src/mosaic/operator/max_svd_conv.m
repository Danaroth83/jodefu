%% Maximum Singular value of a Convolution Matrix

% Description:
% Calculates the maximum Singular Value of a Convolution Matrix, either
% said a matrix in Toeplitz form that can be seen as a circulant matrix,
% given all non zero elements of the first row. The implementation is based
% on the principle that svd of a circ
%
% Usage:
% [maxsvd,svd_val]=max_svd_conv(vect)
%
% Input:
% vect: A vector containing the first row of non zero elements of the 
%       circulant matrix (if not a vector, the first two dimension are
%       concatenated)
%
% Output:
% maxsvd:  The maximum singular value of the circulant matrix
% svd_val: All the singular values of the circulant matrix
%
% Reference:
% [1] http://pwp.gatech.edu/wp-content/uploads/sites/436/2011/04/circsvd-notes.pdf


function [maxsvd,svd_val]=max_svd_conv(vect)

    vect=reshape(vect,size(vect,1)*size(vect,2),[]);

    svd_val=abs(fft(vect));
    maxsvd=max(svd_val(:));

end



% %% Original implementation
% 
% function [maxeig,eig]=max_eig_easy(vect)
% 
% %Estimates the maximum eigenvalue of a A^H*A, where A is a circulant matrix
% %whose non-zero elements are listed in vect.
% % The script exploits the formula for the calculus of eigenvalues with
% % root of unity 
% 
% vect=vect(:);
% 
% %Removing quasi-zero coefficients
% %tol=0.0001;
% %max_vect=max(abs(vect));
% %vect(vect>-tol*max_vect & vect<tol*max_vect)=[];
% 
% 
% N=length(vect);
% a=xcorr(vect);
% a=a(N:end);
% 
% h=1:N-1;
% maxden=1000;
% ratio=(0:ceil(maxden/2))/maxden;
% % ratio=(0:maxden-1)/maxden; %unuseful as cosine is an even function
% 
% [x,y]=meshgrid(ratio,h);
% len_r=length(ratio);
% w=[ones(1,len_r);2*cos(2*pi*x.*y)];
% 
% eig=sum(repmat(a,[1,len_r]).*w,1);
% maxeig=max(eig);
% 
% end

% %% Original implementation
% Taken from http://pwp.gatech.edu/wp-content/uploads/sites/436/2011/04/circsvd-notes.pdf

% function [U, Lambda, V] = conv_svd(h)
%     n = length(h);
%     hhat = fft(h(:));
%     V = zeros(n);
%     V(:,1) = sqrt(1/n)*ones(n,1);
%     V(:,2:n/2) = sqrt(2/n)*cos(2*pi/n*(0:n-1)’*(1:n/2-1));
%     V(:,n/2+1) = masqrt(1/n)*(-1).^((0:n-1)’);
%     V(:,n/2+2:n) = sqrt(2/n)*sin(2*pi/n*(0:n-1)’*(n/2+1:n-1));
%     U = zeros(n);
%     U(:,1) = sqrt(1/n)*sign(hhat(1))*ones(n,1);
%     U(:,2:n/2) = sqrt(2/n)*cos(2*pi/n*(0:n-1)’*(1:n/2-1) + ones(n,1)*angle(hhat(2:n/2))’);
%     U(:,n/2+1) = sqrt(1/n)*sign(hhat(n/2+1))*(-1).^((0:n-1)’);
%     U(:,n/2+2:n) = sqrt(2/n)*sin(2*pi/n*(0:n-1)’*(n/2+1:n-1) + ones(n,1)*angle(hhat(n/2+2:n))’);
%     Lambda = diag(abs(hhat));
% end

