%% S1l1 Proximal operator
%
% Description:
% This function calculates the proximal opertator for the S1l1 norm, where
% S1 defines the Shatten-1 norm (also  known as nuclear norm) and is applied
% on the first two dimensions and the l1 norm is applied on the third
% one.
% If any of the first two dimensions has size 2, this algorithm
% implements the procedure described in section 4.3.2 of [1], else it
% performs soft thresholding of the singular values each slice of the input
%
% Usage:
% y=prox_S1l1(x,gamma);
% Example:
% y=prox_S1l1(randn(2,50,10000),0.1);
%
% Inputs:
% x: An Nd Matrix, the S1 proximal operator will be applied on its first
%    two dimensions, while the l1 is applied on the remaining ones
% gamma: multiplicative non-negative factor for the norm (soft-thresholding 
%        level of the singular values)
%
% Output:
% y: S1l1 norm computation
%
% References:
% [1] J. Duran, M. Moeller, C. Sbert, and D. Cremers, "On the
% Implementation of Collaborative TV Regularization: Application to
% Cartoon+Texture Decomposition",  IPOL, 6 (2016), pp. 27–74.

function y=prox_S1l1(x,gamma)

%% Reshaping matrix to be 3D

sizex=size(x);
x=reshape(x,size(x,1),size(x,2),[]);

flag_perm=0;
if size(x,1) > size(x,2), flag_perm=1; x=permute(x,[2,1,3]); end

if size(x,1)>2
        
    % Implementation for dimensions' sizes above 2 (works for 2 as well, but it's slower than below)
    
    y=zeros(size(x));
    
    for ii=1:size(x,3)
        A=x(:,:,ii);
        [U,Sigma,V]=svd(A,'econ');
        diag_Sigma = diag(Sigma);
        prox_Sigma = prox_l1(diag_Sigma,gamma);
        y(:,:,ii)=U*diag(prox_Sigma)*V';
    end
    
else
    
    fTiny = 0.00000001;

    [num_derivatives,num_channels,num_pixels] = size(x);
    y = zeros(num_derivatives,num_channels,num_pixels);

    for ii = 1:num_pixels

            valuex=squeeze(x(1,:,ii));
            valuey=squeeze(x(2,:,ii));
            M1 = sum(valuex.^2);
            M2 = sum(valuex.*valuey);
            M3 = sum(valuey.^2);

            T = M1 + M3;
            D = M1 * M3 - M2 * M2;
            MyDet = sqrt(max((T^2 / 4) - D, 0));
            eig1 = max((T / 2) + MyDet, 0);
            eig2 = max((T / 2) - MyDet, 0);
            sig1 = sqrt(eig1);
            sig2 = sqrt(eig2);

            V1 = 0; V2 = 0; V3 = 0; V4 = 0; 

            if (M2 ~= 0)

                v0 = M2;
                v1 = eig1 - M3;
                v2 = eig2 - M3;
                mu1 = sqrt(v0^2 + v1^2);
                mu2 = sqrt(v0^2 + v2^2);

                if(mu1 > fTiny)
                    V1 = v1 / mu1;
                    V3 = v0 / mu1;
                end

                if(mu2 > fTiny)
                    V2 = v2 / mu2;
                    V4 = v0 / mu2;
                end
                
            else
                if(M1 > M3)
                    V1 = 1; V4 = 1;
                    V2 = 0; V3 = 0;
                else
                    V1 = 0; V4 = 0;
                    V2 = 1; V3 = 1;
                end
            end

            %% Compute prox_p of the diagonal entries

            sig1_upd = max(sig1 - gamma, 0);
            sig2_upd = max(sig2 - gamma, 0);

            %% Compute the diagonal entries of $\widehat{\Sigma}\Sigma^{\dagger}_0$
            if(sig1 > fTiny)
                sig1_upd = sig1_upd/sig1;
            end
            if(sig2 > fTiny)
                sig2_upd = sig2_upd/sig2;
            end
            %% Compute solution

            t1 = sig1_upd * V1^2 + sig2_upd * V2^2;
            t2 = sig1_upd * V1 * V3 + sig2_upd * V2 * V4;
            t3 = sig1_upd * V3^2 + sig2_upd * V4^2;

            y(1,:,ii)= valuex * t1 + valuey * t2;
            y(2,:,ii)= valuex * t2 + valuey * t3;
    end
    
end

if flag_perm==1, y=ipermute(y,[2,1,3]); end
y=reshape(y,sizex);

end


% %% Alternative implementation by Condat (slower, only valid if size of one of the first two dimension is 2)
% Based on calculating the svd unitary matrix with the method proposed 
% in http://scipp.ucsc.edu/~haber/ph116A/diag2x2_11.pdf
%
% function y=prox_S1l1(x,gamma)
% 
%     sizex=size(x);
%     x=reshape(x,size(x,1),size(x,2),[]);
% 
%     flag_perm=0;
%     if size(x,1) > size(x,2), flag_perm=1; x=permute(x,[2,1,3]); end
%     y=zeros(size(x));
%     
%     for ii=1:size(x,3)
%         in=x(:,:,ii);
%         
%         s = sum(in.^2,2);
%         theta = atan2(2*sum(prod(in,1),2),-diff(s))/2;
%         c = cos(theta);
%         s = sin(theta);
%         v = [c -s; s c];
%         out = v' * in;
%         tmp = max(sqrt(sum(out.^2,2)), gamma);
%         out = v * (out .* (tmp-gamma)./tmp);
%         
%         y(:,:,ii) = out;
%     end
%     if flag_perm==1, y=ipermute(y,[2,1,3]); end
%     y=reshape(y,sizex);
% end