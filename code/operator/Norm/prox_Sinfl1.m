%% S1linf Proximal operator
%
% Description:
% This function calculates the proximal opertator for the S_infl_1 norm,where
% S_inf defines the Shatten-inf norm (choosing the maximum singular value) 
% and is applied on the first two dimensions and the l1 norm is applied on
% the third one.
% 
% If any of the first two dimensions has size 2, this algorithm
% implements the procedure described in section 4.3.2 of [1], else it
% applies the l_inf,1 on the singular values of each slice
%
% Usage:
% y=prox_S1linf(x,gamma);
% Example:
% y=prox_S1linf(randn(2,50,10000),0.1);
%
% Inputs:
% x: An Nd Matrix, the S1 proximal operator will be applied on its first
%    two dimensions, while the linf is applied on the third one
% gamma: multiplicative non-negative factor for the norm
%
% Output:
% y: S1linf norm computation
%
% References:
% [1] J. Duran, M. Moeller, C. Sbert, and D. Cremers, "On the
% Implementation of Collaborative TV Regularization: Application to
% Cartoon+Texture Decomposition",  IPOL, 6 (2016), pp. 27–74.

function y=prox_Sinfl1(x,gamma)

    sizex=size(x);
    x=reshape(x,size(x,1),size(x,2),[]);

    flag_perm=0;
    if size(x,1) > size(x,2), flag_perm=1; x=permute(x,[2,1,3]); end

    if size(x,1)>2
        
        y=zeros(size(x));

        for ii=1:size(x,3)
            [U,Sigma,V]=svd(x(:,:,ii),'econ');
            diag_Sigma=diag(Sigma);
            prox_Sigma=prox_linf1(diag_Sigma,gamma);
            y(:,:,ii)=U*diag(prox_Sigma)*V';
        end

    else

        fTiny = 0.00000001;
        sigma  = 1/gamma;

        [num_derivatives,num_channels,num_pixels] = size(x);
        y = zeros(num_derivatives,num_channels,num_pixels);

        for ii = 1:num_pixels

            valuex = squeeze(x(1,:,ii));
            valuey = squeeze(x(2,:,ii));

            M1=sum(valuex.^2);
            M2=sum(valuex.*valuey);
            M3=sum(valuey.^2);

            T = M1 + M3;
            D = M1 * M3 - M2 * M2;
            MyDet = sqrt(max((T * T / 4) - D, 0));
            eig1 = max((T / 2) + MyDet, 0);
            eig2 = max((T / 2) - MyDet, 0);
            sig1 = sqrt(eig1);
            sig2 = sqrt(eig2);
            V1 = 0; V2 = 0; V3 = 0; V4 = 0; 

            if (M2 ~= 0)

                v0 = M2;
                v1 = eig1 - M3;
                v2 = eig2 - M3;

                mu1 = sqrt(v0 * v0 + v1 * v1);
                mu2 = sqrt(v0 * v0 + v2 * v2);

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

            proj = zeros(2,1);         
            proj(1) = sigma * abs(sig1);
            proj(2) = sigma * abs(sig2);
            proj = l1projection(proj);

            sig1_upd = sig1 - gamma * proj(1);
            sig2_upd = sig2 - gamma * proj(2);

            %% Compute the diagonal entries of $\widehat{\Sigma}\Sigma^{\dagger}_0$

            if(sig1 > fTiny)
                sig1_upd = sig1_upd/sig1;
            end
            if(sig2 > fTiny)
                sig2_upd = sig2_upd/sig2;
            end

            %% Compute solution

            t1 = sig1_upd * V1 * V1 + sig2_upd * V2 * V2;
            t2 = sig1_upd * V1 * V3 + sig2_upd * V2 * V4;
            t3 = sig1_upd * V3 * V3 + sig2_upd * V4 * V4;

            y(1,:,ii) = valuex * t1 + valuey * t2;
            y(2,:,ii) = valuex * t2 + valuey * t3;

        end
    end

    if flag_perm==1, y=ipermute(y,[2,1,3]); end
    y=reshape(y,sizex);

end


function [u] = l1projection(u)

    fLarge = 100000000;

    MySum = fLarge;
    shrinkfactor = 0;

    while(MySum > 1)

        u=max(u-shrinkfactor,0);
        MySum=sum(abs(u));
        num=sum(u~=0);

        if(num > 0)
            shrinkfactor = (MySum - 1) / num;
        else
            break;
        end
    end

end
