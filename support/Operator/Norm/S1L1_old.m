function y=S1L1_new(x,gamma)

[N1,N2,N3]=size(x);
y=zeros([N1,N2,N3]);

if N1==2

    for ii=1:N3

        
        A=x(:,:,ii);
        [~,Sigma,V]=svd(A,'econ');
        diag_Sigma=diag(Sigma);

        % prox_Sigma = diag_Sigma - bsxfun(@rdivide, diag_Sigma,max(abs(diag_Sigma)/gamma,1));
        prox_Sigma = diag_Sigma - diag_Sigma./max(abs(diag_Sigma)/gamma,1);
        
        div_Sigma=prox_Sigma./diag_Sigma;
        div_Sigma(isnan(div_Sigma))=0;
        
        y(:,:,ii)=A*V*diag(div_Sigma)*V';

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

    %%        // Compute prox_p of the diagonal entries


    %         if(p == 1)

            sig1_upd = max(sig1 - gamma, 0);
            sig2_upd = max(sig2 - gamma, 0);
    %             
    %         elseif(p == INFNORM)
    %         
    %             proj[0] = fabs(sig1)/gamma;
    %             proj[1] = fabs(sig2)/gamma;
    %             l1projection(proj, 2);
    %             
    %             sig1_upd = sig1 - gamma * proj[0];
    %             sig2_upd = sig2 - gamma * proj[1];
    %         end

    %%        // Compute the diagonal entries of $\widehat{\Sigma}\Sigma^{\dagger}_0$
            if(sig1 > fTiny)
                sig1_upd = sig1_upd/sig1;
            end
            if(sig2 > fTiny)
                sig2_upd = sig2_upd/sig2;
            end
            %% // Compute solution

            t1 = sig1_upd * V1^2 + sig2_upd * V2^2;
            t2 = sig1_upd * V1 * V3 + sig2_upd * V2 * V4;
            t3 = sig1_upd * V3^2 + sig2_upd * V4^2;

            y(1,:,ii)= valuex * t1 + valuey * t2;
            y(2,:,ii)= valuex * t2 + valuey * t3;
    end
    
end