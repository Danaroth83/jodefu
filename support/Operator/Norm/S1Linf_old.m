function y=S1Linf_new(x,gamma)

[N1,N2,N3]=size(x);
y=zeros([N1,N2,N3]);
len_diag_Sigma=min(N1,N2);

if N1>2

    for ii=1:N3
        A=x(:,:,ii);
        [~,Sigma,V]=svd(A,'econ');
        diag_Sigma=diag(Sigma);

        tmp=abs(diag_Sigma);
        % prox_Sigma = diag_Sigma - max(bsxfun(@minus,tmp,max(max(bsxfun(@rdivide,cumsum(sort(...
        % tmp,1,'descend'),1)-gamma,(1:length(diag_Sigma))'),[],1),0)),0).*sign(diag_Sigma);
        prox_Sigma = diag_Sigma - max(tmp-max(max((cumsum(sort(tmp,1,'descend'),...
            1)-gamma)./(1:len_diag_Sigma)',[],1),0),0).*sign(diag_Sigma);

        div_Sigma=prox_Sigma./diag_Sigma;
        div_Sigma(diag_Sigma==0)=0;

        y(:,:,ii)=A*V*diag(div_Sigma)*V';
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

end