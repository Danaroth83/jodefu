function outfn = MSG_Lukac(inmosin)

inmosin = extendImage(inmosin,20);
[M,N] = size(inmosin);

N = N+1;
inmos = zeros(M,N);
inmos(:,2:N) = inmosin(:,:);
g = inmos;

out = zeros(M,N,3);
out(:,:,2) = inmos(:,:);

for a=1:M
    for b=1:N
        if (mod(a,4)==1 && mod(b,4)==1)
            out(a,b,1) = inmos(a,b);
        elseif (mod(a,4)==1 && mod(b,4)==3)
            out(a,b,1) = inmos(a,b);
        elseif (mod(a,4)==2 && mod(b,4)==1)
            out(a,b,3) = inmos(a,b);
        elseif (mod(a,4)==2 && mod(b,4)==3)
            out(a,b,3) = inmos(a,b);
        elseif (mod(a,4)==3 && mod(b,4)==2)
            out(a,b,1) = inmos(a,b);
        elseif (mod(a,4)==3 && mod(b,4)==0)
            out(a,b,1) = inmos(a,b);
        elseif (mod(a,4)==0 && mod(b,4)==2)
            out(a,b,3) = inmos(a,b);
        elseif (mod(a,4)==0 && mod(b,4)==0)
            out(a,b,3) = inmos(a,b);
        end
    end
end

h = zeros(M,N,3);
h(:,:,:) = out(:,:,:);

fh=[-1/4 1/2 1/2 1/2 -1/4];
fv2=[-1/8 1/4 0 0 -1/8; 0 0 1/2 1/2 0; -1/8 1/4 0 0 -1/8];
fv1=[-1/8 0 0 1/4 -1/8; 0 1/2 1/2 0 0; -1/8 0 0 1/4 -1/8];

Ah=conv2(g,fh);
Ah=Ah(:,3:2+N);

Av=zeros(M,N);
for i=3:2:M-3
    for j=2:N-1
        Av(i,j) = sum(sum(g(i-2:i+2,j-1:j+1) .* fv1'));
        Av(i+1,j) = sum(sum(g(i-1:i+3,j-1:j+1) .* fv2'));
    end
end

dh=zeros(M,N);
dv=dh;
for i=2:4:M
   dh(i,2:2:N)=g(i,2:2:N)-Ah(i,2:2:N);
   dh(i,1:2:N)=Ah(i,1:2:N)-g(i,1:2:N);
   dv(i,2:2:N)=g(i-1,2:2:N)-Av(i-1,2:2:N);
   dv(i,1:2:N)=Av(i,1:2:N)-g(i,1:2:N);
end
for i=1:4:M
   dh(i,2:2:N)=g(i,2:2:N)-Ah(i,2:2:N);
   dh(i,1:2:N)=Ah(i,1:2:N)-g(i,1:2:N);
   dv(i,2:2:N)=g(i+1,2:2:N)-Av(i+1,2:2:N);
   dv(i,1:2:N)=Av(i,1:2:N)-g(i,1:2:N);
end
for i=4:4:M
   dh(i,1:2:N)=g(i,1:2:N)-Ah(i,1:2:N);
   dh(i,2:2:N)=Ah(i,2:2:N)-g(i,2:2:N);
   dv(i,1:2:N)=g(i-1,1:2:N)-Av(i-1,1:2:N);
   dv(i,2:2:N)=Av(i,2:2:N)-g(i,2:2:N);
end
for i=3:4:M
   dh(i,1:2:N)=g(i,1:2:N)-Ah(i,1:2:N);
   dh(i,2:2:N)=Ah(i,2:2:N)-g(i,2:2:N);
   dv(i,1:2:N)=g(i+1,1:2:N)-Av(i+1,1:2:N);
   dv(i,2:2:N)=Av(i,2:2:N)-g(i,2:2:N);
end

vmap = zeros(M,N);
hmap = zeros(M,N);

for i=10:1:M-10
    for j=10:1:N-10
        vmap(i,j) = 1.75*(abs((g(i-2,j)-g(i+2,j))/4 +(g(i+4,j)-g(i-4,j))/10+(g(i-6,j)-g(i+6,j))/24));
        hmap(i,j) = abs((g(i,j-1)-g(i,j+1))/2+(g(i,j+2)-g(i,j-2))/5+(g(i,j-3)-g(i,j+3))/12);
    end
end

% 1.Interpolate green pixels.
gree=zeros(M,N);
gree=double(gree);
gree(:,:) = h(:,:,2);

vmapu = zeros(M,N);
hmapl = zeros(M,N);
vmapd = zeros(M,N);
hmapr = zeros(M,N);

vmapv = zeros(M,N);
hmaph = zeros(M,N);
   
for i=5:1:M-5
    for j=5:1:N-5  
        vmapu(i,j) = 100/((sum(vmap(i-3:i,j))+0.4*sum(vmap(i-3:i,j-1))+0.4*sum(vmap(i-3:i,j+1)))^2+0.01);
        vmapd(i,j) = 100/((sum(vmap(i:i+3,j))+0.4*sum(vmap(i:i+3,j-1))+0.4*sum(vmap(i:i+3,j+1)))^2+0.01);
        hmapl(i,j) = 100/((sum(hmap(i,j-3:j))+0.4*sum(hmap(i-1,j-3:j))+0.4*sum(hmap(i+1,j-3:j)))^2+0.01);
        hmapr(i,j) = 100/((sum(hmap(i,j:j+3))+0.4*sum(hmap(i-1,j:j+3))+0.4*sum(hmap(i+1,j:j+3)))^2+0.01);
        
        vmapv(i,j) = 100/(sum(sum(vmap(i-2:i+2,j-2:j+2)))^3+0.01);
        hmaph(i,j) = 100/(sum(sum(hmap(i-2:i+2,j-2:j+2)))^3+0.01);
    end
end

for i=4:4:M-6
    for j=5:2:N-6
        h(i+1,j,2) = ((dv(i+1,j)+dv(i-1,j)/2+dv(i+3,j)/2)*vmapv(i+1,j) + (dh(i+1,j)+dh(i+1,j-1)/2+dh(i+1,j+1)/2)*hmaph(i+1,j))/(2*(vmapv(i+1,j) + hmaph(i+1,j))) + g(i+1,j);
    end
end

for i=6:4:M-6
    for j=4:2:N-6
        h(i,j+1,2) = ((dv(i,j+1)+dv(i-2,j+1)/2+dv(i+2,j+1)/2)*vmapv(i,j+1) + (dh(i,j+1)+dh(i,j)/2+dh(i,j+2)/2)*hmaph(i,j+1))/(2*(vmapv(i,j+1) + hmaph(i,j+1))) + g(i,j+1);
    end
end

for i=6:4:M-6
    for j=6:2:N-6
        h(i+1,j,2) = ((dv(i+1,j)+dv(i-1,j)/2+dv(i+3,j)/2)*vmapv(i+1,j) + (dh(i+1,j)+dh(i+1,j-1)/2+dh(i+1,j+1)/2)*hmaph(i+1,j))/(2*(vmapv(i+1,j) + hmaph(i+1,j))) + g(i+1,j);
    end
end

for i=8:4:M-6
    for j=5:2:N-6
        h(i,j+1,2) = ((dv(i,j+1)+dv(i-2,j+1)/2+dv(i+2,j+1)/2)*vmapv(i,j+1) + (dh(i,j+1)+dh(i,j)/2+dh(i,j+2)/2)*hmaph(i,j+1))/(2*(vmapv(i,j+1) + hmaph(i,j+1))) + g(i,j+1);
    end
end

% 2. Update green pixels.
for i=4:4:M-6
    for j=5:2:N-6
        vtot = vmapu(i+1,j)+vmapd(i+1,j)+hmapl(i+1,j)+hmapr(i+1,j);
        avv = (h(i+1,j,2) - g(i+1,j))*0.3 + (vmapu(i+1,j)*(h(i-3,j,2) - g(i-3,j) + h(i-1,j-1,2) - g(i-1,j-1) + h(i-1,j+1,2) - g(i-1,j+1))/3 + vmapd(i+1,j)*(h(i+5,j,2) - g(i+5,j) + h(i+3,j-1,2) - g(i+3,j-1) + h(i+3,j+1,2) - g(i+3,j+1))/3 + hmapl(i+1,j)*(h(i+1,j-2,2) - g(i+1,j-2)) + hmapr(i+1,j)*(h(i+1,j+2,2) - g(i+1,j+2)))*0.7/vtot;
        gree(i+1,j) = g(i+1,j) + avv;
    end
end

for i=6:4:M-6
    for j=4:2:N-6
        vtot = vmapu(i,j+1)+vmapd(i,j+1)+hmapl(i,j+1)+hmapr(i,j+1);
        avv = (h(i,j+1,2) - g(i,j+1))*0.3 + (vmapu(i,j+1)*(h(i-4,j+1,2) - g(i-4,j+1) + h(i-2,j,2) - g(i-2,j) + h(i-2,j+2,2) - g(i-2,j+2))/3 + vmapd(i,j+1)*(h(i+4,j+1,2) - g(i+4,j+1) + h(i+2,j,2) - g(i+2,j) + h(i+2,j+2,2) - g(i+2,j+2))/3 + hmapl(i,j+1)*(h(i,j-1,2) - g(i,j-1)) + hmapr(i,j+1)*(h(i,j+3,2) - g(i,j+3)))*0.7/vtot;
        gree(i,j+1) = g(i,j+1) + avv;
    end
end

for i=6:4:M-6
    for j=4:2:N-6
        vtot = vmapu(i+1,j)+vmapd(i+1,j)+hmapl(i+1,j)+hmapr(i+1,j);
        avv = (h(i+1,j,2) - g(i+1,j))*0.3 + (vmapu(i+1,j)*(h(i-3,j,2) - g(i-3,j) + h(i-1,j-1,2) - g(i-1,j-1) + h(i-1,j+1,2) - g(i-1,j+1))/3 + vmapd(i+1,j)*(h(i+5,j,2) - g(i+5,j) + h(i+3,j-1,2) - g(i+3,j-1) + h(i+3,j+1,2) - g(i+3,j+1))/3 + hmapl(i+1,j)*(h(i+1,j-2,2) - g(i+1,j-2)) + hmapr(i+1,j)*(h(i+1,j+2,2) - g(i+1,j+2)))*0.7/vtot;
        gree(i+1,j) = g(i+1,j) + avv;
    end
end

for i=8:4:M-6
    for j=5:2:N-6
        vtot = vmapu(i,j+1)+vmapd(i,j+1)+hmapl(i,j+1)+hmapr(i,j+1);
        avv = (h(i,j+1,2) - g(i,j+1))*0.3 + (vmapu(i,j+1)*(h(i-4,j+1,2) - g(i-4,j+1) + h(i-2,j,2) - g(i-2,j) + h(i-2,j+2,2) - g(i-2,j+2))/3 + vmapd(i,j+1)*(h(i+4,j+1,2) - g(i+4,j+1) + h(i+2,j,2) - g(i+2,j) + h(i+2,j+2,2) - g(i+2,j+2))/3 + hmapl(i,j+1)*(h(i,j-1,2) - g(i,j-1)) + hmapr(i,j+1)*(h(i,j+3,2) - g(i,j+3)))*0.7/vtot;
        gree(i,j+1) = g(i,j+1) + avv;
    end
end

h(:,:,2) = gree(:,:);

out(:,:,:) = h(:,:,:);

% 3. Red & Blue interpolation
% red  interpolation
for i=5:4:M-6
    for j=5:2:N-6
        out(i,j+1,1) = out(i,j+1,2) - ((out(i,j,2) - out(i,j,1))*4 + (out(i,j+2,2) - out(i,j+2,1))*4 + (out(i-2,j+1,2) - out(i-2,j+1,1)) + (out(i+2,j+1,2) - out(i+2,j+1,1)))/10;
        out(i+1,j,1) = out(i+1,j,2) - ((out(i,j,2) - out(i,j,1))*2 + (out(i+2,j-1,2) - out(i+2,j-1,1)) + (out(i+2,j+1,2) - out(i+2,j+1,1)))/4;
        out(i+1,j+1,1) = out(i+1,j+1,2) - ((out(i,j,2) - out(i,j,1)) + (out(i,j+2,2) - out(i,j+2,1)) + (out(i+2,j+1,2) - out(i+2,j+1,1))*2)/4;
    end
end

for i=7:4:M-6
    for j=6:2:N-6
        out(i,j+1,1) = out(i,j+1,2) - ((out(i,j,2) - out(i,j,1))*4 + (out(i,j+2,2) - out(i,j+2,1))*4 + (out(i-2,j+1,2) - out(i-2,j+1,1)) + (out(i+2,j+1,2) - out(i+2,j+1,1)))/10;
        out(i+1,j,1) = out(i+1,j,2) - ((out(i,j,2) - out(i,j,1))*2 + (out(i+2,j-1,2) - out(i+2,j-1,1)) + (out(i+2,j+1,2) - out(i+2,j+1,1)))/4;
        out(i+1,j+1,1) = out(i+1,j+1,2) - ((out(i,j,2) - out(i,j,1)) + (out(i,j+2,2) - out(i,j+2,1)) + (out(i+2,j+1,2) - out(i+2,j+1,1))*2)/4;
    end
end

% blue interpolation
for i=6:4:M-6
    for j=5:2:N-6
        out(i,j+1,3) = out(i,j+1,2) - ((out(i,j,2) - out(i,j,3))*4 + (out(i,j+2,2) - out(i,j+2,3))*4 + (out(i-2,j+1,2) - out(i-2,j+1,3)) + (out(i+2,j+1,2) - out(i+2,j+1,3)))/10;
        out(i+1,j,3) = out(i+1,j,2) - ((out(i,j,2) - out(i,j,3))*2 + (out(i+2,j-1,2) - out(i+2,j-1,3)) + (out(i+2,j+1,2) - out(i+2,j+1,3)))/4;
        out(i+1,j+1,3) = out(i+1,j+1,2) - ((out(i,j,2) - out(i,j,3)) + (out(i,j+2,2) - out(i,j+2,3)) + (out(i+2,j+1,2) - out(i+2,j+1,3))*2)/4;
    end
end

for i=8:4:M-6
    for j=6:2:N-6
        out(i,j+1,3) = out(i,j+1,2) - ((out(i,j,2) - out(i,j,3))*4 + (out(i,j+2,2) - out(i,j+2,3))*4 + (out(i-2,j+1,2) - out(i-2,j+1,3)) + (out(i+2,j+1,2) - out(i+2,j+1,3)))/10;
        out(i+1,j,3) = out(i+1,j,2) - ((out(i,j,2) - out(i,j,3))*2 + (out(i+2,j-1,2) - out(i+2,j-1,3)) + (out(i+2,j+1,2) - out(i+2,j+1,3)))/4;
        out(i+1,j+1,3) = out(i+1,j+1,2) - ((out(i,j,2) - out(i,j,3)) + (out(i,j+2,2) - out(i,j+2,3)) + (out(i+2,j+1,2) - out(i+2,j+1,3))*2)/4;
    end
end

out=floor(out+0.5);
for i=1:M
    for j=1:N
        for k=1:3  
            if (out(i,j,k)>2^8-1)
                out(i,j,k)=2^8-1;
            elseif (out(i,j,k)<0)
                out(i,j,k)=0;
            end
        end
    end
end

outfn = out(21:M-20,22:N-20,:);
outfn = uint8(outfn);

return;