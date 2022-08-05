function [hh] = MSG_Bayer(A)

[M,N,ch] = size(A);

A = extendImage(A,10);
M = M+20;
N = N+20;
g = zeros(M,N);

% Initialize output image.
h=zeros(M,N,3);

if ch==3
    h(1:2:M,1:2:N,2)=A(1:2:M,1:2:N,2);
    h(2:2:M,2:2:N,2)=A(2:2:M,2:2:N,2);
    h(1:2:M,2:2:N,1)=A(1:2:M,2:2:N,1);
    h(2:2:M,1:2:N,3)=A(2:2:M,1:2:N,3);
    g(1:2:M,1:2:N)=A(1:2:M,1:2:N,2);
    g(2:2:M,2:2:N)=A(2:2:M,2:2:N,2);
    g(1:2:M,2:2:N)=A(1:2:M,2:2:N,1);
    g(2:2:M,1:2:N)=A(2:2:M,1:2:N,3);
else
    g = A;
    h(1:2:M,1:2:N,2)=A(1:2:M,1:2:N);
    h(2:2:M,2:2:N,2)=A(2:2:M,2:2:N);
    h(1:2:M,2:2:N,1)=A(1:2:M,2:2:N);
    h(2:2:M,1:2:N,3)=A(2:2:M,1:2:N);
end

f=[-1/4 1/2 1/2 1/2 -1/4];
Ah=conv2(g,f);
Ah=Ah(:,3:2+N);
Av=conv2(g,f');
Av=Av(3:2+M,:);

dh=zeros(M,N);
dv=dh;
for i=1:2:M
   dh(i,1:2:N)=g(i,1:2:N)-Ah(i,1:2:N);
   dh(i,2:2:N)=Ah(i,2:2:N)-g(i,2:2:N);
   dv(i,1:2:N)=g(i,1:2:N)-Av(i,1:2:N);
   dv(i,2:2:N)=Av(i,2:2:N)-g(i,2:2:N);
end
for i=2:2:M
   dh(i,2:2:N)=g(i,2:2:N)-Ah(i,2:2:N);
   dh(i,1:2:N)=Ah(i,1:2:N)-g(i,1:2:N);
   dv(i,2:2:N)=g(i,2:2:N)-Av(i,2:2:N);
   dv(i,1:2:N)=Av(i,1:2:N)-g(i,1:2:N);
end

vmap = zeros(M,N);
hmap = zeros(M,N);

for i=8:1:M-8
    for j=8:1:N-8
        vmap(i,j) = abs((g(i-1,j)-g(i+1,j))/2+(g(i+2,j)-g(i-2,j))/3.5+(g(i-3,j)-g(i+3,j))/7+(g(i+4,j)-g(i-4,j))/14+(g(i-5,j)-g(i+5,j))/28);
        hmap(i,j) = abs((g(i,j-1)-g(i,j+1))/2+(g(i,j+2)-g(i,j-2))/3.5+(g(i,j-3)-g(i,j+3))/7+(g(i,j+4)-g(i,j-4))/14+(g(i,j-5)-g(i,j+5))/28);      
    end
end

% 1. Interpolate green pixels.

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

for i=6:2:M-6
    for j=6:2:N-6
        h(i+1,j,2) = ((dv(i+1,j)+dv(i,j)/2+dv(i+2,j)/2)*vmapv(i+1,j) + (dh(i+1,j)+dh(i+1,j-1)/2+dh(i+1,j+1)/2)*hmaph(i+1,j))/(2*(vmapv(i+1,j) + hmaph(i+1,j))) + g(i+1,j);
    end
end

for i=6:2:M-6
    for j=6:2:N-6
        h(i,j+1,2) = ((dv(i,j+1)+dv(i-1,j+1)/2+dv(i+1,j+1)/2)*vmapv(i,j+1) + (dh(i,j+1)+dh(i,j)/2+dh(i,j+2)/2)*hmaph(i,j+1))/(2*(vmapv(i,j+1) + hmaph(i,j+1))) + g(i,j+1);
    end
end

% 2. Update green pixels.
for i=6:2:M-6
    for j=6:2:N-6
        vtot = vmapu(i+1,j)+vmapd(i+1,j)+hmapl(i+1,j)+hmapr(i+1,j);
        avv = (h(i+1,j,2) - g(i+1,j))*0.3 + (vmapu(i+1,j)*(h(i-1,j,2) - g(i-1,j)) + vmapd(i+1,j)*(h(i+3,j,2) - g(i+3,j)) + hmapl(i+1,j)*(h(i+1,j-2,2) - g(i+1,j-2)) + hmapr(i+1,j)*(h(i+1,j+2,2) - g(i+1,j+2)))*0.7/vtot;
        gree(i+1,j) = g(i+1,j) + avv;
    end
end

for i=6:2:M-6
    for j=6:2:N-6
        vtot = vmapu(i,j+1)+vmapd(i,j+1)+hmapl(i,j+1)+hmapr(i,j+1);
        avv = (h(i,j+1,2) - g(i,j+1))*0.3 + (vmapu(i,j+1)*(h(i-2,j+1,2) - g(i-2,j+1)) + vmapd(i,j+1)*(h(i+2,j+1,2) - g(i+2,j+1)) + hmapl(i,j+1)*(h(i,j-1,2) - g(i,j-1)) + hmapr(i,j+1)*(h(i,j+3,2) - g(i,j+3)))*0.7/vtot;
        gree(i,j+1) = g(i,j+1) + avv;
    end
end

h(:,:,2) = gree(:,:);


% 3. Red & Blue interpolation
% Find the red component of blue pixels.
for i=4:2:M-4
    for j=4:2:N-4
        h(i,j+1,1) = h(i,j+1,2) - ((h(i-1,j,2)-h(i-1,j,1))+(h(i-1,j+2,2)-h(i-1,j+2,1))+(h(i+1,j,2)-h(i+1,j,1))+(h(i+1,j+2,2)-h(i+1,j+2,1)))*1.25/4 ...
            + ((h(i-3,j,2)-h(i-3,j,1))+(h(i-1,j-2,2)-h(i-1,j-2,1))+(h(i-3,j+2,2)-h(i-3,j+2,1))+(h(i-1,j+4,2)-h(i-1,j+4,1))+(h(i+3,j,2)-h(i+3,j,1))+(h(i+1,j-2,2)-h(i+1,j-2,1))+(h(i+3,j+2,2)-h(i+3,j+2,1))+(h(i+1,j+4,2)-h(i+1,j+4,1)))*0.25/8;
    end
end

% Find the blue component of red pixels.
for i=4:2:M-4
    for j=4:2:N-4
        h(i+1,j,3) = h(i+1,j,2) - ((h(i,j-1,2)-h(i,j-1,3))+(h(i,j+1,2)-h(i,j+1,3))+(h(i+2,j-1,2)-h(i+2,j-1,3))+(h(i+2,j+1,2)-h(i+2,j+1,3)))*1.25/4 ...
            + ((h(i-2,j-1,2)-h(i-2,j-1,3))+(h(i,j-3,2)-h(i,j-3,3))+(h(i-2,j+1,2)-h(i-2,j+1,3))+(h(i,j+3,2)-h(i,j+3,3))+(h(i+4,j-1,2)-h(i+4,j-1,3))+(h(i+2,j-3,2)-h(i+2,j-3,3))+(h(i+4,j+1,2)-h(i+4,j+1,3))+(h(i+2,j+3,2)-h(i+2,j+3,3)))*0.25/8;
    end
end

% Find red and blue components of green pixels.
for i=6:2:M-6
    for j=6:2:N-6
        h(i,j,3) = h(i,j,2) - (((h(i,j-1,2)-h(i,j-1,3)) + (h(i,j+1,2)-h(i,j+1,3)))*hmaph(i,j) + ((h(i-1,j,2)-h(i-1,j,3)) + (h(i+1,j,2)-h(i+1,j,3)))*vmapv(i,j))/(2*(vmapv(i,j) + hmaph(i,j)));
        h(i,j,1) = h(i,j,2) - (((h(i,j-1,2)-h(i,j-1,1)) + (h(i,j+1,2)-h(i,j+1,1)))*hmaph(i,j) + ((h(i-1,j,2)-h(i-1,j,1)) + (h(i+1,j,2)-h(i+1,j,1)))*vmapv(i,j))/(2*(vmapv(i,j) + hmaph(i,j)));
        h(i+1,j+1,1) = h(i+1,j+1,2) - (((h(i,j+1,2)-h(i,j+1,1)) + (h(i+2,j+1,2)-h(i+2,j+1,1)))*vmapv(i+1,j+1) + ((h(i+1,j,2)-h(i+1,j,1)) + (h(i+1,j+2,2)-h(i+1,j+2,1)))*hmaph(i+1,j+1))/(2*(vmapv(i+1,j+1) + hmaph(i+1,j+1)));
        h(i+1,j+1,3) = h(i+1,j+1,2) - (((h(i,j+1,2)-h(i,j+1,3)) + (h(i+2,j+1,2)-h(i+2,j+1,3)))*vmapv(i+1,j+1) + ((h(i+1,j,2)-h(i+1,j,3)) + (h(i+1,j+2,2)-h(i+1,j+2,3)))*hmaph(i+1,j+1))/(2*(vmapv(i+1,j+1) + hmaph(i+1,j+1)));
    end
end

h=floor(h+0.5);
for i=1:M
    for j=1:N
        for k=1:3  
            if (h(i,j,k)>2^8-1)
                h(i,j,k)=2^8-1;
            elseif (h(i,j,k)<0)
                h(i,j,k)=0;
            end
        end
    end
end

hh=h(11:M-10,11:N-10,:);

return;