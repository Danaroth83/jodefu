function im3=extendImage(im,n)

[height, width, depth]=size(im);
im2=zeros(height,width+2*n, depth);

im2(1:height,n+1:width+n,:)=im(:,:,:);

for j=1:depth
    for i=1:n
        im2(:,n-i+1,j)=im(:,i,j);
        im2(:,width+n+i,j)=im(:,width-i,j);
    end
end

im3=zeros(height+2*n,width+2*n, depth);
im3(n+1:height+n,1:width+2*n,:)=im2(:,:,:);

for j=1:depth
    for i=1:n
        im3(n-i+1,:,j)=im2(i,:,j);
        im3(height+n+i,:,j)=im2(height-i,:,j);
    end
end