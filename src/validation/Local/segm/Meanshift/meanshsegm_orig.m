function l=meanshsegm(f,hs,hr) ;
% mean shift segmentation of image f, hs is the kernel size in space,
% hr in range (intensity). l are the output labels
  
[ny,nx,nc]=size(f) ;

% form the joint domain
[xi,yi]=meshgrid(1:nx,1:ny) ;
x=[ xi(:) yi(:) reshape(f,[],nc) ] ;


h=[hs hs repmat(hr,1,nc)] ;

g=meanshiftinit(x,h) ;

z=zeros(size(x)) ;
for i=1:size(x,1),
%for i=1:5,
  z(i,:)=meanshift(g,x(i,:)) ; z(i,:)
end ;  

zl=reshape(z(:,:),ny,nx,nc+2) ;

% Region merging using a supergrid structure

s=ones(2*ny+1,2*nx+1,'int8') ;
s(1:2:(2*ny+1),:)=zeros(ny+1,2*nx+1,'int8') ;
s(:,1:2:(2*nx+1))=zeros(2*ny+1,nx+1,'int8') ;
s(2:2:2*ny,3:2:(2*nx-1))=all(abs(zl(:,2:end,:)-zl(:,1:(end-1),:))<...
                             repmat(reshape(h,1,1,[]),[ny nx-1 1]),3) ;
s(3:2:(2*ny-1),2:2:2*nx)=all(abs(zl(2:end,:,:)-zl(1:(end-1),:,:))<...
                             repmat(reshape(h,1,1,[]),[ny-1 nx 1]),3) ;
l=bwlabel(s,4) ; % find connected regions
l=l(2:2:2*ny,2:2:2*nx) ;


  
  