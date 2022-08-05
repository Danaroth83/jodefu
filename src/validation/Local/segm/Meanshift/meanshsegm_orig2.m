function l=meanshsegm(f,hs,hr) ;
% mean shift segmentation of image f, hs is the kernel size in space,
% hr in range (intensity). l are the output labels
% Use Epanechnikov kernel in the joint domain  
  
[ny,nx,nc]=size(f) ;

h=[hs hs repmat(hr,1,nc)] ;

z=zeros(ny,nx,nc+2) ;
for ix=1:nx,
  for  iy=1:ny,
    % initial point for the mean-shift
    y=double([ix iy reshape(f(iy,ix,:),1,nc) ]) ;
    while true,
      % find bounds
      xl=max(ix-hs,1) ;  xh=min(ix+hs,nx) ; nxw=xh-xl+1 ;
      yl=max(iy-hs,1) ;  yh=min(iy+hs,ny) ; nyw=yh-yl+1 ;
      nw=nxw*nyw ; iw=(0:(nw-1))' ;
      % window
      fw=[ fix(iw/nyw+xl) mod(iw,nyw)+yl reshape(f(yl:yh,xl:xh,:),[],nc) ] ;
      fw=double(fw) ;
      % squared distance
      r=(fw-repmat(y,nw,1))./repmat(h,nw,1) ; 
      r=sum(r.*r,2) ;
      % find indices of points smaller than 1
      ind=(r<1.0) ;
      % new point
      y0=y ;
      y=mean(fw(ind,:),1) ;
      if norm(y-y0)<1e-3,
        break ;
      end ;
    end ; % while loop
    z(iy,ix,:)=y ;
    % reshape(z(iy,ix,:),1,nc+2) 
  end ;
end ;


% Region merging using a supergrid structure

s=ones(2*ny+1,2*nx+1,'int8') ;
s(1:2:(2*ny+1),:)=zeros(ny+1,2*nx+1,'int8') ;
s(:,1:2:(2*nx+1))=zeros(2*ny+1,nx+1,'int8') ;
s(2:2:2*ny,3:2:(2*nx-1))=all(abs(z(:,2:end,:)-z(:,1:(end-1),:))<...
                             repmat(reshape(h,1,1,[]),[ny nx-1 1]),3) ;
s(3:2:(2*ny-1),2:2:2*nx)=all(abs(z(2:end,:,:)-z(1:(end-1),:,:))<...
                             repmat(reshape(h,1,1,[]),[ny-1 nx 1]),3) ;
l=bwlabel(s,4) ; % find connected regions
l=l(2:2:2*ny,2:2:2*nx) ;


  
  