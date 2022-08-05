function l=meanshsegm(im,hs,hr) ;
% MEANSHSEGM Mean shift segmentation
% CMP Vision Algorithms http://visionbook.felk.cvut.cz
%  
% Segment a grayscale or color image using mean shift image segmentation
% Warning: Segmenting a medium-size image (200x 300 pixels)
% takes a couple of minutes and the time increases with
% the image and kernel size.  
%
% Usage: l = meanshsegm(im,hs,hr)
%  im   [m x n x d]  Scalar (d=1) or color (d=3) input image.
%  hs  (default 10)
%      Spatial kernel size, in pixels.
%  hr  (default 20)
%      Range kernel size.
% l  [m x n]  Output labeling. Each pixel position contains an integer
%  1... N corresponding to an assigned region number; N is the
%  number of regions.
%   keyboard
% Set default parameters, if needed
if nargin<3,
  hr=20 ;
end ;

if nargin<3,
  hs=10 ;
end ;

% Prepare the joint space-range scaling vector 
% h
% Matrix z will store the filtered pixels z after the mean shift
% discontinuity preserving filtering .
[m,n,d] = size(im);
h  = [hs hs repmat(hr,1,d)];
z  = zeros( m, n, d+2 );
im = double(im);

% For each pixel, set up an initial value y concatenating spatial
% coordinates and the pixel value (scalar or vector).
for ix = 1:n
  for iy = 1:m
    y = double( [ix iy reshape(im(iy,ix,:),1,d)] );
% Perform mean-shift filtering until convergence. We only need to consider
% a spatial neighborhood of 2 h_s x 2 h_s pixels, since pixels
% further away are guaranteed to fall outside the kernel.  The contents of
% this neighborhood (coordinates and pixel values) are rearranged into
% a matrix fw. Note how the x and y coordinates are generated using
% integer division and the modulo function. The same operation could be
% accomplished using ndgrid which would be more elegant but
% unfortunately much slower.
    xl = max(ix-hs,1);  xh = min(ix+hs,n);  nw = xh-xl+1;
    yl = max(iy-hs,1);  yh = min(iy+hs,m);  mw = yh-yl+1;
    nw = nw*mw;  iw = (0:(nw-1))';
    fw = [fix(iw/mw+xl) mod(iw,mw)+yl reshape(im(yl:yh,xl:xh,:),[],d)];
% We are now ready to perform the mean shift
% steps . We calculate the  
% squared scaled distance r^2 with respect to the current point y.
% Indices of all points that are sufficiently close to fall into the
% support of the kernel (r<1) are stored to ind.
% Since the derivative profile g_(r) is constant for r<1,
% the new position y is simply a mean of the points involved.
% If convergence is detected, we exit the while-loop and store
% the result to z, otherwise we continue to iterate.
    while true
      r = (fw-repmat(y,nw,1)) ./ repmat(h,nw,1);
      r = sum( r.*r, 2 );
      ind = r<1.0;
      y0 = y;
      y = mean( fw(ind,:), 1 );
      if norm(y-y0)<1e-3*max(hs,hr), break; end
    end
    z(iy,ix,:) = y;
  end % for iy
end % for ix

% The second part of mean shift segmentation is clustering of the basins
% of attraction:
% we use the supergrid technique 
% Supergrid elements corresponding
% to pixels are all set to one. Supergrid elements corresponding to
% an edge between two pixels are set to one if their filtered
% values in z are closer than h_s in the spatial domain
% and closer than h_r in the range domain. Function bwlabel
% is used to find the connected components that define the
% final labelling.
s = ones( 2*m+1, 2*n+1, 'int8' );
s(1:2:(2*m+1),:) = zeros( m+1, 2*n+1, 'int8' );
s(:,1:2:(2*n+1)) = zeros( 2*m+1, n+1, 'int8' );
s(2:2:2*m,3:2:(2*n-1)) = all(cat(3, ...                % horizontal edges
  abs(z(:,2:end,1:2)-z(:,1:(end-1),1:2)) < hs, ...
  abs(z(:,2:end,3:end)-z(:,1:(end-1),3:end)) < hr ),3);
s(3:2:(2*m-1),2:2:2*n) = all(cat(3, ...                % vertical edges
  abs(z(2:end,:,1:2)-z(1:(end-1),:,1:2)) < hs, ...
  abs(z(2:end,:,3:end)-z(1:(end-1),:,3:end)) < hr ),3);
l = bwlabel( s, 4 );                                   % find connected regions
l = l( 2:2:2*m, 2:2:2*n );                             % extract labeling

