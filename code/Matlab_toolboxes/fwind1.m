function h = fwind1(f1,f2,hd,w1,w2)

method1 = 'linear';
method2 = 'bilinear';

if nargin==2 % Uniform spacing case with Huang's method
  hd = f1; w1 = f2;
  n = length(w1); m = n;
elseif nargin==3 % Uniform spacing with separable window
  w2 = hd;  hd = f1; w1 = f2;
  n = length(w1); m = length(w2);
elseif nargin==4 % Non-uniform spacing with Huang's method
  n = length(w1); m = n;
else
  n = length(w1); m = length(w2);
end


%
% Create 2-D window using either Huang's method or separable method.
%
if nargin==2 || nargin==4 % Huang's method: Create 2-D circular window
  if any(abs(w1-rot90(w1,2))>sqrt(eps))
    error('1d window must be simmetric')
  end
  if length(w1)<2
    error('Window length must be greater than one')
  end

  t = (-(n-1)/2:(n-1)/2)*(2/(n-1));
  [t1,t2] = meshgrid(t,t);
  t12 = sqrt(t1.*t1 + t2.*t2);
  d = find(t12<t(1) | t12>t(length(w1)));
  if ~isempty(d), t12(d) = zeros(size(d)); end
  w = zeros(size(t12)); w(:) = interp1(t,w1,t12(:),method1);
  if ~isempty(d), w(d) = zeros(size(d)); end

else % Create separable window
  w = w2(:)*w1(:).';

end

%
% Design filter using fsamp2 and apply window
%
if nargin<4, % Uniformly spaced data
  % Interpolate Hd to be the same size as W, if necessary
  if any([m n]~=size(hd)), 
    if any(size(hd)<[2 2]),
        error(message('images:fwind1:hdMustHaveAtLeast2rowsAnd2cols'))
    end

    [f1,f2] = freqspace(size(hd));
    % Extrapolate hd so that interpolation is never out of range.
    [mh,nh] = size(hd);
    if floor(nh/2)==nh/2 % if even
      hd = [hd,hd(:,1)]; f1 = [f1 1];
    else
      hd = [zeros(mh,1) hd zeros(mh,1)]; 
      df = f1(2)-f1(1); f1 = [f1(1)-df f1 f1(nh)+df];
    end
    [mh,nh] = size(hd);
    if floor(mh/2)==mh/2 % if even
      hd = [hd;hd(1,:)]; f2 = [f2 1];
    else
      hd = [zeros(1,nh);hd;zeros(1,nh)]; 
      df = f2(2)-f2(1); f2 = [f2(1)-df f2 f2(mh)+df];
    end
    [t1,t2] = freqspace([m n],'meshgrid');
    
    % Promote to double for call to interp2
    if ~isa(hd,'double')
       hd = double(hd);
    end
    
    hd = interp2(f1,f2,hd,t1,t2,method2);
    d = find(isnan(hd)); if ~isempty(d), hd(d) = zeros(size(d)); end
  end
  h = fsamp2(hd) .* w;

else % Non-uniformly spaced data
  h = fsamp2(f1,f2,hd,size(w)) .* w;

end