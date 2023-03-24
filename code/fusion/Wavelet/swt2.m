function varargout = swt2(x,n,varargin)
%SWT2 Discrete stationary wavelet transform 2-D.
%   SWT2 performs a multilevel 2-D stationary wavelet analysis
%   using either a specific orthogonal wavelet ('wname' see 
%   WFILTERS for more information) or specific orthogonal wavelet 
%   decomposition filters.
%
%   SWC = SWT2(X,N,'wname') or [A,H,V,D] = SWT2(X,N,'wname') 
%   compute the stationary wavelet decomposition of the 
%   matrix X at level N, using 'wname'.
%   N must be a strictly positive integer (see WMAXLEV for more
%   information). 2^N must divide size(X,1) and size(X,2).
%
%   Outputs [A,H,V,D] are 3-D Arrays which contain the 
%   coefficients: for 1 <= i <= N,
%   A(:,:,i) contains the coefficients of approximation 
%   of level i. 
%   H(:,:,i), V(:,:,i) and D(:,:,i) contain the coefficients 
%   of details of level i, (Horiz., Vert., Diag.). 
%   SWC = [H(:,:,1:N) ; V(:,:,1:N); D(:,:,1:N); A(:,:,N)].
%
%   SWC = SWT2(X,N,Lo_D,Hi_D) or  [A,H,V,D] = SWT2(X,N,Lo_D,Hi_D)
%   compute the stationary wavelet decomposition as above,
%   given these filters as input: 
%     Lo_D is the decomposition low-pass filter and
%     Hi_D is the decomposition high-pass filter.
%     Lo_D and Hi_D must be the same length.
%
%   NOTE: When X represents an indexed image, then X is an  
%   m-by-n matrix, the output arrays SWC or CA, CH, CV, CD
%   are m-by-n-by-p arrays.
%   When X represents a truecolor image, then it becomes an 
%   m-by-n-by-3 array. This arrays consist of three m-by-n matrices
%  (representing the red, green, and blue color planes) concatenated 
%   along the third dimension. The output arrays SWC or CA, CH, 
%   CV, CD are m-by-n-by-p-by-3 arrays.
%   For more information on image formats, see the reference 
%   pages of IMAGE and IMFINFO functions.
%
%   See also DWT2, WAVEDEC2.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 02-Oct-95.
%   Last Revision 26-Aug-2013.
%   Copyright 1995-2013 The MathWorks, Inc.
%   $Revision: 1.6.4.11 $  $Date: 2013/09/14 19:39:03 $

% Check arguments.
narginchk(3,4)
if errargt(mfilename,n,'int')
    error(message('Wavelet:FunctionArgVal:Invalid_Input'));
end

% Check type of float.
varInfo = whos('x');
if isequal(varInfo.class,'double')
    typeFloat = 'double';
else
    typeFloat = 'single';
end

% Preserve initial size.
s = size(x);
a3d_Flag = length(s)>2;
size2D = s(1:2);
pow = 2^n;
if any(rem(size2D,pow))
    sOK = ceil(size2D/pow)*pow;
    oriStr = ['(' int2str(s(1))   ',' int2str(s(2)) ')'];
    sugStr = ['(' int2str(sOK(1)) ',' int2str(sOK(2)) ')'];
    msg = getWavMSG('Wavelet:moreMSGRF:SWT_size_MSG',n,oriStr,sugStr);  
    errargt(mfilename,msg,'msg');
    varargout = cell(1,nargout);
    return
end        

% Compute decomposition filters.
if ischar(varargin{1})
    [lo,hi] = wfilters(varargin{1},'d');
else
    lo = varargin{1};   hi = varargin{2};
end
lo = lo(:)';
hi = hi(:)';

% Set DWT_Mode to 'per'.
old_modeDWT = dwtmode('status','nodisp');
modeDWT = 'per';
dwtmode(modeDWT,'nodisp');

% Compute non-decimate wavelet coefficients.
a = zeros([s,n],typeFloat);
h = zeros([s,n],typeFloat);
v = zeros([s,n],typeFloat);
d = zeros([s,n],typeFloat);

for k = 1:n
    % Extension.
    lf    = length(lo);
    first = [lf+1,lf+1];
    x  = wextend('2D',modeDWT,x,[lf/2,lf/2]);

    % Decomposition.
    if ~a3d_Flag
        [a(:,:,k),h(:,:,k),v(:,:,k),d(:,:,k)] = decomposeLOC(x);
        x = a(:,:,k);
    else
        for j=1:3
            [a(:,:,j,k),h(:,:,j,k),v(:,:,j,k),d(:,:,j,k)] = ...
                decomposeLOC(x(:,:,j));
        end
        x = a(:,:,:,k);
    end

    % upsample filters.
    lo(2,length(lo)) = 0;
    lo = lo(:)';
    hi(2,length(hi)) = 0;
    hi = hi(:)';    
end

if nargout==4    
    if a3d_Flag && n==1
        varargout = {a,h,v,d};
    else
        varargout = {a,h,v,d};
    end

elseif nargout==1
    if ~a3d_Flag
        varargout{1} = cat(3,h,v,d,a(:,:,n));
    else
        varargout{1} = cat(4,h,v,d,a(:,:,:,n));
    end
end

% Restore DWT_Mode.
dwtmode(old_modeDWT,'nodisp');

    %---------------------------------------------------------
    function [ca,ch,cv,cd] = decomposeLOC(x)
        % Modification for geck G980781
        y = conv2(double(x),lo);  % 'single' replaced by 'double'
        z = (conv2(y',lo))';
        ca = keepLOC(z,size2D);
        z = (conv2(y',hi))';
        ch = keepLOC(z,size2D);
        y = conv2(double(x),hi);  % 'single' replaced by 'double'
        z = (conv2(y',lo))';
        cv = keepLOC(z,size2D);
        z = (conv2(y',hi))';
        cd = keepLOC(z,size2D);
    end
    %---------------------------------------------------------
    function y = keepLOC(z,siz)
        sz = size(z);
        siz(siz>sz) = sz(siz>sz);
        last = first+siz-1;
        y = z(first(1):last(1),first(2):last(2),:);
    end
    %---------------------------------------------------------

end
