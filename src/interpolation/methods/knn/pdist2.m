function [D,I] = pdist2(X,Y,dist,varargin)
%PDIST2 Pairwise distance between two sets of observations.
%   D = PDIST2(X,Y) returns a matrix D containing the Euclidean distances
%   between each pair of observations in the MX-by-N data matrix X and
%   MY-by-N data matrix Y. Rows of X and Y correspond to observations,
%   and columns correspond to variables. D is an MX-by-MY matrix, with the
%   (I,J) entry equal to distance between observation I in X and
%   observation J in Y.
%
%   D = PDIST2(X,Y,DISTANCE) computes D using DISTANCE.  Choices are:
%
%       'euclidean'        - Euclidean distance (default)
%       'squaredeuclidean' - Squared Euclidean distance
%       'seuclidean'       - Standardized Euclidean distance. Each
%                            coordinate difference between rows in X and Y
%                            is scaled by dividing by the corresponding
%                            element of the standard deviation computed
%                            from X, S=NANSTD(X). To specify another value
%                            for S, use
%                            D = PDIST2(X,Y,'seuclidean',S).
%       'cityblock'        - City Block distance
%       'minkowski'        - Minkowski distance. The default exponent is 2.
%                            To specify a different exponent, use
%                            D = PDIST2(X,Y,'minkowski',P), where the
%                            exponent P is a scalar positive value.
%       'chebychev'        - Chebychev distance (maximum coordinate
%                            difference)
%       'mahalanobis'      - Mahalanobis distance, using the sample
%                            covariance of X as computed by NANCOV.  To
%                            compute the distance with a different
%                            covariance, use
%                            D = PDIST2(X,Y,'mahalanobis',C), where the
%                            matrix C is symmetric and positive definite.
%       'cosine'           - One minus the cosine of the included angle
%                            between observations (treated as vectors)
%       'correlation'      - One minus the sample linear correlation
%                            between observations (treated as sequences of
%                            values).
%       'spearman'         - One minus the sample Spearman's rank
%                            correlation between observations (treated as
%                            sequences of values)
%       'hamming'          - Hamming distance, percentage of coordinates
%                            that differ
%       'jaccard'          - One minus the Jaccard coefficient, the
%                            percentage of nonzero coordinates that differ
%       function           - A distance function specified using @, for
%                            example @DISTFUN
%
%   A distance function must be of the form
%
%         function D2 = DISTFUN(ZI,ZJ),
%
%   taking as arguments a 1-by-N vector ZI containing a single observation
%   from X or Y, an M2-by-N matrix ZJ containing multiple observations from
%   X or Y, and returning an M2-by-1 vector of distances D2, whose Jth
%   element is the distance between the observations ZI and ZJ(J,:).
%
%   For built-in distance metrics, the distance between observation I in X
%   and observation J in Y will be NaN if observation I in X or observation
%   J in Y contains NaNs.
%
%   D = PDIST2(X,Y,DISTANCE,'Smallest',K) returns a K-by-MY matrix D
%   containing the K smallest pairwise distances to observations in X for
%   each observation in Y. PDIST2 sorts the distances in each column of D
%   in ascending order. D = PDIST2(X,Y,DISTANCE, 'Largest',K) returns the K
%   largest pairwise distances sorted in descending order. If K is greater
%   than MX, PDIST2 returns an MX-by-MY distance matrix. For each
%   observation in Y, PDIST2 finds the K smallest or largest distances by
%   computing and comparing the distance values to all the observations in
%   X.
%
%   [D,I] = PDIST2(X,Y,DISTANCE,'Smallest',K) returns a K-by-MY matrix I
%   containing indices of the observations in X corresponding to the K
%   smallest pairwise distances in D. [D,I] = PDIST2(X,Y,DISTANCE,
%   'Largest',K) returns indices corresponding to the K largest pairwise
%   distances.
%
%   Example:
%      % Compute the ordinary Euclidean distance
%      X = randn(100, 5);
%      Y = randn(25, 5);
%      D = pdist2(X,Y,'euclidean');         % euclidean distance
%
%      % Compute the Euclidean distance with each coordinate difference
%      % scaled by the standard deviation
%      Dstd = pdist2(X,Y,'seuclidean');
%
%      % Use a function handle to compute a distance that weights each
%      % coordinate contribution differently.
%      Wgts = [.1 .3 .3 .2 .1];            % coordinate weights
%      weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
%      Dwgt = pdist2(X,Y, @(Xi,Xj) weuc(Xi,Xj,Wgts));
%
%   See also PDIST, KNNSEARCH, CREATENS, KDTreeSearcher,
%            ExhaustiveSearcher.

%   An example of distance for data with missing elements:
%
%      X = randn(100, 5);     % some random points
%      Y = randn(25, 5);      % some more random points
%      X(unidrnd(prod(size(X)),1,20)) = NaN; % scatter in some NaNs
%      Y(unidrnd(prod(size(Y)),1,5)) = NaN; % scatter in some NaNs
%      D = pdist2(X, Y, @naneucdist);
%
%      function D = naneucdist(XI, YJ) % euclidean distance, ignoring NaNs
%      [m,p] = size(YJ);
%      sqdxy = bsxfun(@minus,XI,YJ) .^ 2;
%      pstar = sum(~isnan(sqdxy),2); % correction for missing coordinates
%      pstar(pstar == 0) = NaN;
%      D = sqrt(nansum(sqdxy,2) .* p ./ pstar);
%
%
%   For a large number of observations, it is sometimes faster to compute
%   the distances by looping over coordinates of the data (though the code
%   is more complicated):
%
%      function D = nanhamdist(XI, YJ) % hamming distance, ignoring NaNs
%      [m,p] = size(YJ);
%      nesum = zeros(m,1);
%      pstar = zeros(m,1);
%      for q = 1:p
%          notnan = ~(isnan((XI(q)) | isnan(YJ(:,q)));
%          nesum = nesum + (XI(q) ~= YJ(:,q)) & notnan;
%          pstar = pstar + notnan;
%      end
%      D = nesum ./ pstar;

%   Copyright 2009-2018 The MathWorks, Inc.


if nargin > 0
    X = convertStringsToChars(X);
end

if nargin > 1
    Y = convertStringsToChars(Y);
end

if nargin > 2
    dist = convertStringsToChars(dist);
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error(message('stats:pdist2:TooFewInputs'));
end

if ~ismatrix(X) || ~ismatrix(Y)
    error(message('stats:pdist2:UnsupportedND'));
end

[nx,p] = size(X);
[ny,py] = size(Y);
if py ~= p
    error(message('stats:pdist2:SizeMismatch'));
end

additionalArg = [];

if nargin < 3
    dist = 'euc';
else %distance is provided
    if ischar(dist)
        methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
            'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
            'spearman'; 'hamming'; 'jaccard'; 'squaredeuclidean'};
        i = find(strncmpi(dist,methods,length(dist)));
        if length(i) > 1
            error(message('stats:pdist2:AmbiguousDistance', dist));
        elseif isempty(i)
            error(message('stats:pdist2:UnrecognizedDistance', dist));
        else
            if i == 12 %'squaredeuclidean'
                dist = 'sqe'; % represent squared Euclidean
            else
                dist = methods{i}(1:3);
            end
            
            if ~isempty(varargin)
                arg = varargin{1};
                
                % Get the additional distance argument from the inputs
                if isnumeric(arg)
                    switch dist
                        case {'seu' 'mah' 'min'}
                            additionalArg = arg;
                            varargin = varargin(2:end);
                    end
                end
            end
        end
    elseif isa(dist, 'function_handle')
        distfun = dist;
        dist = 'usr';
    else
        error(message('stats:pdist2:BadDistance'));
    end
end

pnames = {'smallest' 'largest' 'radius'};
dflts =  {        []        []  []};
[smallest,largest,radius] = parseyarg(pnames, dflts, varargin{:});

smallestLargestFlag = [];
if sum([~isempty(largest) ~isempty(smallest) ~isempty(radius)]) > 1
    error(message('stats:pdist2:SmallestAndLargest'));
end
if ~isempty(smallest)
    if ~(isscalar(smallest) && isnumeric(smallest) && smallest >= 1 && round(smallest) == smallest)
        error(message('stats:pdist2:BadSmallest'));
    end
    smallestLargestFlag = min(smallest,nx);
elseif ~isempty(largest)
    if ~(isscalar(largest) && isnumeric(largest) && largest >= 1 && round(largest) == largest)
        error(message('stats:pdist2:BadLargest'));
    end
    smallestLargestFlag = -min(largest,nx);
elseif ~isempty(radius)
    if ~(isscalar(radius) && isnumeric(radius) && radius >= 0 )
        error(message('stats:pdist2:BadRadius'));
    end
elseif nargout > 1
    error(message('stats:pdist2:TooManyOutputs'));
end

% For a built-in distance, integer/logical/char/anything data will be
% converted to float. Complex floating point data can't be handled by
% a built-in distance function.
%if ~strcmp(dist,'usr')

try
    outClass = superiorfloat(X,Y);
catch
    if isfloat(X)
        outClass = class(X);
    elseif isfloat(Y)
        outClass = class(Y);
    else
        outClass = 'double';
    end
end

if ~strcmp(dist,'usr')
    if ~strcmp(class(X),outClass) || ~strcmp(class(Y),outClass)
        warning(message('stats:pdist2:DataConversion', outClass));
    end
    X = cast(X,outClass);
    Y = cast(Y,outClass);
    if  ~isreal(X) || ~isreal(Y)
        error(message('stats:pdist2:ComplexData'));
    end
end

% Degenerate case, just return an empty of the proper size.
if (nx == 0) || (ny == 0)
    if ~isempty(radius)
        D = repmat({zeros(1,0, outClass)},1,ny);
        I = repmat({zeros(1,0)},1,ny);
    else
        if ~isempty(smallestLargestFlag)
            nD = abs(smallestLargestFlag);
        else
            nD = nx;
        end
        D = zeros(nD,ny,outClass); % X and Y were single/double, or cast to double
        I = zeros(nD,ny);
    end
    return;
end

switch dist
    case 'seu' % Standardized Euclidean weights by coordinate variance
        if isempty(additionalArg)
            additionalArg =  nanvar(X,[],1);
            if any(additionalArg == 0)
                warning(message('stats:pdist2:ConstantColumns'));
            end
            additionalArg = 1./ additionalArg;
        else
            if ~(isvector(additionalArg) && length(additionalArg) == p...
                    && all(additionalArg >= 0))
                error(message('stats:pdist2:InvalidWeights'));
            end
            if any(additionalArg == 0)
                warning(message('stats:pdist2:ZeroInverseWeights'));
            end
            additionalArg = 1./ (additionalArg .^2);
        end
        
    case 'mah' % Mahalanobis
        if isempty(additionalArg)
            if nx == 1
                error(message('stats:pdist2:tooFewXRowsForMah'));
            end
            additionalArg = nancov(X);
            [T,flag] = chol(additionalArg);
        else %provide the covariance for mahalanobis
            if ~isequal(size(additionalArg),[p,p])
                error(message('stats:pdist2:InvalidCov'));
            end
            %cholcov will check whether the covariance is symmetric
            [T,flag] = cholcov(additionalArg,0);
        end
        
        if flag ~= 0
            error(message('stats:pdist2:InvalidCov'));
        end
        
        if ~issparse(X) && ~issparse(Y)
            additionalArg = T \ eye(p); %inv(T)
        end
        
    case 'min' % Minkowski
        if isempty(additionalArg)
            additionalArg = 2;
        elseif ~(isscalar(additionalArg) && additionalArg > 0)
            error(message('stats:pdist2:BadMinExp'));
        end
    case 'cos' % Cosine
        [X,Y,flag] = normalizeXY(X,Y);
        if flag
            warning(message('stats:pdist2:ZeroPoints'));
        end
        
    case 'cor' % Correlation
        X = bsxfun(@minus,X,mean(X,2));
        Y = bsxfun(@minus,Y,mean(Y,2));
        [X,Y,flag] = normalizeXY(X,Y);
        if flag
            warning(message('stats:pdist2:ConstantPoints'));
        end
        
    case 'spe'
        X = tiedrank(X')'; % treat rows as a series
        Y = tiedrank(Y')';
        X = X - (p+1)/2; % subtract off the (constant) mean
        Y = Y - (p+1)/2;
        [X,Y,flag] = normalizeXY(X,Y);
        if flag
            warning(message('stats:pdist2:TiedPoints'));
        end
        
    otherwise
        
end

% Note that if the above switch statement is modified to include the
% 'che', 'euc', or 'cit' distances, that code may need to be repeated
% in the corresponding block below.
if strcmp(dist,'min') % Minkowski distance
    if isinf(additionalArg) %the exponent is inf
        dist = 'che';
        additionalArg = [];
    elseif additionalArg == 2 %the exponent is 2
        dist = 'euc';
        additionalArg = [];
    elseif additionalArg == 1 %the exponent is 1
        dist = 'cit';
        additionalArg = [];
    end
end

% Call a mex file to compute distances for the build-in distance measures
% on non-sparse real float (double or single) data.
if ~strcmp(dist,'usr') && (~issparse(X) && ~issparse(Y))
    additionalArg = cast(additionalArg,outClass);
    if strcmp(dist,'sqe')
       radius = sqrt(radius); 
    end
    if nargout < 2
        D = internal.stats.pdist2mex(X',Y',dist,additionalArg,smallestLargestFlag,radius);
    else
        [D,I] = internal.stats.pdist2mex(X',Y',dist,additionalArg,smallestLargestFlag,radius);
    end
    
    % The following MATLAB code implements the same distance calculations as
    % the mex file. It assumes X and Y are real single or double.  It is
    % currently only called for sparse inputs, but it may also be useful as a
    % template for customization.
elseif ~strcmp(dist,'usr')
    if any(strcmp(dist, {'ham' 'jac' 'che'}))
        xnans = any(isnan(X),2);
        ynans = any(isnan(Y),2);
    end
    
    if ~isempty(radius)
        D = cell(1,ny);
        I = cell(1,ny);
    elseif isempty(smallestLargestFlag)
        D = zeros(nx,ny,outClass);
    else
        D = zeros(abs(smallestLargestFlag),ny,outClass);
        I = zeros(abs(smallestLargestFlag),ny);
        
    end
    
    switch dist
        case 'euc'  % Euclidean
            if isempty(radius) && isempty(smallestLargestFlag) && (issparse(X) && issparse(Y))
                D = internal.stats.pdist2SparseMEX(X',Y',dist,additionalArg);
            else
                for i = 1:ny
                    dsq = zeros(nx,1,outClass);
                    for q = 1:p
                        dsq = dsq + (X(:,q) - Y(i,q)).^2;
                    end
                    dsq = sqrt(dsq);
                    if ~isempty(radius)
                        [D{i},I{i}] = radiusSort(dsq,radius);
                    elseif isempty(smallestLargestFlag)
                        D(:,i) = dsq;
                    else
                        [D(:,i),I(:,i)] = partialSort(dsq,smallestLargestFlag);
                    end
                end
            end
        case 'sqe'  % Squared Euclidean
            if isempty(radius) && isempty(smallestLargestFlag) && (issparse(X) && issparse(Y))
                D = internal.stats.pdist2SparseMEX(X',Y',dist,additionalArg);
            else
                for i = 1:ny
                    dsq = zeros(nx,1,outClass);
                    for q = 1:p
                        dsq = dsq + (X(:,q) - Y(i,q)).^2;
                    end
                    
                    if ~isempty(radius)
                        [D{i},I{i}] = radiusSort(dsq,radius);
                    elseif isempty(smallestLargestFlag)
                        D(:,i) = dsq;
                    else
                        [D(:,i),I(:,i)] = partialSort(dsq,smallestLargestFlag);
                    end
                end
            end
        case 'seu'    % Standardized Euclidean
            wgts = additionalArg;
            for i = 1:ny
                dsq = zeros(nx,1,outClass);
                for q = 1:p
                    dsq = dsq + wgts(q) .* (X(:,q) - Y(i,q)).^2;
                end
                dsq = sqrt(dsq);
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(dsq,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) = dsq;
                else
                    [D(:,i),I(:,i)] = partialSort(dsq,smallestLargestFlag);
                end
            end
            
        case 'cit'    % City Block
            for i = 1:ny
                dsq = zeros(nx,1,outClass);
                for q = 1:p
                    dsq = dsq + abs(X(:,q) - Y(i,q));
                end
                
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(dsq,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) = dsq;
                else
                    [D(:,i),I(:,i)] = partialSort(dsq,smallestLargestFlag);
                end
            end
            
        case 'mah'    % Mahalanobis
            for i = 1:ny
                del = bsxfun(@minus,X,Y(i,:));
                dsq = sum((del/T) .^ 2, 2);
                dsq = sqrt(dsq);
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(dsq,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) = dsq;
                else
                    [D(:,i),I(:,i)] = partialSort(dsq,smallestLargestFlag);
                end
            end
            
        case 'min'    % Minkowski
            expon = additionalArg;
            for i = 1:ny
                dpow = zeros(nx,1,outClass);
                for q = 1:p
                    dpow = dpow + abs(X(:,q) - Y(i,q)).^expon;
                end
                dpow = dpow .^ (1./expon);
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(dpow,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) = dpow;
                else
                    [D(:,i),I(:,i)] = partialSort(dpow,smallestLargestFlag);
                end
                
            end
            
        case {'cos' 'cor' 'spe'}   % Cosine, Correlation, Rank Correlation
            % This assumes that data have been appropriately preprocessed
            for i = 1:ny
                d = zeros(nx,1,outClass);
                for q = 1:p
                    d = d + (X(:,q).*Y(i,q));
                end
                d(d>1) = 1; % protect against round-off, don't overwrite NaNs
                d = 1 - d;
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(d,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) = d;
                else
                    [D(:,i),I(:,i)] = partialSort(d,smallestLargestFlag);
                end
                
            end
        case 'ham'    % Hamming
            for i = 1:ny
                nesum = zeros(nx,1,outClass);
                for q = 1:p
                    nesum = nesum + (X(:,q) ~= Y(i,q));
                end
                nesum(xnans|ynans(i)) = NaN;
                nesum = (nesum ./ p);
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(nesum,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) = nesum;
                else
                    [D(:,i),I(:,i)] = partialSort(nesum,smallestLargestFlag);
                end
                
            end
        case 'jac'    % Jaccard
            for i = 1:ny
                nzsum = zeros(nx,1,outClass);
                nesum = zeros(nx,1,outClass);
                for q = 1:p
                    nz = (X(:,q) ~= 0 | Y(i,q) ~= 0);
                    ne = (X(:,q) ~= Y(i,q));
                    nzsum = nzsum + nz;
                    nesum = nesum + (nz & ne);
                end
                nesum(xnans | ynans(i)) = NaN;
                d = (nesum ./ nzsum);
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(d,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) = d;
                else
                    [D(:,i),I(:,i)] = partialSort(d,smallestLargestFlag);
                end
                
            end
        case 'che'    % Chebychev
            for i = 1:ny
                dmax = zeros(nx,1,outClass);
                for q = 1:p
                    dmax = max(dmax, abs(X(:,q) - Y(i,q)));
                end
                dmax(xnans | ynans(i)) = NaN;
                if ~isempty(radius)
                    [D{i},I{i}] = radiusSort(dmax,radius);
                elseif isempty(smallestLargestFlag)
                    D(:,i) =  dmax;
                else
                    [D(:,i),I(:,i)] = partialSort(dmax,smallestLargestFlag);
                end
            end
    end
    
    % Compute distances for a caller-defined distance function.
else % if strcmp(dist,'usr')
    try
        D = feval(distfun,Y(1,:),X(1,:));
    catch ME
        if strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                && ~isempty(strfind(ME.message, func2str(distfun)))
            error(message('stats:pdist2:DistanceFunctionNotFound', func2str( distfun )));
        end
        % Otherwise, let the catch block below generate the error message
        D = [];
    end
    
    if ~isnumeric(D)
        error(message('stats:pdist2:OutputBadType'));
    end
    
    if ~isempty(radius)
        D = cell(1,ny);
        if nargout >= 2
            I = cell(1,ny);
        end
        
        for i = 1:ny
            try
                temp = feval(distfun,Y(i,:),X);
            catch ME
                if isa(distfun, 'inline')
                    m = message('stats:pdist2:DistanceInlineError');
                    ME2 = MException(m.Identifier,'%s',getString(m));
                    throw(addCause(ME2,ME));
                else
                    m = message('stats:pdist2:DistanceFunctionError',func2str(distfun));
                    ME2 = MException(m.Identifier,'%s',getString(m));
                    throw(addCause(ME2,ME));
                end
            end
            
            if nargout < 2
                D{i} = radiusSort(temp,radius);
            else
                [D{i},I{i}] = radiusSort(temp, radius);
            end
        end
    elseif ~isempty(smallestLargestFlag)
        D = zeros(abs(smallestLargestFlag),ny,class(D));
        if nargout > 1
            I = zeros(abs(smallestLargestFlag),ny,class(D));
        end
        
        for i = 1:ny
            
            try
                temp = feval(distfun,Y(i,:),X);
            catch ME
                if isa(distfun, 'inline')
                    m = message('stats:pdist2:DistanceInlineError');
                    ME2 = MException(m.Identifier,'%s',getString(m));
                    throw(addCause(ME2,ME));
                else
                    m = message('stats:pdist2:DistanceFunctionError',func2str(distfun));
                    ME2 = MException(m.Identifier,'%s',getString(m));
                    throw(addCause(ME2,ME));
                end
            end
            
            if nargout < 2
                D(:,i) = partialSort(temp,smallestLargestFlag);
            else
                [D(:,i),I(:,i)] = partialSort(temp,smallestLargestFlag);
            end
        end
        
    else  %compute all the pairwise distance
        % Make the return have whichever numeric type the distance function
        % returns.
        D = zeros(nx,ny,class(D));
        
        for i = 1:ny
            try
                D(:,i) = feval(distfun,Y(i,:),X);
                
            catch ME
                if isa(distfun, 'inline')
                    m = message('stats:pdist2:DistanceInlineError');
                    ME2 = MException(m.Identifier,'%s',getString(m));
                    throw(addCause(ME2,ME));
                else
                    m = message('stats:pdist2:DistanceFunctionError',func2str(distfun));
                    ME2 = MException(m.Identifier,'%s',getString(m));
                    throw(addCause(ME2,ME));
                end
            end
        end
    end
    
end

%---------------------------------------------
% Normalize the data matrices X and Y to have unit norm
function [X,Y,flag] = normalizeXY(X,Y)
Xmax = max(abs(X),[],2);
X2 = bsxfun(@rdivide,X,Xmax);
Xnorm = sqrt(sum(X2.^2, 2));

Ymax = max(abs(Y),[],2);
Y2 = bsxfun(@rdivide,Y,Ymax);
Ynorm = sqrt(sum(Y2.^2, 2));
% Find out points for which distance cannot be computed.

% The norm will be NaN for rows that are all zeros, fix that for the test
% below.
Xnorm(Xmax==0) = 0;
Ynorm(Ymax==0) = 0;

% The norm will be NaN for rows of X that have any +/-Inf. Those should be
% Inf, but leave them as is so those rows will not affect the test below.
% The points can't be normalized, so any distances from them will be NaN
% anyway.

% Find points that are effectively zero relative to the point with largest norm.
flag =  any(Xnorm <= eps(max(Xnorm))) || any(Ynorm <= eps(max(Ynorm)));
Xnorm = Xnorm .* Xmax;
Ynorm = Ynorm .* Ymax;
X = bsxfun(@rdivide,X,Xnorm);
Y = bsxfun(@rdivide,Y,Ynorm);


function [D,I] = partialSort(D,smallestLargest)
if smallestLargest > 0
    n = smallestLargest;
else
    %sort(D,'descend') puts the NaN values at the beginning of the sorted list.
    %That is not what we want here.
    D = D*-1;
    n = -smallestLargest;
end

if nargout < 2
    D = sort(D,1);
    D = D(1:n,:);
else
    [D,I] = sort(D,1);
    D = D(1:n,:);
    I = I(1:n,:);
end

if smallestLargest < 0
    D = D * -1;
end

function [D,I] = radiusSort(D,radius)
I = find (D <= radius);
D = D(I);
if nargout < 2
    D = sort(D,1)'; %return a row vector
else
    [D,I2] = sort(D,1);
    D= D';
    I = I(I2)';
end


