function [W, lambda] = pca2(X, no_dims)
%PCA Perform pricinpal components analysis algorithm
%
%   [W, Fe] = PCA(X, no_dims)
%
%         Input:
%                 X                       - Input Data matrix. Each row vector of data is a data point
%                 no_dims                 - Sets the number of dimensions of the feature points in the embedded feature space
%                
%
%         Output:
%                 W                       - Each column is an embedding function, for a new
%                                           data point (row vector) x,  y = x*eigvector
%                                           will be the embedding result of x.
%                 Fe                      - The coordinates of the low-dimensional data
%
%
%       Copyright notes
%       Author: Wenzhi Liao IPI, Telin, Ghent University, Belgium
%       Date: 25/2/2010


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
	
	% Make sure data is zero mean
    mapping.mean = mean(X, 1);
	X = bsxfun(@minus, X, mapping.mean);

	% Compute covariance matrix
    if size(X, 2) < size(X, 1)
        C = cov(X);
    else
        C = (1 / size(X, 1)) * (X * X');        % if N>D, we better use this matrix for the eigendecomposition
    end
	
	% Perform eigendecomposition of C
	C(isnan(C)) = 0;
	C(isinf(C)) = 0;
    [W, lambda] = eig(C);
    
    % Sort eigenvectors in descending order
    [lambda, ind] = sort(diag(lambda), 'descend');
    if no_dims > size(W, 2)
        no_dims = size(W, 2);
        warning(['Target dimensionality reduced to ' num2str(no_dims) '.']);
    end
	W = W(:,ind(1:no_dims));
    %lambda = lambda(1:no_dims);
    %ttt=zeros(no_dims);
    %for i=1:no_dims
    %    ttt(i,i)=lambda(i);
    %end
    
	
	% Apply mapping on the data
    if ~(size(X, 2) < size(X, 1))
        W = bsxfun(@times, X' * W, (1 ./ sqrt(size(X, 1) .* lambda))');     % normalize in order to get eigenvectors of covariance matrix
    end
   