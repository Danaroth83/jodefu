%% Projection on the l_infinity,1 norm ball of radius gamma
% 
% Description:
% This function implements the algorithm described in [1] for calculating
% the projection  of x over the l_inf,1 ball of radius gamma>=0:
%   projball_linf1(x,gamma) = argmin_{x' : linf1norm(x')<=gamma} ||x'-x||_2
% where:
%   linf1norm= @(x) sum(max(abs(x),[],1),2:ndims(x));
% The algorithm is noniterative and exact. Its average complexity, for y
% of size N x M, is O(NM.log(N))
%
% Usage:
% y = projball_l1inf(x,gamma);
%
% Input:
% x: Input Matrix
% gamma: Radius of the projection ball
% 
% Output:
% y: Evaluation of the projection
%
% References:
% [1] L. Condat, "Fast projection onto the simplex and the l1 ball",
% Mathematical Programming, vol. 158, no. 1, pp. 575-585, 2016.

function x = projball_linf1(y, tau)

    sizey=size(y);
    y=reshape(y,size(y,1),[]);
    y=y.';

	x = abs(y);
	v = sum(max(x,[],2));
	if v<=tau, x = y; return; end
	x = sort(x,2,'descend');
	[N,M] = size(y);
	S = cumsum(x,2);
	idx = ones(N,1);
	theta = (v-tau)/N;
	mu = zeros(N,1);
	active = ones(N,1);
	thetaold = 0;
	while thetaold~=theta
		for n=1:N
			if active(n)
				j = idx(n);
				while (j<M)&&((S(n,j)-theta)/j<x(n,j+1)), j=j+1; end
				idx(n) = j;
				mu(n) = S(n,j)/j;
				if (j==M)&&(mu(n)-theta/j<=0), active(n)=0; mu(n)=0; end
			end
		end
		thetaold = theta;
		theta = (sum(mu)-tau)/sum(active./idx);
	end
	x = min(abs(y),(mu-theta./idx).*active).*sign(y);
    
    x=x.';
    x=reshape(x,sizey); 
    
end
