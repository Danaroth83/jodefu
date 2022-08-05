function [d,dist] = disimilitud(E1,E2,method1,method2,alfa)


%% Parametros
if (nargin < 3)
	error('parametros incorrectos');
end
if nargin == 3
    method2 = 'miguel';
    alfa = 0.5;
else
    if nargin == 4
        alfa = 0.5;
    end
end
if nargin >= 5
    if alfa < 0.0 || alfa > 1.0
        error('alfa must be [0,1] valued');
    end
end
[m1,n1] = size(E1);
[m2,n2] = size(E2);
if (m1 ~= m2)
	error('dimensiones incorrectas');
end

%% matriz de distancias
dist = distanceMatrix(E1,E2,method1);

%% disimilitud
switch lower(method2)
    case {'orlando'}
        min_r = mean(min(dist));
        min_c = mean(min(dist'));
        d = (min_r + min_c)*(abs(n1-n2)+1);
    case {'grana'}
        min_r = norm(min(dist));
        min_c = norm(min(dist'));
        d = (min_r + min_c);
    case {'grana2'}
        min_r = norm(min(dist));
        min_c = norm(min(dist'));
        d = (min_r + min_c)*(abs(n1-n2)+1);
    case {'hausdorff'}
        min_r = min(dist);
        min_c = min(dist');
        d = max(max(min_r),max(min_c));
    case {'hausdorff_robust'}
        min_r = sort(min(dist),'ascend');
        min_c = sort(min(dist'),'ascend');
        h_r = ceil(alfa*size(min_r,1));
        h_c = ceil(alfa*size(min_c,1));
        d_r = sum(min_r(1:h_r));
        d_c = sum(min_c(1:h_c));
        d = max(d_r,d_c);
    otherwise % miguel
        min_r = mean(min(dist));
        min_c = mean(min(dist'));
        d = alfa*min_r + (1.0-alfa)*min_c;
end
