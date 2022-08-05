function [d,dist] = disimilitud_ab(E1,E2,P1,P2,method)


%% Parametros
if (nargin < 4)
	error('parametros incorrectos');
end
if nargin < 5
    method = 'sam';
end

[n1,m1] = size(E1);
[n2,m2] = size(E2);
k1 = max(size(P1));
k2 = max(size(P2));
if (m1 ~= m2) || (n1 ~= k1) || (n2 ~= k2)
	error('dimensiones incorrectas');
end

%% matriz de distancias
dist = distanceMatrix(E1',E2',method);

%% disimilitud
% calculate significant matrix
s = zeros(size(dist));
l = dist;
while (max(P1) > 0.001 && max(P2) > 0.001)
    [a,b] = size(l);
    if (a == 1)
        i = 1;
        [m,j] = min(l);
    elseif (b == 1)
        [m,i] = min(l);
        j = 1;
    else
        [m,i] = min(l);
        [m,j] = min(m);
        i = i(j);
    end
    k = min(P1(i),P2(j));
    s(i,j) = k;
    P1(i) = P1(i) - k;
    P2(j) = P2(j) - k;
    l(i,j) = Inf;
end
% calculate dissimilarities
d = sum(sum(dist.*s));

