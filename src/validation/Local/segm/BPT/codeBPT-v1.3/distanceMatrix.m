function [dist_matrix] = distanceMatrix(E1,E2,method)


%% Parametros
if (nargin < 3)
	error('parametros incorrectos');
end
[m1,n1] = size(E1);
[m2,n2] = size(E2);
if (m1 ~= m2)
	error('dimensiones incorrectas');
end

%% matriz de distancias
dist_matrix = zeros(n1,n2);

%% calcular
for i=1:n1
	for j=1:n2
		dist_matrix(i,j) = distance(E1(:,i),E2(:,j),method);
	end
end
