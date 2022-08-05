function [distance] = distance(e1,e2,method)
%% [distance] = distance(e1,e2,method)
%  Calculates the distance between two given endmembers e1 and e2
%  e1 and e2 are vectors with L-dimensionality, they can be row or column vectors but both have to eb of the same type.
%  method indicates the metric used to calculate the distance. The possible options are:
%    * 'ed': the Euclidean Distance (default)
%    * 'cbd': the City Block Distance
%    * 'td': the Tchebyshev Distance
%    * 'sam': the Spectral Angle Map (in degrees)


%% parameters
if nargin < 2
    error('Insuficient parameters');
end
if nargin == 2
    method = 'ed';
end
[n1 m1] = size(e1);
[n2 m2] = size(e2);
if n1 ~= n2 || m1 ~= m2
    error('Incorrect endmembers size');
end

%% Distance
switch lower(method)
    case {'cbd'}
        distance = sum(abs(e1-e2));
    case {'td'}
        distance = max(abs(e1-e2));
    case {'sam'}
        distance = real(acosd(sum(e1.*e2)/(sqrt(sum(power(e1,2))).*sqrt(sum(power(e2,2))))));
        if isnan(distance)
             distance = sqrt(sum(power(e1-e2,2)));
        end
    otherwise % euclidean
        distance = sqrt(sum(power(e1-e2,2)));
end