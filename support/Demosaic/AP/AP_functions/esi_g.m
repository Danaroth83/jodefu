%--------------------------------------------------------------------------
%
%  Yue M. Lu
%  Ecole Polytechnique Federale de Lausanne (EPFL)
%
%--------------------------------------------------------------------------
%
%  esi_g.m
%
%  First created: 02-27-2009
%  Last modified: 06-09-2009
%
%--------------------------------------------------------------------------


function g = esi_g(x)

%  Edge-Sensitive Interpolation for the G channel
%
%  INPUT
%    x: the CFA measurement
%
%  OUTPUT
%    g: the estimated full-resolution g channel

[M, N] = size(x);

g = x;

% Interpolation over B pixels
for m = 4 : 2 : M - 2
    for n = 3 : 2: N - 3
        dx_H = 2 * x(m, n) - x(m, n - 2) - x(m, n + 2);
        delta_H = abs(g(m, n - 1) - g(m, n + 1)) + abs(dx_H);
        
        dx_V = 2 * x(m, n) - x(m + 2, n) - x(m - 2, n);
        delta_V = abs(g(m - 1,n) - g(m + 1, n)) + abs(dx_V);
        
        if delta_V > delta_H
            g(m, n) = (g(m, n - 1) + g(m, n + 1)) / 2 + dx_H / 4;
        elseif delta_V < delta_H
            g(m, n) = (g(m - 1,n) + g(m + 1, n)) / 2 + dx_V / 4;
        else
            g(m, n) = (g(m, n - 1) + g(m, n + 1) + g(m - 1,n) + g(m + 1, n)) / 4 ...
                + (dx_H + dx_V) / 8;
        end
    end 
end

% Interpolation over R pixels
for m = 3 : 2 : M - 3
    for n = 4 : 2 : N - 2
        
        dx_H = 2 * x(m, n) - x(m, n - 2) - x(m, n + 2);
        delta_H = abs(g(m, n - 1) - g(m, n + 1)) + abs(dx_H);
        
        dx_V = 2 * x(m, n) - x(m + 2, n) - x(m - 2, n);
        delta_V = abs(g(m - 1,n) - g(m + 1, n)) + abs(dx_V);
        
        if delta_V > delta_H
            g(m, n) = (g(m, n - 1) + g(m, n + 1)) / 2 + dx_H / 4;
        elseif delta_V < delta_H
            g(m, n) = (g(m - 1,n) + g(m + 1, n)) / 2 + dx_V / 4;
        else
            g(m, n) = (g(m, n - 1) + g(m, n + 1) + g(m - 1,n) + g(m + 1, n)) / 4 ...
                + (dx_H + dx_V) / 8;
        end
    end
end

