%--------------------------------------------------------------------------
%
%  Yue M. Lu
%  Ecole Polytechnique Federale de Lausanne (EPFL)
%
%--------------------------------------------------------------------------
%
%  get_L_polymtx.m
%
%  First created: 01-12-2009
%  Last modified: 06-11-2009
%
%--------------------------------------------------------------------------

function L_polymtx = get_L_polymtx(L, cc, M, N)

%  Generate the full polyphase matrix of the lowpass filter L
%
%  INPUT
%    L: the 2-D lowpass filter
%
%    cc: the coordinates of the center of the filter
%
%    M: the number of rows in the image
%
%    N: the number of columns in the image
%
%  OUTPUT
%    L_polymtx: the full polyphase matrix of L

if any([M N] < size(L)) || any(mod([M N], 2) ~= 0)
    error('N must be an even number larger than the dimension of the filter L.');
end

L_flt = zeros(M, N);
L_flt(1:size(L,1), 1:size(L,2)) = L;
L_flt = circshift(L_flt, [1 1] - [cc(1) cc(2)]);

% Obtain the four polyphase components (downsampling by 2 along each
% dimension)
L_poly = polyphase_components(L_flt, 2, 2, 's', 0);

% L_polymtx = 
%
% | L00(z)  L01(z)z_2^{-1}  L10(z)z_1^{-1}  L11(z)z_1^{-1)z_2^{-1} |
% |                                                                |
% | L01(z)  L00(z)          L11(z)z_1^{-1}  L10(z)z_1^{-1}         |
% |                                                                |  
% | L10(z)  L11(z)z_2^{-1}  L00(z)          L01(z)z_2^{-1}         |
% |                                                                |
% | L11(z)  L10(z)          L01(z)          L00(z)                 |
%

L00 = fftn(L_poly{1,1});

L_polymtx = zeros(4, 4, numel(L00));

L_polymtx(1, 1, :) = L00(:);
L_polymtx(2, 2, :) = L00(:);
L_polymtx(3, 3, :) = L00(:);
L_polymtx(4, 4, :) = L00(:);

L10 = fftn(L_poly{2, 1});
L_polymtx(3, 1, :) = L10(:);
L_polymtx(4, 2, :) = L10(:);

L10s = fftn(circshift(L_poly{2, 1}, [1 0]));
L_polymtx(1, 3, :) = L10s(:);
L_polymtx(2, 4, :) = L10s(:);

L01 = fftn(L_poly{1, 2});
L_polymtx(2, 1, :) = L01(:);
L_polymtx(4, 3, :) = L01(:);

L01s = fftn(circshift(L_poly{1, 2}, [0 1]));
L_polymtx(1, 2, :) = L01s(:);
L_polymtx(3, 4, :) = L01s(:);

L11 = fftn(L_poly{2, 2});
L_polymtx(4, 1, :) = L11(:);

L11s = fftn(circshift(L_poly{2, 2}, [1 1]));
L_polymtx(1, 4, :) = L11s(:);

L11s = fftn(circshift(L_poly{2, 2}, [1 0]));
L_polymtx(2, 3, :) = L11s(:);

L11s = fftn(circshift(L_poly{2, 2}, [0 1]));
L_polymtx(3, 2, :) = L11s(:);
