%--------------------------------------------------------------------------
%
%  Yue M. Lu
%  Ecole Polytechnique Federale de Lausanne (EPFL)
%
%--------------------------------------------------------------------------
%
%  generate_poly_filters.m
%
%  First created: 01-12-2009
%  Last modified: 06-09-2009
%
%--------------------------------------------------------------------------

function [r_flt_mtx, b_flt_mtx] = generate_poly_filters(L, cc, niter, szF)

%  Generate the polyphase filters used in the one-step implementation
%
%  INPUT
%    L: the lowpass filter used by the AP algorithm
%    
%    cc: coordinate of the center of L
% 
%    niter: number of iterations
%
%    szF: the size of the images (not essential, as long as szF is chosen
%    to be large enough)
%
%  OUTPUT
%    r_flt_mtx: a 3-by-4 cell matrix. The submatrix formed by the first,
%    third, and fourth column contains the filter for the MIMO block, while
%    the second column contains the three channel filters F_00, F_10, F_10.
%    
%    b_flt_mtx: a 3-by-4 cell matrix. The submatrix formed by the first,
%    second, and fourth column contains the filter for the MIMO block, while
%    the third column contains the three channel filters.


if any(rem(szF, 2))
    error('szF must contain even numbers only.');
end

M = szF(1);
N = szF(2);

% obtain the full polyphase matrix of L
L_polymtx = get_L_polymtx(L, cc, M, N);

nfreq = size(L_polymtx, 3);

flt_mtx = cell(3);
flt_vec = cell(3, 1);

for ch = 1 : 2
    % red: ch = 1; blue: ch = 2
    if ch == 1
        T_mtx = L_polymtx([1 3 4], [1 3 4], :);
        V = L_polymtx([1 3 4], 2, :);
    else
        T_mtx = L_polymtx([1 2 4], [1 2 4], :);
        V = L_polymtx([1 2 4], 3, :);
    end
    
    A_mtx = repmat(eye(3), [1 1 nfreq]);
    B_mtx = zeros(size(T_mtx));
    C_vec = zeros(size(V));
    
    % A = T^{niter}
    % B = T^{niter -1} + T^{niter - 2} + ... + T + I
    % C = B * V
    
    for n = 1 : niter
        B_mtx = B_mtx + A_mtx;
        
        for nf = 1 : nfreq
            A_mtx(:,:,nf) = A_mtx(:,:,nf) * T_mtx(:,:,nf);
        end
    end
    
    for nf = 1 : nfreq
        C_vec(:, 1, nf) = B_mtx(:,:, nf) * V(:, 1, nf);
    end
    
    for m = 1 : 3
        for n = 1 : 3
            % getting rid of those coefficients that are extremely small
            flt_mtx{m, n} = truncate_flt(real(ifftn(reshape(A_mtx(m, n, :), [M/2, N/2]))));
        end
    end
    
    for m = 1 : 3
        % getting rid of those coefficients that are extremely small
        flt_vec{m} = truncate_flt(real(ifftn(reshape(C_vec(m, 1, :), [M/2, N/2]))));
    end
    
    if ch == 1
        r_flt_mtx = cell(3, 4);
        r_flt_mtx(:, 1) = flt_mtx(:, 1);
        r_flt_mtx(:, 2) = flt_vec;
        r_flt_mtx(:, 3:4) = flt_mtx(:, 2:3);
    else
        b_flt_mtx = cell(3, 4);
        b_flt_mtx(:, 1:2) = flt_mtx(:, 1:2);
        b_flt_mtx(:, 3) = flt_vec;
        b_flt_mtx(:, 4) = flt_mtx(:, 3);
    end
    
end

function y = truncate_flt(x)

% Get rid of extremely small coefficients in the filter, to save hard disk
% space for storing the filter.

y = x;
y(abs(x) < 1e-10) = 0;