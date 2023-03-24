%--------------------------------------------------------------------------
%
%  Yue M. Lu
%  Ecole Polytechnique Federale de Lausanne (EPFL)
%
%--------------------------------------------------------------------------
%
%  simplify_filters.m
%
%  First created: 03-06-2009
%  Last modified: 06-09-2009
%
%--------------------------------------------------------------------------

function [flts_r, flts_b] = simplify_filters(niter, K)

%  Retrieve and simplify the pre-computed polyphase filters used in the
%  one-step implementation of the AP algorithm
%
%  INPUT
%    niter: the number of iterations used by the AP algorithm
%
%    K: We truncate all filters to the size of K-by-K
%
%  OUTPUT
%    flts_r: the obtained polyphase filters for the red channel. See NOTE
%    for more details about the data structure.
%
%    flts_b: the obtained polyphase filters for the blue channel
%
%  NOTE
%    Both flts_r and flts_b are 3-by-4 cell arrays. Each cell corresponds
%    to one of the polyphase filters (see generate_poly_filters for
%    details). Every such polyphase filter can be represented by a K-by-K
%    matrix (after truncation), and is stored as follows.
%
%    {1st left singular vector of the matrix, 1st right singular vector,
%    2nd left, 2nd right, ..., kth left, kth right, and finally, the 2-D
%    matrix itself}
%
%    The singular vectors can be used to construct the separable (or in
%    general,low-rank) approximations of the original nonseparable filters.

% Construct the full path of the mat file containing the pre-computed
% filter coefficients
p = fileparts(mfilename('fullpath'));
fname = fullfile(p,'Filters',['poly_flts_niter_',num2str(niter),'.mat']);

if exist(fname, 'file') ~= 2
    % error(['First run precompute_filters.m (with parameter "niter" set to ' num2str(niter) ') to generate the polyphase filters.']);
    precompute_filters(niter);
end

load(fname,'r_flt_mtx','b_flt_mtx'); % obtain two filter structures: r_flt_mtx and b_flt_mtx

flt_tmp = cell(3, 4);

for ch = 1 : 2
    % red: ch = 1; blue: ch = 2
    if ch == 1
        flt_mtx = r_flt_mtx;
    else
        flt_mtx = b_flt_mtx;
    end
    
    for m = 1 : 3
        for n = 1 : 4
            flt = flt_mtx{m, n};
            % truncation and low-rank (separable) approximation
            flt = circshift(flt, floor(K/2)*[1 1]);
            flt_tmp{m, n} = lowrank_approx(flt(1:K, 1:K));
        end
    end
    
    if ch == 1
        flts_r = flt_tmp;
    else
        flts_b = flt_tmp;
    end
end


function flt_lowrank = lowrank_approx(flt)

K = size(flt, 1);

[u, s, v] = svd(flt);

flt_lowrank = cell(K * 2 + 1, 1);
for n = 1 : K
    % 1st left eigenvector, 1st right eigenvector, 2nd left, 2nd right, ...
    flt_lowrank{2*n-1} = u(:, n) * s(n, n);
    flt_lowrank{2*n} = v(:, n);
end

% the original 2-D nonseparable filter itself
flt_lowrank{2*K+1} = flt;

