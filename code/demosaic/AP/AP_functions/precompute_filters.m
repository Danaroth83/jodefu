%--------------------------------------------------------------------------
%
%  Yue M. Lu
%  Ecole Polytechnique Federale de Lausanne (EPFL)
%
%--------------------------------------------------------------------------
%
%  precompute_filters.m
%
%  First created: 03-05-2009
%  Last modified: 06-09-2009
%
%--------------------------------------------------------------------------

% Pre-Compute the polyphase filters used in the one-step implementation of
% the AP algorithm.
%
% Parameters to be specified:
%
% niter: the number of iterations used by the AP algorithm
% L (optional): the lowpass filter used by the AP algorithm


function precompute_filters(niter)

    % number of iterations
    if nargin<=0 || isempty(niter), niter = 20; end

    % lowpass filter used in the original AP algorithm
    L = conv([1 2 1]/4, [-1 2 6 2 -1]/8);
    L = kron(L, L'); % get the 2-D separable filter

    cc = (size(L) + 1)/2;
    szF = [512 768];

    % Compute the polyphase filters
    [r_flt_mtx, b_flt_mtx] = generate_poly_filters(L, cc, niter, szF);

    % Saving the computed filter coefficients in a mat file
    p = fileparts(mfilename('fullpath'));

    fname = fullfile(p,'Filters',['poly_flts_niter_',num2str(niter)]);
    % [p(1:length(p)-25) 'Filters/poly_flts_niter_' num2str(niter)];
    save(fname, 'r_flt_mtx', 'b_flt_mtx');
end
