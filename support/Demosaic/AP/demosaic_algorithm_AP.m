%% ALTERNATING PROJECTIONS DEMOSAICKING
%
% Description:
% Implements the Alternating projections algorithm described in [1]
%
% Usage:
% [I_out,time]=(I_in,approx_rank,niter,K)
%
% Input:
% I_in:        Image mosaicked with a Bayer mask (band order: [G,R;B,G])
% approx_rank: if 0, use the original nonseparable filters (default);
%              if 1, use separable approximations (i.e. rank-one
%              approximations);
%              if 2, use rank-two approximations; and so on;
% niter:       Iteration used to compute the decomposition filters (default: 20)
% K:           Truncated sizes of the filters (default: 6)
%
% Output:
% I_out:       Demosaicked image
% time:        Computation time (excluding the comuputation of the filters)

function [I_out,time]=demosaic_algorithm_AP(I_in,approx_rank,niter,K)

    if nargin<=1 || isempty(approx_rank), approx_rank=0; end
    if nargin<=2 || isempty(niter), niter=20; end
    if nargin<=3 || isempty(K), K=6; end
    
    current_folder=pwd;
    algorithm_folder=fileparts(mfilename('fullpath'));
    cd(fullfile(algorithm_folder,'AP_functions'));
    
    [flts_r, flts_b] = simplify_filters(niter, K);
    [I_out, time] = onestep_demosaicking_conv(I_in, flts_r, flts_b, approx_rank);
    
    cd(current_folder);
    
end