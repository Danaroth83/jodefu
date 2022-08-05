%--------------------------------------------------------------------------
%
%  Yue M. Lu
%  Ecole Polytechnique Federale de Lausanne (EPFL)
%
%--------------------------------------------------------------------------
%
%  polyphase_components.m
%
%  First created: 03-26-2008
%  Last modified: 06-09-2009
%
%--------------------------------------------------------------------------

function y = polyphase_components(x, pM, pN, type, IsOutFreq)

%  Get the polyphase components (periodicity pM, pN) of an image (or filter)
% 
%  INPUT
%    x: a 2-D image or filter
%
%    pM: the downsampling factor along the vertical axis
%
%    pN: the downsampling factor along the horizontal axis
%
%    type: the type of polyphase decomposition we use. There are two
%    different types in the multirate signal processing literature.
%
%    IsOutFreq: Give the output in the spatial (0) or frequency (1) domain
%
%  OUTPUT
%    y: a cell array containing the polyphase components


% sanity check
[M, N] = size(x);
if rem(M, pM) || rem(N, pN)
    error('pM and pN must be able to divide M and N, respectively.');
end

y = cell(pM, pN);

if strcmp(type, 's')
    % the default option
    shift = -1;
else
    shift = 1;
end

for m = 1 : pM
    for n = 1 : pN
        % shifting
        x2 = circshift(x, shift * [m-1, n-1]);
        % followed by downsampling
        y{m,n} = x2(1:pM:end,1:pN:end);
        
        if IsOutFreq
            y{m,n} = fftn(y{m,n});
        end
    end
end
