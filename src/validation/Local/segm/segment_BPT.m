function [S, T] = segment_BPT(img, val)

if nargin<2
    val = 1000;
end

global tree


if isa(img, 'struct')
    tree = img;
    % prunedtree = pruneBPTheight(5);
    % prunedtree = pruneBPTnbiterations(2850);
    prunedtree = pruneBPTnbregions(val);
    Stmp = retrievesegmentation(prunedtree,img);
    
    S = zeros(size(Stmp));
    val = unique(Stmp);
    for i=1:length(val)
        S(Stmp==val(i)) = i;
    end
else
%       keyboard
    addpath('segm/BPT/codeBPT-v1.3')
    [m n p] = size(img);
    global I
    I = reshape(img,m*n,p);
%     keyboard
    global mapp
    mapp = grad(img,'supremum');
    global mapwhed
    mapwhed = whed(mapp);
%     mapwhed = whed(mapp,I);
%     mapwhed=reshape(1:m*n,[m,n]);
    
%   mapwhed=round(img-min(img(:))+1);

    %% creation of the tree
    tic
    % global tree 
    tree = initstructarray(@R_mean,@O_SAM);
    % tree = initstructarray(@R_mean,@O_diff);
    toc
    
    tic
    updatestructarray(@O_SAM,@merging_mean,@priority_size);
    % updatestructarray(@O_diff,@merging_mean,@priority_size);
    toc
    
    completestructarray;
    
    
    fprintf('BPT with %d nodes\n', length(tree));
    
    %% pruning of the tree
%     keyboard
    % prunedtree = pruneBPTheight(5);
    % prunedtree = pruneBPTnbiterations(2850);
    prunedtree = pruneBPTnbregions(val);
    Stmp = retrievesegmentation(prunedtree,img);
    
    S = zeros(size(Stmp));
    val = unique(Stmp);
    for i=1:length(val)
        S(Stmp==val(i)) = i;
    end
    
    if nargout > 1
        T = tree;
    end
    
end