function varargout = retrievesegmentation(prunedtree,varargin)

global mapwhed
global mapp

segmap = zeros(size(mapwhed));

PosNodes = [prunedtree.pruning];
node = prunedtree(logical(PosNodes));

for i=1:length(node)
    
    if isempty(node(i).nodeinfo.leaves)
        % the considered node is a leaf
        segmap(mapwhed==node(i).label) = node(i).label;
    else
        ind = ismember(mapwhed,node(i).nodeinfo.leaves);
        segmap(ind) = node(i).label;
    
    end
end

if nargout==1 % desired output = segmentation map
    
    varargout(1) = {segmap};

elseif nargout==0 % desired output = display of the segmentation map
    
    regionlabels = [node.label];
    borders = zeros(size(mapwhed));
    for i=1:length(regionlabels)
        reg = segmap==regionlabels(i);
        reg = bwmorph(reg,'remove');
        borders(mapp==0&reg) = 1;
    end
    borders = bwmorph(borders,'spur');
    borders = imcomplement(borders);
    
    if nargin == 1 % displays only the borders
        
        figure, imshow(borders)
        xlabel('corresponding segmentation of the pruned tree')
        
    elseif nargin == 2 % superimpose the borders to the original image
        
        im = varargin{1};
        if size(im,1)==size(borders,1)&&size(im,2)==size(borders,2)
            
            ind = borders~=0;
            indn = repmat(ind,[1 1 size(im,3)]);
            im = 1.2*im; im(im>1) = 1;
            im(indn==0) = 0;
            figure, imshow(im)
            
        else
            error('The dimensions of the input image don''t match the segmentation map')
        end
        
    else
        error('Wring number of input arguments. Please check')
    end
else
    error('Wrong number of output arguments. Please check')
end