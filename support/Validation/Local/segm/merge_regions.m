function Snew = merge_regions(S, I, minsize)
% suppress regions in a segmentation map having size smaller that minsize

Snew = S;

labels = unique(Snew);
nreg = length(labels);

regsize = zeros(1, nreg);
for i=1:nreg
    regsize(i) = nnz(Snew==labels(i));
end

labelsmallregs = labels(regsize < minsize);

while ~isempty(labelsmallregs)
    for i=1:length(labelsmallregs)
        curlab = labelsmallregs(i);
        
        % find pixels 8-connected to the border
        adjpixels = logical(imdilate(double(Snew==curlab), strel('square',3)) - double(Snew==curlab));
        adjpixellabels = Snew(adjpixels);
        
        adjlabels = unique(adjpixellabels);
        adjval = I(adjpixels);
        
        % merge the region to the adjacent one with closest mean on the
        % adjacent pixels
        adjmeans = zeros(1, length(adjlabels));
        for j=1:length(adjmeans)
            if nnz(~labelsmallregs==adjlabels(j))==0
                adjmeans(j) = mean(adjval(adjpixellabels==adjlabels(j)));
            else
                adjmeans(j) = Inf;
            end
        end
        if nnz(~(isinf(adjmeans)))>0
            regmean = mean(I(Snew==curlab));
            [~, idx] = min(abs(adjmeans-regmean*ones(size(adjmeans))));
            Snew(Snew==curlab) = adjlabels(idx);
        end
    end
    
    labels = unique(Snew);
    nreg = length(labels);
    
    regsize = zeros(1, nreg);
    for i=1:nreg
        regsize(i) = nnz(Snew==labels(i));
    end
    labelsmallregs = labels(regsize < minsize);
end

S = Snew;
for i=1:nreg
    Snew(S==labels(i)) = i;
end

