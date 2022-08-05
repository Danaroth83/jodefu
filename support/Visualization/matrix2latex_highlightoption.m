function highlight=matrix2latex_highlightoption(columnLabels)

% support function for highlighting best results in matrix
% takes columnLabels as input

highlight='';
highlight(1)='c';

for ii=1:numel(columnLabels)
    if strncmpi(columnLabels{ii},'D_',2) || strncmpi(columnLabels{ii},'SAM',3) || strncmpi(columnLabels{ii},'ERGAS',5) || strncmpi(columnLabels{ii},'MSE',3)
        highlight(ii+1)='s';
    else
        highlight(ii+1)='g';
    end
end
