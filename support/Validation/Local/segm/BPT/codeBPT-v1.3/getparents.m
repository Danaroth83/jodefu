function branch = getparents(node,str)

global tree

if strcmp(str,'index')
    h = tree(node).nodeinfo.height;
    branch = zeros(1,h-1);
elseif strcmp(str,'structure')
    h = node.nodeinfo.height;
    branch(h-1) = struct;
else
    error('wrong input string parameter. Please check')
end