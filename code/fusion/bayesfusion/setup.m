%  Set the path of necessary m file.
%  Date: April 2nd, 2015
%  QI WEI, University of Toulouse
clear mypath;   
st = pwd;
mypath = [ ...
    [st,';'], ...    
    [st,'\func_global;'],...
    [st,'\Quality_Indices;'],... 
    [st,'\func_Sparse']    
    ];
addpath(mypath);
clear mypath;


