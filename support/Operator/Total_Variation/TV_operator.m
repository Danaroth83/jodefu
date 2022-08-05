%% Choice of Total Variation Operator
%
% Author:
% Daniele Picone
%
% Description:
% This function allows to select one operator function that allows to
% perform one of the Total Variation variants described at the end
%
% Usage:
% Op=TV_operator(choice)
% Example: Op=TV_operator('TV');
%
% Input:
% - choice: A string describing the desired operator,among the choices below:
%       'TV'  : Classical Total Variation
%       'TVu' : Upwind Total Variation
%       'TVsn': Shannon Total Variation (where 'n' is the oversampling parameter, eg:'TVs2')
%       'TVh' : Hessian operator (double gradient)
%       'TVg' : Gradient applicable  to a gradient operator (to create a Hessian)
%       'TVp' : New Total Variation proposed by Condat (takes a TV as input)
%       'TVm' : Variant of the above
%
% Output:
% - Op: A struct defining the desired operator f(x) with the following fields:
%       direct  : The function handle for the direct operation f(x)
%       adjoint : The function handle of the adjoint of the above 
%       norm    : The operator norm (non squared)


function Op=TV_operator(choice)
    
    if nargin<=0 || isempty(choice), choice='TV'; end

    if any(strcmpi(choice,{'c','d','TV','Gradient','Total_Variation'}))

        Op.direct  = @(x) TV_direct(x);
        Op.adjoint = @(x) TV_adjoint(x);
        Op.norm    = sqrt(8);
        
    elseif any(strcmpi(choice,{'g','TVg','TGV','Total_Generalized_Variation'}))    
    
        Op.direct  = @(x) Hessian_direct(x);
        Op.adjoint = @(x) Hessian_adjoint(x);
        Op.norm    = sqrt(8);
        
    elseif any(strcmpi(choice,{'u','TVu','uTV','Upwind_Total_Variation'}))    
    
        Op.direct  = @(x) TVu_direct(x);
        Op.adjoint = @(x) TVu_adjoint(x);
        Op.norm    = 4;
        
    elseif any(strcmpi(choice,{'h','TVh','Hessian','Hess'}))
        
        Op.direct  = @(x) Hessian_direct(TV_direct(x));
        Op.adjoint = @(x) TV_adjoint(Hessian_adjoint(x));
        Op.norm    = 8;
    
    elseif any(strcmpi(choice,{'p','TVp','newTV'}))
        
        Op.direct  = @(x) TVp_direct(x);
        Op.adjoint = @(x) TVp_adjoint(x);
        Op.norm    = 1;  % Wrong, double check
        
    elseif any(strcmpi(choice,{'m','TVm','TVpm','newTVm'}))
        
        Op.direct  = @(x) TVpm_direct(x);
        Op.adjoint = @(x) TVpm_adjoint(x);
        Op.norm    = 1;  % Wrong, double check
    
    elseif (strncmpi(choice,'s',1)) || (strncmpi(choice,'TVs',3))
        
        m=1; tail=0;
        while ~isnan(m), n=m; m=str2double(choice(end-tail:end)); tail=tail+1; end
        
        Op.direct  = @(x) TVs_direct(x,n);
        Op.adjoint = @(x) TVs_adjoint(x,n);
        Op.norm    = sqrt(2)*pi/n;  % According to Moisan code


    else
        error('Choice is unknown');     
    end
end