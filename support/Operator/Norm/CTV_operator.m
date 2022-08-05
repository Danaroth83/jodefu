%% Collaborative Total Variation Proximity Operator
%
% (Copyright) Daniele Picone
% ver. 3.0
%
% Description:
% This function returns the proximity operator  for a Collaborative Total
% Variation Norm of a 4D matrix x, whose dimensions are:
% [Spatial Columns,Spatial Rows,Bands,Derivatives]
% This implementation refers to the article [1]
%
% Usage:
% Op_handle=CTV_operator('Field',Value)
% Example:
% Op_handle=CTV_operator('norm','l221','order','dbp');
%
% Fields:
% 'norm:' the kind of norm that is used for the calculation (table 1 in [1])
% 'order': the order in which the pixels, bands and derivatives are considered in the norm
%        eg: 'dbp' means the order is derivative, bands, pixels
%
% Output:
% Op_handle: A struct whose fields are:
%     - direct: The function handle with a two inputs (x,lambda) for
%             calculating lambda*f(x), where f(x) is the requested norm
%     - prox: The proximal operator handle with two inputs (x,gamma)
%             relative to the gamma*f(x) function, where gamma>0 is a
%             generic scalar and f(x) is the selected norm
%     - proxconj: The proximal operator handle with three inputs
%             (x,gamma,lambda), relative to the function
%             g(x)= gamma [lambda f(x)]*, where gamma>0 and lambda>0 are
%             two positive scalars, f(.) is the requested norm and [.]* 
%             is the Fenchel conjugate operator (the convex conjugate of a 
%             function)
%
% Allowed values for the 'norm' field:
% 'l221': Norm l_221 (Vector Total Variation)
% 'l211': Norm l_211 (l_2 over the first dimension, l_1 over the rest)
% 'l2inf1': Norm l_2(inf)1 (l_2 over the first dimension, l_inf over the second, l_1 over the third)
% 'linfinf1': Norm l_inf over the first two dimension, l_1 over the third
% 'linf11': Norm l_inf over the first dimension, l_1 over the rest
% 'linf21': Norm l_inf over the first dimension, l_2 over the second, l_1 over the third
% 'S1l1': Shatten Norm over the first two dimensions, l_1 over the third
% 'S1Linf': Shatten Norm over the first two dimensions, l_inf over the third
% 'l1': Norm l_1 over all dimensions
% 'l2': Norm l_2 over all dimensions
% 'l2sq': Half norm l_2 squared over all dimensions
%
% References:
% [1] "Collaborative Total Variation: A General Framework for Vectorial TV
% Models" by Duran et al.

function Op=CTV_operator(varargin)

    norm_label='l221';
    norm_order=[];
    
    for ii=1:2:numel(varargin)
        pname=varargin{ii};
        pval=varargin{ii+1}; 
        if any(strcmpi(pname,{'norm','norm_label'}))  % Type of norm
            norm_label=pval; 
        elseif any(strcmpi(pname,{'order','norm_order'})) % Order of variables to consider in the norm
            norm_order=pval;
        end
    end
    
    norm_label=lower(norm_label);
    
    if any(strcmpi(norm_label,{'linf11b','linf11c'})) % *Special case as shortcut
        norm_label='linf11';
        norm_order='cdp';  % Columns (bands), derivative, pixel order
    end
    
    if isempty(norm_order)
         norm_order='dcp'; % Derivative, columns (bands), pixel order
    end

    norm_order=strrep(norm_order,'b','c');
    perm=zeros(1,3);
    for ii=1:3
        switch norm_order(ii)
            case 'p'
                perm(ii)=1;
            case 'c'
                perm(ii)=2;
            case 'd'
                perm(ii)=3;
            otherwise
                error('Norm order can''t be recognized');
        end
    end

    switch norm_label
        case 'l2sq'
            Op.prox     = @(x,gamma) 1/(1+gamma)*x;
            Op.proxconj = @(x,gamma,lambda) 1/(1+gamma*lambda)*x;
            Op.direct   = @(x,lambda) lambda*sum(x(:).^2)/2;
            Op.direconj = @(x,lambda) lambda*sum(x(:).^2)/2;
            Op.grad     = @(x,lambda) lambda*x;
            return; % This operator follows its own protocol
        case 's1l1'
            opprox  = @(x,gamma) prox_S1l1(x,gamma);
            opnorm  = @(x) norm_S1l1(x);
            oppconj = [];
            opnconj = [];
        case 'sinfl1'
            opprox  = @(x,gamma) prox_Sinfl1(x,gamma);
            opnorm  = @(x) norm_Sinfl1(x);
            oppconj = [];
            opnconj = [];
        case 'l2inf1'              
            opprox  = @(x,gamma) prox_l2inf1(x, gamma);
            oppconj = @(x,lambda) projball_l21inf(x, lambda);
            opnorm  = @(x) sum(sqrt(max(sum(x.^2,1),[],2)),3); % Sqrt external to save computation time
            opnconj = @(x) max(sum(sqrt(sum(x.^2,1)),2),[],3);
        case 'linf21'
            opprox  = @(x,gamma) ipermute(prox_l2inf1(permute(x,[2,1,3:ndims(x)]),gamma),[2,1,3:ndims(x)]);
            oppconj = @(x,lambda) ipermute(projball_l21inf(permute(x,[2,1,3:ndims(x)]),lambda),[2,1,3:ndims(x)]);
            opnorm  = @(x) sum(sqrt(sum(max(abs(x),[],1).^2,2)),3);
            opnconj = @(x) sqrt(max(sum(sum(abs(x),1).^2,2),[],3));  % Sqrt external to save computation time
        case 'linf11' 
            opprox  = @(x,gamma) prox_linf1(x,gamma);
            oppconj = @(x,lambda) projball_l1inf(x,lambda);
            opnorm  = @(x) sum(max(abs(x),[],1),2:3);
            opnconj = @(x) max(sum(abs(x),1),[],2:3);
        case 'linfinf1'
            opprox  = @(x,gamma) reshape(prox_linf1(conc(x),gamma),size(x));
            oppconj = @(x,lambda) reshape(projball_l1inf(conc(x),lambda),size(x));
            opnorm  = @(x) sum(max(abs(x),[],1:2),3);
            opnconj = @(x) max(sum(abs(x),1),[],2:3);
        case 'l111'
            opprox  = @(x,gamma) prox_l1(x,gamma);
            oppconj = @(x,lambda) projball_linf(x,lambda);
            opnorm  = @(x) sum(abs(x),1:3);
            opnconj = @(x) max(abs(x),[],1:3);
        case 'l1'
            opprox  = @(x,gamma) prox_l1(x,gamma);
            oppconj = @(x,lambda) projball_linf(x,lambda);
            opnorm  = @(x) sum(abs(x(:)));
            opnconj = @(x) max(abs(x(:)));
        case 'l112'
            opprox  = @(x,gamma) ipermute(prox_l21(permute(x,[3,1,2,4:ndims(x)]),gamma),[3,1,2,4:ndims(x)]);
            oppconj = @(x,lambda) ipermute(projball_l2inf(permute(x,[3,1,2,4:ndims(x)]),lambda),[3,1,2,4:ndims(x)]);
            opnorm  = @(x) sqrt(sum(sum(abs(x),1:2).^2,3));    
            opnconj = @(x) sqrt(sum(max(abs(x),[],1:2).^2,3));
        case 'l121'
            opprox  = @(x,gamma) ipermute(prox_l21(permute(x,[2,1,3:ndims(x)]),gamma),[2,1,3:ndims(x)]);
            oppconj = @(x,lambda) ipermute(projball_l2inf(permute(x,[2,1,3:ndims(x)]),lambda),[2,1,3:ndims(x)]);
            opnorm  = @(x) sum(sqrt(sum(sum(abs(x),1).^2,2)),3);    
            opnconj = @(x) sqrt(max(sum(max(abs(x),[],1).^2,2),[],3)); %Sqrt external to save computation time
        case 'l211'
            opprox  = @(x,gamma) prox_l21(x,gamma);
            oppconj = @(x,lambda) projball_l2inf(x,lambda);
            opnorm  = @(x) sum(sqrt(sum(x.^2,1)),2:3);    
            opnconj = @(x) sqrt(max(sum(x.^2,1),[],2:3)); %Sqrt external to save computation time
        case 'l221'
            opprox  = @(x,gamma) reshape(prox_l21(conc(x),gamma),size(x));
            oppconj = @(x,lambda) reshape(projball_l2inf(conc(x),lambda),size(x));
            opnorm  = @(x) sum(sqrt(sum(x.^2,1:2)),3);
            opnconj = @(x) sqrt(max(sum(x.^2,1:2),[],3));
        case 'l222'
            opprox  = @(x,gamma) reshape(prox_l21(conc(conc(x)),gamma),size(x));
            oppconj = @(x,lambda) reshape(projball_l2inf(conc(conc(x)),lambda),size(x));
            opnorm  = @(x) sqrt(sum(x.^2,1:3));
            opnconj = @(x) sqrt(sum(x.^2,1:3));
        case 'l2'
            opprox  = @(x,gamma) reshape(prox_l21(x(:),gamma),size(x));
            oppconj = @(x,lambda) reshape(projball_l2inf(x(:),lambda),size(x));
            opnorm  = @(x) sqrt(sum(x(:).^2));
            opnconj = @(x) sqrt(sum(x(:).^2));
        otherwise
            error('Requested norm is not available');
    end

    Op.prox   = @(x,gamma) reshape(ipermute(opprox(permute(conc(x),perm),gamma),perm),size(x));
    Op.direct = @(x,lambda) lambda*opnorm(permute(conc(x),perm));
    if isempty(oppconj)
        Op.proxconj = @(x,gamma,lambda) x-gamma*Op.prox(x/gamma,lambda/gamma);
    else
        Op.proxconj = @(x,~,lambda) reshape(ipermute(oppconj(permute(conc(x),perm),lambda),perm),size(x));
    end
    if ~isempty(opnconj)
        Op.direconj = @(x,lambda) indicator(opnconj(permute(conc(x),perm))<=lambda);
    end
    
end


%% SUPPORTING FUNCTIONS

% Concatenation of the first two dimensions into one
function x=conc(x)
    sizes=size(x); if length(sizes)==2,sizes(3)=1; end
    x=reshape(x,[sizes(1)*sizes(2),sizes(3:end)]);
end

% Indicator function
function y=indicator(cond)
    if cond, y=0; else, y=Inf; end  
end

% Norm S_1-l_1
function out=norm_S1l1(in)
    sizes=size(in);
    in=reshape(in,size(in,1),size(in,2),[]);
    out=zeros([1,1,size(in,3)]);
    for ii=1:size(in,3)
        out(ii)=sum(svd(in(:,:,ii)));
    end
    out=sum(reshape(out,[1,1,sizes(3:end)]),3);
end

% Norm S_{inf}-l_1
function out=norm_Sinfl1(in)
    sizes=size(in);
    in=reshape(in,size(in,1),size(in,2),[]);
    out=zeros([1,1,size(in,3)]);
    for ii=1:size(in,3)
        out(ii)=max(svd(in(:,:,ii)));
    end
    out=sum(reshape(out,[1,1,sizes(3:end)]),3);
end
