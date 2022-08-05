%% Collaborative Total Variation Proximity Operator
%
% Description:
% This function returns the proximity operator  for a Collaborative Total
% Variation Norm of a 3D matrix x, whose dimensions are [pixels,Bands, Derivatives]
% This implementation refers to the article [1]
%
% Usage:
% Op_handle=CTV_operator('Field',Value)
%
% Example:
% Op_handle=CTV_operator('norm','l221','order','dbp');
% x=double(imread('peppers.png'));
% x_tv=cat(4,[diff(x,1,1);zeros(1,size(x,2),size(x,3))],[diff(x,1,2) zeros(size(x,1),1,size(x,3))]);
% y=Op_handle(x_tv,0.1);
%
% Optional fields are the following:
% norm: the kind of norm that is used for the calculation (table 1 in [1])
% order: the order in which the pixels, bands and derivatives are considered in the norm
%        eg: 'dbp' means the order is derivative, bands, pixels
% size: sizes of the desired output
%
% Allowed fields for the norm:
% 'L221': Norm l_221 (Vector Total Variation)
% 'L211': Norm l_211 (l_2 over the first dimension, l_1 over the rest)
% 'L2inf1': Norm l_2(inf)1 (l_2 over the first dimension, l_inf over the second, l_1 over the third)
% 'Linfinf1': Norm l_inf over the first two dimension, l_1 over the third
% 'Linf11': Norm l_inf over the first dimension, l_1 over the rest
% 'Linf21': Norm l_inf over the first dimension, l_2 over the second, l_1 over the third
% 'S1L1': Shatten Norm over the first two dimensions, l_1 over the third
% 'S1Linf': Shatten Norm over the first two dimensions, l_inf over the third
%
% References:
% [1] "Collaborative Total Variation: A General Framework for Vectorial TV
% Models" by Duran et al.

function Op_handle=CTV_operator_old(varargin)

    norm_label='l221';
    norm_order=[];
    sizes=[];

    for ii=1:2:numel(varargin)
        pname=varargin{ii};
        pval=varargin{ii+1}; 
        if any(strcmpi(pname,{'norm','norm_label'}))  % Type of norm
            norm_label=pval; 
        elseif any(strcmpi(pname,{'order','norm_order'})) % Order of variables to consider in the norm
            norm_order=pval;
        elseif any(strcmpi(pname,{'size','size_img'}))    % Size of the image (only considered if input matrix is 3D)
            sizes=pval;
        end
    end
    
    norm_label=lower(norm_label);

    if isempty(norm_order)
        if any(strcmpi(norm_label,{'l221','l211','l111','l121','l2inf1','l1inf1'}))
            norm_order='dcp'; % Derivative, columns (bands), pixel order
        elseif any(strcmpi(norm_label,{'l212','l112'}))
            norm_order='dpc'; % Derivative, pixel, columns (bands) order
        elseif any(strcmpi(norm_label,{'linf11','linf21','linfinf1','s1l1','s2l1','sinfl1'}))
            norm_order='cdp'; % Columns (bands), derivative, pixel order
        end
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
        case 's1l1'
            opy= @(x,gamma) ipermute(S1L1(permute(x,perm),gamma),perm);
        case 's1linf'
            opy= @(x,gamma) ipermute(S1Linf(permute(x,perm),gamma),perm); 
        case 'l2inf1'              
            opy= @(x,gamma) l2inf1(x, gamma, perm);
        case 'l2inf1_extra'
            opy= @(x,gamma) l2inf1_alt(x, gamma, perm);
        case 'linf21'
            opy= @(x,gamma) linf21(x,gamma,perm);
        case 'linf11'
            opy = @(x,gamma) x - projball_l1_alt(x,gamma,perm);
        case 'linf11_extra' 
            opy = @(x,gamma) x - projball_l1(x,gamma,perm);
        case 'linfinf1'
            opy = @(x,gamma) linfinf1(y, gamma, perm);
        case 'l111'
            opy = @(x,gamma) x - x ./ max(abs(x)/gamma,1);
        case 'l211'
            opy = @(x,gamma) x - x ./ max(sqrt(sum(x.^2,perm(1)))/gamma,1);
        case 'l221'
            opy = @(x,gamma) x - x ./ max(sqrt(sum(x.^2,perm(1:2)))/gamma,1);
        case 'l222'
            opy = @(x,gamma) x - x ./ max(sqrt(sum(x.^2,perm(1:3)))/gamma,1);
        otherwise
            error('Requested norm is not available');
    end

    if isempty(sizes)
        Op_handle= @(x,gamma) reshape(opy(reshape(x,[],size(x,3),size(x,4)),gamma),size(x));
    else
        Op_handle= @(x,gamma) reshape(opy(reshape(x,[],sizes(3:end)),gamma),sizes);
    end

end


%% SUPPORTING FUNCTIONS

% Projection over a l1 ball of radius tau (over dimension dim)

function x = projball_l1(y, tau, dim) 
    tmp = abs(y);
    % if sum(tmp)<=tau, x = y; return; end
	lambda = max((cumsum(sort(tmp,dim(1),'descend'))-tau)./shiftdim((1:size(tmp,dim(1))),2-dim(1)),[],dim(1));
	x = y - max(min(y, lambda),-lambda);
end

function x = projball_l1_alt(y, tau, dim)
    tmp = abs(y);
	x= max(tmp-max(max((cumsum(sort(tmp,1,'descend'),1)-tau)./shiftdim((1:size(tmp,dim(1))),2-dim(1)),[],dim(1)),0),0).*sign(y);    
end

% Projection over a l12 ball of radius tau (over dimension [dim(1),dim(2)])

function x = projball_l12(y, tau, dim)
    tmp = sqrt(sum(y.^2,dim(1)));
	% if sum(tmp)<=tau, x = y; return; end
	lambda = max((cumsum(sort(tmp,dim(1),'descend'))-tau)./shiftdim((1:size(y,dim(2))),2-dim(2)),[],dim(1));
	tmp = max(tmp, lambda);
	x = y.*(tmp-lambda)./tmp;
end

% Norms l_(inf)21,  l_2(inf)1 l_(inf)(inf)1

function x = linf21(y, tau, dim)
    x=y-sign(y).*projball_l12(y,tau,dim);
    x(isnan(x))=0;
end

function x = l2inf1(y, tau, dim)
    tmp1 = sqrt(sum(y.^2,dim(1)));
    tmp2 = min(max(max((cumsum(sort(tmp1,dim(2),'descend'),2)-tau)./shiftdim((1:size(y,dim(2))),2-dim(2)),[],dim(2)),0)./tmp1,1);
    x = y .* tmp2;
end

function x = l2inf1_alt(y, tau, dim)
    tmp1=sqrt(sum(y.^2,dim(1)));
    v=projball_l1(tmp1,tau,dim(2));
    x=y.*max(tmp1-v,0)./tmp1;
    x(isnan(x))=0;
end

function x = linfinf1(y, tau, dim)
    y=permute(y,dim);
    [N1,N2,N3]=size(y);
    y=reshape(y,[],size(y,3));
    tmp=abs(y);
    x = y - max(tmp-max(max((cumsum(sort(tmp,1,'descend'),1)-tau)./(1:N1*N2)',[],1),0),0).*sign(y);
    x=ipermute(x,dim);
    x=reshape(x,[N1,N2,N3]);
end



