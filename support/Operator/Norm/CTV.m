function y=CTV(x,varargin)

% Calculates the proximity operator Collaborative Total Variation Norm of 
% a 4D matrix x, whose dimensions are [Vertical pixels, Horizontal pixels,
% Bands, Derivatives]
% This implementation refers to the article:
% [1] "Collaborative Total Variation: A General Framework for Vectorial TV
% Models" by Duran et al.
% Example of usage:
% y=CTV(x,'gamma',0.1,'norm','l221','order','dbp')
% Optional fields are the following:
% gamma: the CTV norm is multiplied by gamma  (default: 1)
% norm: the kind of norm that is used for the calculation (table 1 in [1])
% order: the order in which the pixels, bands and derivatives are considered in the norm
%        eg: 'dbp' means the order is derivative, bands, pixels
% size_img: needed if the image is passed as 3D to reconstruct horizontal and vertical pixels


gamma=1;
norm_label='l221';
norm_order=[];
size_img=[];

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if strcmpi(pname,'gamma')                          % Gamma parameter
        gamma=pval;   
    elseif any(strcmpi(pname,{'norm','norm_label'}))  % Type of norm
        norm_label=pval; 
    elseif any(strcmpi(pname,{'order','norm_order'})) % Order of variables to consider in the norm
        norm_order=pval;
    elseif any(strcmpi(pname,{'size','size_img'}))    % Size of the image (only considered if input matrix is 3D)
        size_img=pval;
    end
end

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

dim=length(size(x));

if dim==4
    size_img=[size(x,1),size(x,2)];
    x=reshape(x,[prod(size_img),size(x,3),size(x,4)]);
end

x=permute(x,perm);
[N1,N2,N3]=size(x);


switch norm_label
    case 's1l1'
        y=S1L1_new(x,gamma);

    case 's1linf'
        y=S1Linf_new(x,gamma); 
       
    case 'l2inf1'        
        tmp1 = sqrt(sum(x.^2,1));
        % projection of tmp1 on the l1-norm ball of radius divsigma :
        % max(bsxfun(@minus,tmp1,max(max(bsxfun(@rdivide,cumsum(...
        %	sort(tmp1,1,'descend'),1)-divsigma,(1:size(tmp1,1))'),[],1),0)),0);

        % tmp2 = min(bsxfun(@rdivide,max(max(bsxfun(@rdivide,cumsum(...
        %    sort(tmp1,2,'descend'),2)-gamma,(1:N2)),[],2),0),tmp1),1);
        tmp2 = min(max(max((cumsum(sort(tmp1,2,'descend'),2)-gamma)./(1:N2),[],2),0)./tmp1,1);
        y = x .* tmp2;
       
    case 'l2inf1_extra'
        tmp1=sqrt(sum(x.^2,1));
        v=shiftdim(projball_l1(squeeze(tmp1),gamma),-1);
        y=x.*bsxfun(@rdivide,max(tmp1-v,0),tmp1);
        y(isnan(y))=0;
        
    case 'linf21'
        % y=x-sign(x).*projball_l12(abs(x),gamma);
        y=x-projball_l12(x,gamma);
        y(isnan(y))=0;
       
    case 'linf11'
        tmp=abs(x);
        % y = x - max(bsxfun(@minus,tmp,max(max(bsxfun(@rdivide,cumsum(sort(...
		% tmp,1,'descend'),1)-gamma,(1:N1)'),[],1),0)),0).*sign(x);
        y = x - max(tmp-max(max((cumsum(sort(...
		 tmp,1,'descend'),1)-gamma)./(1:N1)',[],1),0),0).*sign(x);
     
    case 'linf11_extra'
        y = x - projball_l1(x,gamma);
    
    case 'linfinf1'
        x=reshape(x,[N1*N2,N3]);
        y = x - max(bsxfun(@minus,abs(x),max(max(bsxfun(@rdivide,cumsum(sort(...
		abs(x),1,'descend'),1)-gamma,(1:N1*N2)'),[],1),0)),0).*sign(x);
        y=reshape(y,[N1,N2,N3]);
        
    case 'l111'
        % y = x - bsxfun(@rdivide, x, max(abs(x)/gamma,1));
        y = x - x ./ max(abs(x)/gamma,1);
    case 'l211'
        % y = x - bsxfun(@rdivide, x, max(sqrt(sum(x.^2,1))/gamma,1));
        y = x - x ./ max(sqrt(sum(x.^2,1))/gamma,1);
    case 'l221'
        % y = x - bsxfun(@rdivide, x, max(sqrt(sum(sum(x.^2,1),2))/gamma,1));
        y = x - x ./ max(sqrt(sum(x.^2,[1,2]))/gamma,1);
end

y=ipermute(y,perm);

if ~isempty(size_img)
   y=reshape(y,[size_img(1),size_img(2),size(y,2),size(y,3)]);
end

end



function x = projball_l1(y, tau)
    dim=length(size(y));
    N1=size(y,1);
    tmp = abs(y);
	% if sum(tmp)<=tau, x = y; return; end
	lambda = repmat(max((cumsum(sort(tmp,1,'descend'))-tau)./(1:N1)',[],1),[N1,ones(1,dim-1)]);
	x = y - max(min(y, lambda),-lambda);
end

function x = projball_l12(y, tau)
	dim=length(size(y));
    N2=size(y,2);
    tmp = squeeze(sqrt(sum(y.^2,1)));
	% if sum(tmp)<=tau, x = y; return; end
	lambda = repmat(max((cumsum(sort(tmp,1,'descend'))-tau)./(1:N2)',[],1),[N2,ones(1,dim-2)]);
	tmp = max(tmp, lambda);
	x = bsxfun(@times, y, shiftdim((tmp-lambda)./tmp,-1));
end
