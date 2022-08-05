function varargout=size_check(ratio,varargin)

% This function cuts a series of images, so that when downscaling the first
% image by a factor "ratio", resulting images have integer sizes.
% External edges inputted in the specific field will be ignored in this
% calculation and will be adjusted for the desired resize.
% Usage:
% [I_out1,I_out2,...,edge]=crop_borders([ratio1,ratio2,...],I_in1,I_in2,...,'edge',[edge1,edge2,...]);

tol=10^(-6);
idx_edge=0;
for ii=1:nargin-1
    if strncmpi(varargin{ii},'ed',2)
        idx_edge=ii;
    end
end

if ((nargout>nargin  || nargin<2) && idx_edge==0) || ((nargout>nargin-2 || nargin<4) && idx_edge~=0)
    error('Usage: [I_out1,I_out2,...,edge]=crop_borders(ratio,I_in1,I_in2,...,''edge'',[edge1,edge2,...])');
end

if idx_edge==0
    N=nargin;
    edges=zeros(1,N);
else
    N=nargin-3;
    edges=varargin{N+2};
end
if length(ratio)>N, error('too many inputted scaling ratios'); end
ratio=[ratio,ratio(1)*ones(1,N-length(ratio))];
input_images=varargin(1:N);

ysize=zeros(1,N);
xsize=zeros(1,N);
ratio_in=zeros(1,N);

for ii=1:N
    ysize(ii)=size(input_images{ii},1)-2*edges(ii);
    xsize(ii)=size(input_images{ii},2)-2*edges(ii);
    ratio_in(ii)=xsize(ii)/xsize(1); 
end
if any(ysize<=0) || any(xsize<=0), error('Image sizes are zero'); end

ratio_in2=ratio_in./ratio;
limit=min(xsize(1),ysize(1))/ratio_in(1);
multiplier=1;
while all(rem(ratio_in2*multiplier,1)<tol | rem(ratio_in2*multiplier,1)>1-tol)==0
    multiplier=multiplier+1;
    if multiplier>limit
        error('Some error occurred involving non integer scale ratios');
    end
end

ratio_multiplied=ratio_in*multiplier;
% ysize_new=floor(ysize./ratio_multiplied).*ratio_multiplied;
% xsize_new=floor(xsize./ratio_multiplied).*ratio_multiplied;
% edges_new=floor(edges./ratio_multiplied).*ratio_multiplied;
edges_new=zeros(size(edges)); ysize_new=zeros(size(ysize)); xsize_new=zeros(size(xsize));
for ii=1:numel(edges)
    ysize_new(ii)=find(rem(0:1:ysize(ii),ratio_multiplied(ii))<tol | rem(0:1:ysize(ii),ratio_multiplied(ii))>ratio_multiplied(ii)-tol,1,'last')-1;
    xsize_new(ii)=find(rem(0:1:xsize(ii),ratio_multiplied(ii))<tol | rem(0:1:xsize(ii),ratio_multiplied(ii))>ratio_multiplied(ii)-tol,1,'last')-1;
    edges_new(ii)=find(rem(0:1:edges(ii),ratio_multiplied(ii))<tol | rem(0:1:edges(ii),ratio_multiplied(ii))>ratio_multiplied(ii)-tol,1,'last')-1;
end
edges_cut=edges-edges_new;

varargout=cell(1,N);
for ii=1:N
    image_output=input_images{ii};
    varargout{ii}=image_output(edges_cut(ii)+(1:ysize_new(ii)+2*edges_new(ii)),...
        edges_cut(ii)+(1:xsize_new(ii)+2*edges_new(ii)),:);
end
varargout{N+1}=edges_new;

end

    
    