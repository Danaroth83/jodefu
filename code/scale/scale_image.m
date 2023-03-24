function [I_out,edge_out,time]=scale_image(I_in,ratio,flag_upscale,flag_downscale,varargin)

% Usage for upscaling by a scale factor "ratio":
%       I_out=scale_image(I_in,ratio,flag_upscale,flag_downscale,'edge',edgesize)
%       for downscaling by a scale factor "ratio":
%       I_out=scale_image(I_in,1/ratio,flag_upscale,flag_downscale,'GNyq',GNyq)
% Resizes images according to flag_scale option
% If flag_upscale=1 or flag_downscale=1, it uses imresize
% For upsampling, if flag_upscale=0 it uses an interp. kernel generated with fir1
%                 if flag_upscale=2 it uses a Lagrange polynomial interpolator
%                 
% For downsampling, if flag_downscale=0 it uses an MTF-matched filter and decimation
%                   if flag_downscale=2, the MTF-matching is done with a
%                   different set of gains at Nyquist frequency (GNyq)
% edgesize are the sizes of the edges of the image are in the format
% [up,left;down,right]

[L1,L2,Nb]=size(I_in);
if numel(ratio)==1, ratio=[ratio,ratio]; end
if isequal(size(ratio),[2,1]), ratio=ratio.'; end
if nargin<=2, flag_upscale=2; end
if nargin<=3, if Nb==1, flag_downscale=1; else, flag_downscale=0; end; end

if Nb==1, GNyq=0.15; else, GNyq=0.29.*ones(1,Nb); end
GNyqMS=GNyq;
edge=0;
flag_removeedge=1;
flag_disalign=1;

if nargin<=1 || (nargin>=5 && mod(nargin,2)~=0)
    error('Amount of inputs is incorrect');
end

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    
    if strcmpi(pname,'edge')
        edge=pval;
    elseif strcmpi(pname,'GNyq')
        GNyq=pval;
        if length(GNyq)~=Nb
            error('Nyquist gain must be a vector with length equal to the number of bands');
        end
    elseif strcmpi(pname,'GNyq2')
        GNyqMS=pval;
        GNyqMS(isnan(GNyqMS))=mean(GNyqMS(~isnan(GNyqMS)));
    elseif strcmpi(pname,'removeedge')
        flag_removeedge=pval;
        if flag_removeedge~=0  && flag_removeedge~=1
            error('remove edge flag can only be 1 or 0');
        end
    elseif strcmpi(pname,'disalign')
        flag_disalign=pval;
        if flag_disalign~=0  && flag_disalign~=1
            error('disalign flag can only be 1 or 0');
        end
    else
        error('Unrecognized input field');
    end
end

if (ratio(1)<1 && ratio(2)>1) || (ratio(1)>1 && ratio(2)<1)
    [I_in,edge,time1]=scale_image(I_in,[ratio(1),1],flag_upscale,flag_downscale,...
        'edge',edge,'removeedge',flag_removedge,'GNyq',GNyq,'GNyq2',GNyqMS,...
        'disalign',flag_disalign);
    [I_out,edge_out,time2]=scale_image(I_in,[1,ratio(2)],flag_upscale,flag_downscale,...
        'edge',edge,'removeedge',flag_removedge,'GNyq',GNyq,'GNyq2',GNyqMS,...
        'disalign',flag_disalign);
    time=time1+time2;
    return;
end

size_edge=size(edge);
edge=repmat(edge,3-size_edge);

tol=1e-10;
if any(rem(L1*ratio,1)>=tol) || any(rem(L2*ratio,1)>=tol)
    error('Output image sizes are not integer');
end

%Removing extra edge
edge_old=edge;
for ii=1:4
    [i1,i2]=ind2sub([2,2],ii);    
    edge(i1,i2)=find(rem((0:edge_old(i1,i2))*ratio(i2),1)<tol,1,'last')-1;
end
edge_rem=edge_old-edge;
I_in=I_in(edge_rem(1,1)+1:end-edge_rem(2,1),edge_rem(1,2)+1:end-edge_rem(2,2),:);

if nargout==3, t1=tic; end

if isequal(ratio,[1,1])
    I_out=I_in;
elseif ratio(1)>=1 && ratio(2)>=1
    if flag_upscale==0
        if flag_disalign==1
            I_out=interpCascade_offcenter(I_in,ratio,0);
        else
            I_out=interpCascade_mod(I_in,ratio,0);
        end
    elseif flag_upscale==1
        if flag_disalign==1 %unworking
            error('This combination of upscaling methods is currently unavailable');
            % decimation=[1,1];
            % if rem(ratio(1),2)==0, a=2*L1*ratio(1); decimation(1)=2; else a=L1*ratio(1); end
            % if rem(ratio(2),2)==0, b=2*L2*ratio(2); decimation(2)=2; else b=L2*ratio(2); end
            % I_out=imresize(I_in,[a,b]);
            % I_out=I_out(decimation(1):decimation(1):end,decimation(2):decimation(2):end,:);
        else
            I_out=imresize(I_in,round([L1,L2].*ratio));
        end
    elseif flag_upscale==2
        if flag_disalign==1
            I_out=interpCascadeAlparone_offcenter(I_in,ratio,0);
        else
            I_out=interpCascadeAlparone_mod(I_in,ratio,0);
        end
    end
else
    if flag_downscale==0 || flag_downscale==1
        I_out=resize_images_MS(I_in,1./ratio,GNyq,flag_downscale,1,flag_disalign);
    elseif flag_downscale==2
        I_out=resize_images_MS(I_in,1./ratio,GNyqMS,0,1,flag_disalign);
    end
end

if nargout==3, time=toc(t1); end

rer=round(edge.*repmat(ratio,[2,1]));
if flag_removeedge==1
    I_out=I_out(rer(1,1)+1:end-rer(2,1),rer(1,2)+1:end-rer(2,2),:);
end

if nargout>=2
    if flag_removeedge==1
        edge_out=zeros(size_edge);
    else
        if isequal(size_edge,[1,1])
            edge_out=min(rer(:));
        elseif isequal(size_edge,[1,2])
            edge_out=min(rer,[],1);
        elseif isequal(size_edge,[2,1])
            edge_out=min(rer,[],2);
        else
            edge_out=rer;
        end
        edge_rem=rer-repmat(edge_out,3-size_edge);
        I_out=I_out(edge_rem(1,1)+1:end-edge_rem(2,1),edge_rem(1,2)+1:end-edge_rem(2,2),:);
    end
end
    