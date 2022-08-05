function Op=load_degmaskoperator(varargin)

% Author:
% Daniele Picone
% v2.0
%
% Description:
% This function generates the operator defining the direct model of the
% transfer function for operation involving image fusion (pansharpening)
% and mosaicking.
% This function assumes that the input of the direct model is a fused 
% image (sizes: L1 x L2 x Nb) of two sources and the output is a mosaicked
% generated from a stack of sizes L1 x L2 x (Nb1+Nb), where the first Nb
% bands are from the high resolution image and the remaining Nb bands are
% an interpolated version of the low resolution image.
% 
% Usage:
% Eg of usage for image fusion without mosaicking:
% Op=load_degmaskoperator('lpfilter',ones(4)/16,'spectralweights',ones(4,1)/4,'sizes',[512,512,4]);
% Eg of usage for demosaicking without fusion:
% Op=load_degmaskoperator('mask',randn([512,512,4]));
% Eg of usage for image fusion with mosaicking:
% Op=load_degmaskoperator('lpfilter',ones(4)/16,'spectralweights',ones(4,1)/4,'mask',randn([512,512,5]));


% Output:
% Op: A struct for the masking operator whose fields are:
%   - direct: function handle for the direct model operator
%   - adjoint: function handle of the adjoint operator of the above

% This function takes as input a variable set of inputs:
% mask: mask necessary to generate the compressed image  (sizes: L1 x L2 x (Nb1+Nb),
%       default: [] = no mask applied, sizes field necessary)
% lpfilter: degradation spatial filter for the generating the low
%           resolution image from the ideal image
%           (sizes: any x any x Nb, default: [] = no degradation)
% spectralweights: weights necessary to generate a simulated HR image from
%                   the ideal one (sizes: Nb x Nb1)
% shift: shifting of the masked image (sizes: (Nb1+Nb) x 2, default: [] = no shifting)
% sizes: sizes of the fused image (needed format: [L1,L2,Nb])

current_folder=fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'mask'));
addpath(fullfile(current_folder,'operator'));

sizes=[];
mask=[];
shift=[];
lpfilter=[];
spectralweights=[];
flag_fft=[];
band_HR=0;
padding=[];

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if any(strcmpi(pname,{'cfa','mask'}))
        if isnumeric(pval), mask=pval; shift=[]; band_HR=0; end
        if isfield(pval,'data'), mask=pval.data; end
        if isfield(pval,'shift'), shift=pval.shift; end
        if isfield(pval,'band_HR'), band_HR=pval.band_HR; end
        sizes=size(mask); sizes(3)=size(mask,3)-band_HR; 
    elseif strcmpi(pname,{'lpfilter'})  % Degradation spatial filter
            lpfilter=pval;
    elseif strcmpi(pname,{'GNyq'})  % GNyq: Gain at Nyquist frequency of spectral sensor
            error('This function requires the full expression of the degradation filter');
    elseif strcmpi(pname,'spectralweights') % Percentage pectral coverage of MS over PAN
        spectralweights=pval;
        if size(spectralweights,1)<size(spectralweights,2), spectralweights=spectralweights.'; end
    elseif any(strcmpi(pname,{'flag_fft','fft'})) % Percentage pectral coverage of MS over PAN
        flag_fft=pval;
    end
end

%% Check for conditions
if isempty(sizes), error('The field sizes is required'); end
if isempty(padding), padding=find_zeropad(mask,shift); end

%% Operator for the mask
if isempty(mask)
    op_dirmask = @(x) x;
    op_adjmask = @(y) y;
    norm_mask  = 1;
elseif isempty(shift) || all(shift(:)==0)
    op_dirmask = @(x) sum(x.*mask,3);
    op_adjmask = @(y) y.*mask;
    norm_mask  = sqrt(max(sum(mask.^2,3),[],1:2));
else
    
    % padding=[-min(min(shift(:,1)),0),-min(min(shift(:,2)),0),...
    %               max(max(shift(:,1)),0),max(max(shift(:,2)),0)];
    % padding=find_zeropad(mask,shift);
    op_dirmask = @(x) opdir_mosaicshift(x,mask,shift,padding);
    op_adjmask = @(y) opadj_mosaicshift(y,mask,shift,padding);
    maskshifted= opdir_mosaicshift(ones(size(mask)),mask,shift,padding);
    norm_mask  = sqrt(max(sum(maskshifted.^2,3),[],1:2));
end

%% Operator for the degradation
if isempty(lpfilter) || numel(lpfilter)==1
    op_dirdeg = @(x) x;
    op_adjdeg = @(y) y;
    norm_deg  = 1;
else
    
    if isempty(spectralweights) && band_HR>0, spectralweights=ones(size(lpfilter,3),band_HR)/size(lpfilter,3); end
    
    norm_lpfilter = max_svd_conv(lpfilter);
    norm_spectralweights = sqrt(sum(spectralweights.^2,1:2));
    norm_deg = sqrt(norm_lpfilter.^2+norm_spectralweights.^2);
    
    if isempty(flag_fft), if max(size(lpfilter,1),size(lpfilter,2))>32, flag_fft=1; else, flag_fft=0; end; end
    if flag_fft==1
        lpfilter_pad=padarray(lpfilter,[ceil((sizes(1)-size(lpfilter,1))/2),ceil((sizes(2)-size(lpfilter,2))/2)],0,'pre');
        lpfilter_pad=padarray(lpfilter_pad,[floor((sizes(1)-size(lpfilter,1))/2),floor((sizes(2)-size(lpfilter,2))/2)],0,'post');
        fft2_lpfilter=fft2(lpfilter_pad);
        conj_fft2_lpfilter=conj(fft2_lpfilter);
        op_dirdeg = @(x) opdir_fusionfft(x,fft2_lpfilter,spectralweights);
        op_adjdeg = @(y) opadj_fusionfft(y,conj_fft2_lpfilter,spectralweights);
    else
        op_dirdeg = @(x) opdir_fusion(x,lpfilter,spectralweights);
        op_adjdeg = @(y) opadj_fusion(y,lpfilter,spectralweights);
    end
end


%% Combine operators

Op.mosaic   = @(x) op_dirmask(x);
Op.demosaic = @(x) op_adjmask(x);
Op.direct   = @(x) op_dirmask(op_dirdeg(x));
Op.adjoint  = @(y) op_adjdeg(op_adjmask(y));
Op.norm     = norm_mask*norm_deg;

end


