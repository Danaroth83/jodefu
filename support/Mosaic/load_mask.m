%% Mask Loading
%
% Description:
% This function loads a mask from a database of masks used in the
% literature. The mask also allow to be generated for multimodal signals,
% such as an High Resolution (HR) and and upsampled Low Resolution (LR); in
% this case the two masks are stacked together over the third dimension
% (HR first bands, LR remaining bands) 
%
% Usage:
% mask=load_mask('Field',Value);
% Example:
% mask=load_mask('sizes',[512,512,4],'mask','random','band_HR',1,'ratio',2,'alpha',0.5);
% 
% Input Fields:
% 'sizes': The sizes of the desired mask (format: [L1,L2,Nb], where L1 and L2
% are the vertical and the horizontal sizes, Nb is the amount of bands, default: [512,512,1])
% 'type': Generate a mask from a given database (default: 'ICIP2018')
% 'ratio': Scale ratio between HR and LR image
% 'bands': Number of bands assigned to the HR image
% 'alpha': Factor deciding percentage of mask pixels assigned to the HR
%          over the LR (default: 'same', special: 'same' assigns the same
%          value to the bands of the HR image and LR); valid for 'weighted'
%          and 'random' masks
% 
% Output:
% mask: A cell whose fields are
%    data: The requested mask (sizes: [L1,L2,Nb])
%    shift: shift of the input pixels before masking (sizes: [Nb,2])
%    band_HR: The amount of bands assigned to the PAN
%    padding: The amount of padding of the final masked image (format: [left, up, right, down])
%    mosaic:  Function handle to generate the mosaicked image.
%        - Usage: I_out=mask.mosaic(I_in,mask), where I_in is a struct
%        whose fields HR and LR contain the HR and LR image, respectively 
%    demosaic: Function handle to separate the two masks into the one for
%        HR and LR, in the fields mask_HR and mask_LR respectively.
%        Also returns the adjoint masking operator of the HR and LR,
%        in the fields image_HR and image_LR, respectively
%        - Usage: out=mask.demosaic(I_mosaic,mask); where I_mosaic is the mosaicked image
%
% Available mask types:
% 'random': Each pixel is randomly assigned to each available band
% 'weighted': Non binary mask; a random weight is assigned to each band
%             according to the Dirichlet distribution
% 'ICIP2018': One pixel every (scale ratio)^2 are assigned to the LR image
% 'Kodakver1','Kodakver2','Kodakver3': Standard periodic patented masks 
%       by Kodak
% 'BinaryTreeU','BinaryTreeD','BinaryTreeS': Binary mask splitting with 
%       uniform, dominant and hierarchical assignation [1]
% 'Bayer','Bayermod','QuadBayer': Bayer based patterns
% 'Yamanaka','Lukac': Lukac [2]/Yamanaka [3] pattern mask
% 'CASSI': Stack of random binary masks, shifted one pixel to the right for
%          each acquisition [4]
%
% References:
% [1] Miao L., Qi H., Ramanath R., Member, and Snyder W.E., "Binary 
% Tree-based Generic Demosaicking Algorithm for Multispectral Filter 
% Arrays", IEEE TIP, 15(11), 3550-3558.
% [2] Lukac R. and Plataniotis K.N., "Color filter arrays: design and 
% performance analysis", 
% [3] Yamanaka, S., "Solid state color camera", 1977
% [4] Arce G. R., Brady D. J., Carin L., Arguello H. and Kittle, D. S., 
% "Compressive coded aperture spectral imaging: An introduction.",
% IEEE SPM, 31(1), 105-115.


function mask=load_mask(varargin)

current_folder=fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'Mask'));
addpath(fullfile(current_folder,'Operator'));

%% Parsing inputs

sizes=[512,512,1];
start_pos=[];
masktype_PAN=[];
masktype_MS=[];
type=[];
ratio=[1,1];
flag_PAN=0;
band_HR=[];
alpha='same';
flag_chooseoriginalMS=1;

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1}; 
    if any(strcmpi(pname,{'sizes'}))                    % sizes of the final mask for the LR image
        if length(pval)==1, sizes=[pval,pval,1]; end
        if length(pval)==2, sizes=[pval,pval,1]; end
        if length(pval)>2, sizes=pval; end
    elseif any(strcmpi(pname,{'mask','type'}))          % Mask type
        type=pval;
    elseif any(strcmpi(pname,{'mask_LR','mask_MS','masktype_MS','masktype_LR'}))    % Mask type for the LR image
        masktype_MS=pval;
    elseif any(strcmpi(pname,{'mask_HR','mask_PAN','masktype_PAN','masktype_LR'}))   % Mask type for the HR image
        masktype_PAN=pval;
    elseif any(strcmpi(pname,{'band_HR','band_PAN'}))   % Amount of bands for the HR image
        band_HR=pval;
    elseif any(strcmpi(pname,{'scale','ratio'}))        % Scale ratio between HR and LR image
        if length(pval)==1, ratio=[pval,pval]; end 
        if length(pval)>1, ratio=pval; end
    elseif any(strcmpi(pname,{'alpha'}))                % alpha is the fraction of pixels assigned to the HR image
        alpha=pval;
    elseif any(strcmpi(pname,{'start_pos'}))
        start_pos=pval;
    elseif any(strcmpi(pname,{'flag_chooseoriginalMS'})) %
        flag_chooseoriginalMS=pval;
    end
end

if isempty(start_pos), start_pos=[floor(ratio(1)/2)+1,floor(ratio(2)/2)+1]; end
if isempty(band_HR) && flag_PAN==1, band_HR=1; end
if isempty(band_HR) && flag_PAN==0, band_HR=0; alpha=0; end
L1=sizes(1); L2=sizes(2); Nb=sizes(3);
Nb_PAN=band_HR;
if isempty(alpha), alpha='same'; end
if ischar(alpha)
    if any(strcmpi(alpha,{'even','same'})), alpha=Nb_PAN/(Nb+Nb_PAN);
    elseif strcmpi(alpha,'half'), alpha=0.5;
    else, warning('Unknown value for alpha'); alpha=Nb_PAN/(Nb+Nb_PAN);
    end
end


%% Inizializing masks

if ~isempty(type) && (~isempty(masktype_MS) || ~isempty(masktype_PAN))
    disp('Warning: The general mask overrides the separate mask assignment');
end

if strcmpi(type, 'CASSI'), masktype_MS='CASSI'; masktype_PAN='CASSI';
elseif strcmpi(type, 'none'), masktype_MS='none'; masktype_PAN='none';
elseif strncmpi(type, 'Kodak',5), masktype_MS=type; masktype_PAN=type;
elseif strcmpi(type, 'weighted'), masktype_MS='weighted'; masktype_PAN='weighted';
elseif any(strcmpi(type, {'random','Amba'})), masktype_MS='random'; masktype_PAN='random';
elseif strcmpi(type, 'binary'), masktype_MS='binary'; masktype_PAN='binary';
elseif strcmpi(type, 'KverCon'), masktype_MS='mindis'; masktype_PAN='vertical';
elseif strcmpi(type, 'KverUni'), masktype_MS='period'; masktype_PAN='vertical';
elseif strcmpi(type, 'KdiagCon'), masktype_MS='mindis'; masktype_PAN='diagonal';
elseif strcmpi(type, 'KdiagUni'), masktype_MS='period'; masktype_PAN='diagonal';
elseif strcmpi(type, 'ICIP2018')
    masktype_PAN='coverage'; Nb_PAN=1;
    if Nb==3 || Nb==4, masktype_MS='mindis'; else, masktype_MS='period'; end
% elseif strcmpi(type, 'Bayer'), masktype_MS='Bayer'; masktype_PAN='none';
% elseif strcmpi(type, 'Bayermod'), masktype_MS='Bayermod'; masktype_PAN='none';
% elseif strcmpi(type, 'Yamanaka'), masktype_MS='Yamanaka'; masktype_PAN='none';
elseif ischar(type), masktype_MS=type;
else, mask_MS=type;  % Capability of assigning pattern to the LR image mask
end

if isempty(masktype_PAN)
    if band_HR==0, masktype_PAN='none';
    else, masktype_PAN='coverage';
    end
end
if isempty(masktype_MS), masktype_MS='period'; end

%% Assigning extra fields to the mask

mask.shift=zeros(Nb+Nb_PAN,2);
mask.band_HR=Nb_PAN;

%% Assigning mask to the LR
if strcmpi(masktype_MS,'period')
    for jj=floor(sqrt(Nb)):-1:1
        if mod(Nb,jj)==0, break; end
    end
    mask_MS=reshape(1:Nb,[jj,Nb/jj]);
elseif strcmpi(masktype_MS,'mindis')
    if Nb>=4
        mask_MS=[1,2,3,4;3,4,1,2];
    elseif Nb==3
        mask_MS = double(imread('CFA_Condat_18.tif'))+1;
    elseif Nb==2
        mask_MS=[1,2;2,1];
    end
    if Nb>4, disp('Warning: Unavailable mask for given amount of bands'); end
elseif strcmpi(masktype_MS,'vertical')
        mask_MS=(1:Nb)';
elseif strcmpi(masktype_MS,'horizontal')
        mask_MS=1:Nb;
elseif any(strcmpi(masktype_MS,{'diagonal','Aggarwal'}))
    mask_MS=[];
    for ii=1:Nb
        mask_MS=cat(1,mask_MS,circshift(1:Nb,[0,ii-1]));
    end
elseif strcmpi(masktype_MS,'Bayer')
    mask_MS=[2,1;3,2];
    if Nb~=3, disp('Warning: Unavailable mask for given amount of bands'); end
elseif strcmpi(masktype_MS,'QuadBayer')
    mask_MS=[2,2,1,1;2,2,1,1;3,3,2,2;3,3,2,2];
    if Nb~=3, disp('Warning: Unavailable mask for given amount of bands'); end
elseif strcmpi(masktype_MS,'Bayermod')
    mask_MS=[1,2,3,2;2,1,2,3];
    if Nb~=3, disp('Warning: Unavailable mask for given amount of bands'); end
elseif strcmpi(masktype_MS,'Yamanaka')
    mask_MS=[2,1;2,3];
    if Nb~=3, disp('Warning: Unavailable mask for given amount of bands'); end
elseif strcmpi(masktype_MS,'Lukac')
    mask_MS=[2,1;2,3;1,2;3,2];
    if Nb~=3, warning('Unavailable mask for given amount of bands'); end
elseif strncmpi(masktype_MS,'Kodak',5)
    if Nb==2, extra_var=1; else, extra_var=3; end
    if Nb==3, extra_var2=2; else, extra_var2=4; end
    if strcmpi(masktype_MS,'Kodakver3')
        mask_MS=[2,1;extra_var,extra_var2];
    else
        mask_MS=[2,1;2,1;extra_var,extra_var2;extra_var,extra_var2];
    end
    if Nb>4, disp('Warning: Unavailable mask for given amount of bands'); end
elseif strncmpi(masktype_MS,'BinaryTree',10)
    if strcmpi(masktype_MS,'BinaryTreeD')
        mask_MS=binarytree_mask(Nb,'d');
    elseif strcmpi(masktype_MS,'BinaryTreeS')
        mask_MS=binarytree_mask(Nb,'s');
    else
        mask_MS=binarytree_mask(Nb,'u');
    end
elseif strcmpi(masktype_MS,'random')
    mask_MS=0;
elseif any(strcmpi(masktype_MS,{'binary','CASSI'}))
    mask_MS=-1;
elseif strcmpi(masktype_MS,'weighted')
    mask_MS=-2;
elseif strcmpi(masktype_MS,'none')
    mask_MS=[];
end

init_mask=[];

if ~isempty(Nb_PAN)
    if strcmpi(masktype_PAN,'none')
        mask_PAN=0;
        flag_chooseoriginalMS=0;
        if Nb_PAN>0 && ~isempty(Nb_PAN), fprintf('Warning: No mask was assigned to the HR image\n'); end
    elseif any(strcmpi(masktype_PAN,{'default','dft','coverage'}))
        mask_PAN=ones(ratio);
        mask_PAN(start_pos(1),start_pos(2))=0;
        if Nb_PAN>1, disp('Warning: Unavailable mask for given amount of bands'); end  
    elseif strcmpi(masktype_PAN,'vertical')
        mask_PAN=[(1:Nb_PAN).',zeros(Nb_PAN,1)];
    elseif strcmpi(masktype_PAN,'horizontal')
        mask_PAN=[1:Nb_PAN;zeros(1,Nb_PAN)];
    elseif strcmpi(masktype_PAN,'diagonal')
        mask_PAN=[];
        row=zeros(1,2*Nb_PAN); row(1:2:end)=1:Nb_PAN;
        for ii=1:2*Nb_PAN
            mask_PAN=cat(1,mask_PAN,circshift(row,[0,ii-1]));
        end
    elseif strncmpi(masktype_PAN,'Kodak',5)
        if ~strcmpi(masktype_PAN,'Kodakver1')
            mask_PAN=[0,1];
        else
            mask_PAN=[1,0;0,1];
        end
        if Nb_PAN>1, disp('Warning: Unavailable mask for given amount of bands'); end
    elseif strcmpi(masktype_PAN,'random')
        mask_PAN=zeros([L1,L2]);
        num_zeros=round(alpha*numel(mask_PAN));
        idx_zeros=1:numel(mask_PAN);
        sample_length=[round(num_zeros/Nb_PAN)*ones(1,Nb_PAN-1),num_zeros-round(num_zeros/Nb_PAN)*(Nb_PAN-1)];
        for jj=1:Nb_PAN
            idx_pan=randsample(length(idx_zeros),sample_length(jj));
            mask_PAN(idx_zeros(idx_pan))=jj;
            idx_zeros(idx_pan)=[];
        end
        flag_chooseoriginalMS=0;
    elseif any(strcmpi(masktype_PAN,{'binary','CASSI'}))
        init_mask=zeros([L1*L1,Nb_PAN]);
        num_minusone=L1*L2;
        alpha_val2=0.5;
        sample_length=round(alpha_val2*L1*L2);
        for jj=1:Nb_PAN
            if strcmpi(masktype_PAN,'binary')
                init_mask(randsample(num_minusone,sample_length),jj)=1;
            else
                init_mask(randsample(num_minusone,sample_length),jj)=1/(Nb_PAN+Nb);
            end
        end
        init_mask=reshape(init_mask,[L1,L2,Nb_PAN]);
        mask_PAN=0;
        Nb_PAN=0;
        flag_chooseoriginalMS=0;
    elseif any(strcmpi(masktype_PAN,{'weighted'}))
        if ~strcmpi(masktype_MS,{'weighted'}), error('This mask just works in combo with weighted masks'); end
        alpha_val=alpha/(1-alpha)*Nb/Nb_PAN; alpha_val(isnan(alpha_val))=1;
        if abs(alpha_val-1)<1E-6
            init_mask=simplex_uniform(Nb+Nb_PAN,L1*L2);
            init_mask=reshape(init_mask,L1,L2,[]);
        else
            init_mask=reshape(drchrnd([alpha_val*ones(1,Nb_PAN),ones(1,Nb)],L1*L2),[L1,L2,Nb+Nb_PAN]);
        end
        mask_PAN=0;
        Nb=0; Nb_PAN=0;
        flag_chooseoriginalMS=0;
    end
end

mask.data=create_combo_mask(mask_PAN,mask_MS,[L1,L2,Nb],Nb_PAN);
mask.data=cat(3,init_mask,mask.data);

if strcmpi(masktype_PAN,'CASSI'), mask.shift(1:mask.band_HR,2)=(0:mask.band_HR-1).'; end
if strcmpi(masktype_MS,'CASSI'), mask.shift(mask.band_HR+1:end,2)=(max(mask.shift(:))+(1:Nb)).'; end

%% Chooses the closest true MS value (doesn't consider interpolated MS)
if flag_chooseoriginalMS==1  
    originalgrid=zeros([L1,L2]);
    originalgrid(start_pos(1):ratio(1):end,start_pos(2):ratio(2):end)=1;
    for kk=Nb_PAN+(1:Nb)
        maxcorr=0;
        for ii=-start_pos(1)+(1:ratio(1))
            for jj=-start_pos(2)+(1:ratio(2))
                mask_temp=circshift(mask.data(:,:,kk),[ii,jj]);
                corry=sum(originalgrid.*abs(mask_temp),1:2);
                % fprintf('%d %d %f\n ',ii,jj,corry);
                if corry>maxcorr
                    maxcorr=corry;
                    mask.data(:,:,kk)=mask_temp;
                    mask.shift(kk,:)=[ii,jj]; 
                end
            end
        end
    end
end

mask.label='mask';
mask.type=type;
mask.type_LR=masktype_MS;
mask.type_HR=masktype_PAN;
mask.padding=find_zeropad(mask.data,mask.shift);
mask.ratio=ratio;
mask.start_pos=start_pos;
mask.mosaic=   @(I_in,mask_in) mask_join(I_in,mask_in);
mask.demosaic= @(I_mosaic,mask_in) mask_separate(I_mosaic,mask_in);

end