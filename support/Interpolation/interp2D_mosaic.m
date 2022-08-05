%% SCATTERED 2D INTERPOLATION FUNCTION

% Description:
% This script allows to perform a scattered interpolation for an image over%
% a grid with missing samples:
%
% Usage:
% I_out=interp2D_mosaic(I_in,'Field',Value)
% Eg:
% I_in=double('cameraman.tif');
% mask=ones(size(I_in));
% mask(2:2:end,2:2:end)=0;
% I_out=interp2D_mosaic(I_in,'mask',mask,'method','RBF_gaussian','param',0.7,'nn',100)
%
% Input:
% I_in: Input image on a N1xN2 grid with missing samples
%
% Output: Interpolated image
%
% Fields:
% 'mask': A matrix of the same sizes of I_in, where zeros mark missing
%         samples. If it has more bands, I_in is considered to be a
%         mosaicked version various images with missing samples
% 'nn': Number of neighbors around the target sample used to perform the interpolation
% 'method': The interpolation method, to choose amongst this list
%       - 'RBF': Gaussian Radial Basis Function
%       - 'RBF_spline': Thin Plate Radial Basis Function
%       - 'RBF_iq': Inverse Quadratic Radial Basis Function
%       - 'RBF_mq': Multi Quadratic Radial Basis Function
%       - 'RBF_imq': Inverse Multi Quadratic Radial Basis Function;
%       - 'IDW': Inverse weighting interpolation (also known as Shepard interpolation)
%       - 'WB': Bilinear over Delauney triangulated points
%       - 'NN': Nearest Neighbour
% 'param': For RBF, it is the scaling distance parameter (default: 1).
%          For IDW it is the power assigned to the distance (default: 2)



function I_out=interp2D_mosaic(I_in,varargin)

    current_folder=fileparts(mfilename('fullpath'));
    addpath(fullfile(current_folder,'methods'));
    
    % Parser of input variables   
    mask=[];
    method='bilinear';
    param=[];
    Nsamples=[];
    
    for ii=1:2:numel(varargin)
        pname=varargin{ii};
        pval=varargin{ii+1};
        if strcmpi(pname,'mask')                                                % Position of samples to interplate (marked with zeros)
            mask=pval;
        elseif any(strcmpi(pname,{'type','method'}))                            % Interpolation method
            method=pval;
        elseif any(strcmpi(pname,{'parameter','param','sigma','k'}))            % Method parameter variable
            param=pval;
        elseif any(strcmpi(pname,{'Nsamples','neighbors','neighbours','nn'}))   % Amount of neighbours to consider for interpolation
            Nsamples=pval;
        end
    end
    
    
    if isempty(mask), mask=ones(size(I_in)); mask(I_in==0)=0; end
    if isempty(Nsamples), Nsamples=size(I_in,1)*size(I_in,2); end
   
    
    % Fix for 4D+ images
    sizes=size(mask);
    I_in=reshape(I_in,size(mask,1),size(mask,2),[]);
    mask=reshape(mask,size(mask,1),size(mask,2),[]);
    
    % Find positions of zeros in the mask
    toler=10E-9*max(mask(:));
    mask_idxnonzero=mask>toler;
    % if any(sum(mask_idxnonzero,3)>=2), fprintf('Masks are supposed to be non overlapping\n'); end
    
    % Fix for mosaicked images
    I_in=I_in.*mask_idxnonzero;
    
    % Fix for methods not using param
    if any(strcmpi(method,{'WB','bilinear'})) || any(strcmpi(method,{'NN','nearest','nearestneighbour'})), param=1; end
    Nparam=max(length(param),1);
    
    % Initialization and memory allocation
    I_in(mask_idxnonzero)=I_in(mask_idxnonzero)./mask(mask_idxnonzero);
    I_out=reshape(I_in,[size(mask,1)*size(mask,2),size(mask,3)]);
    I_out=repmat(I_out,[1,1,Nparam]);
    
    % Creation of the sample grid
    [y,x]=meshgrid(1:size(mask,2),1:size(mask,1));
    x=x(:); y=y(:);
    
    for kk=1:size(mask,3)
        I_current=reshape(I_in(:,:,kk),[],1);
        mask_current=reshape(mask_idxnonzero(:,:,kk),[],1);
        xi=x(mask_current);
        yi=y(mask_current);
        fi=I_current(mask_current);
        xo=x(~mask_current);
        yo=y(~mask_current);
        
        if any(strcmpi(method,{'WB','bilinear'}))
            % addpath(fullfile(method_folder,'bilinear'));
            period=period_custom(mask_idxnonzero(:,:,kk));
            if sum(mask_idxnonzero(1:period(1),1:period(2),kk),1:2)==1
                H=[1:period(1),period(1)-1:-1:1]'*[1:period(2),period(2)-1:-1:1]/period(1)/period(2);
                I_temp=imfilter(I_in(:,:,kk),H,'circular');
                I_temp=I_temp(:);
                fo=I_temp(~mask_current);
            elseif isequal(mask_idxnonzero(1:period(1),1:period(2),kk),[1,0;0,1]) || isequal(mask_idxnonzero(1:period(1),1:period(2),kk),[0,1;1,0]) 
                I_temp=imfilter(I_in(:,:,kk),[0,1,0;1,4,1;0,1,0]/4,'circular');
                I_temp=I_temp(:);
                fo=I_temp(~mask_current);
            else % isempty(period) || isequal(period,[size(mask,1),size(mask,2)])
                F=scatteredInterpolant(xi,yi,fi,'linear','linear');
                fo=F(xo(:),yo(:));
            end
        elseif any(strcmpi(method,{'NN','nearest','nearestneighbour'}))
            F=scatteredInterpolant(xi,yi,fi,'nearest','nearest');
            fo=F(xo(:),yo(:));
        elseif strncmpi(method,'RBF',3)
            idx_string=strfind(method,'_');
            if isempty(idx_string), RBF_shape=[]; else, RBF_shape=method(idx_string+1:end); end
            fo=interp2D_RBF(xi,yi,fi,xo,yo,Nsamples,param,RBF_shape);
        elseif strncmpi(method,'Shepard',7) || strncmpi(method,'IDW',3)
            idx_string=strfind(method,'_');
            if isempty(idx_string), IDW_method=[]; else, IDW_method=method(idx_string+1:end); end
            fo=interp2D_Shepard(xi,yi,fi,xo,yo,Nsamples,param,IDW_method);
        elseif strncmpi(method,'Kriging',7)
            error('Kriging is currently not supported');              
        end
        I_out(~mask_current,kk,:)=fo;
    end
    
    % Reshaping the output in image form
    I_out=reshape(I_out,[size(mask,1),size(mask,2),size(mask,3),Nparam]);
    I_out=reshape(I_out,[sizes,Nparam]);
end