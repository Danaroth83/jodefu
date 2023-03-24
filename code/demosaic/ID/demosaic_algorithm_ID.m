%% DEMOSAICK INTENSITY DIFFERENCE
% 
% Author:
% Daniele Picone
%
% Description:
% This function realizes the demosaic of a certain mosaic with the
% procedure described in [1] and [2]. It also introduces a more extensive
% method to manage linear interpolation for masks which don't contain a 
% single non-zero element within its periodicity.
%
% Usage:
% I_out = demosaic_algorithm_ID(I_MSFA,mask,type,period,central_wv)
%
% Input:
% I_MSFA: Mosaicked image (sizes: N1 x N2)
% mask: Binary mask (sizes: N1 x N2 x Nb)
% period: period of the mask
% central_wv: Central wavelengths of the bands in nm (sizes: Nb x 1)
%
% Output:
% I_out: Demosaicked image
%
% References:
% [1] Mihoubi S., Losson O., Mathon B. and Macaire L., "Multispectral
% Demosaicing using Intensity-based Spectral Correlation" IPTA 2015, 
% pp. 461-466
% [2] Mihoubi S., Losson O., Mathon B. and Macaire L.,  "Multispectral 
% demosaicing using pseudo-panchromatic image", IEEE TCI, IEEE, 2017, 3 (4),
% pp.982-995.

function I_hat = demosaic_algorithm_ID(I_MSFA,cfa,type,period,central_wv)

    current_folder=fileparts(mfilename('fullpath'));
    addpath(fullfile(current_folder,'..','..','interpolation'));

    if isfield(cfa,'data'),cfa=cfa.data; end
    if nargin<=4 && isfield(I_MSFA,'wavelength'), central_wv=I_MSFA.wavelength; end
    if isfield(I_MSFA,'data'), I_MSFA=I_MSFA.data; end
    if nargin<=2, type='ID'; end
    
    [L1,L2,Nb]=size(cfa);
    if strcmpi(type,'ISD') && isempty(central_wv)
        error('ISD needs sensor spectral response to operate');
    end
    
    % mask_check=sum(cfa~=0,3);
    % if max(mask_check(:))>=1, fprintf('Warning: Demosaic is supposed to work on non overlapping masks\n'); end
    
    % current_folder=fileparts(mfilename('fullpath'));
    % addpath(fullfile(current_folder,'..','Quality_indices'));
    
    % Find periodicity in the mask
    % [H,period]=period_custom(cfa);
    % H=[1:period(2),period(2)-1:-1:1]*[1:period(1),period(1)-1:-1:1].'/period(1)/period(2)

    if strcmpi(type,'SD')
        
        I_tilde=I_MSFA.*cfa;
        % I_hat=imfilter(I_tilde,H,'replicate');
        I_hat=interp2D_mosaic(I_tilde,'mask',cfa,'method','WB','nn',100);
        cfa_rep=repmat(permute(cfa,[1,2,4,3]),[1,1,Nb,1]);
        I_tilde_rep=repmat(permute(I_tilde,[1,2,4,3]),[1,1,Nb,1]);
        I_hat_rep=repmat(I_hat,[1,1,1,Nb]);
        Delta=I_hat_rep.*cfa_rep-I_tilde_rep;
        % Delta_tilde=imfilter(Delta,H,'replicate');
            Delta_tilde=zeros(size(Delta));
            for kk=1:Nb
                Delta_tilde=interp2D_mosaic(Delta(:,:,:,kk),'mask',cfa_rep(:,:,:,kk),'method','WB','nn',100);
            end
        I_hat=permute(sum((I_hat_rep-Delta_tilde).*repmat(cfa,[1,1,1,Nb]),3),[1,2,4,3]);
        
    elseif strcmpi(type,'ISD')
        
        central_wv_rep=repmat(central_wv(:),[1,Nb]);
        central_wv_rep2=repmat(central_wv(:)',[Nb,1]);
        sigma=1.74;
        
        N=floor(exp(-(abs(central_wv_rep-central_wv_rep2)-100)/20/sigma));        
        I_tilde=I_MSFA.*cfa;
        % I_hat=imfilter(I_tilde,H,'replicate');
        I_hat=interp2D_mosaic(I_tilde,'mask',cfa,'method','WB','nn',100);

        cfa_rep=repmat(permute(cfa,[1,2,4,3]),[1,1,Nb,1]);
        I_tilde_rep=repmat(permute(I_tilde,[1,2,4,3]),[1,1,Nb,1]);
        I_hat_rep=repmat(I_hat,[1,1,1,Nb]);
        Delta=I_hat_rep.*cfa_rep-I_tilde_rep;
        % Delta_tilde=imfilter(Delta,H,'replicate');
            Delta_tilde=zeros(size(Delta));
            for kk=1:Nb
                Delta_tilde=interp2D_mosaic(Delta(:,:,:,kk),'mask',cfa_rep(:,:,:,kk),'method','WB','nn',100);
            end
        I_hat=permute(sum((I_hat_rep-Delta_tilde).*repmat(cfa,[1,1,1,Nb]),3),[1,2,4,3]);
        I_hat_rep=repmat(I_hat,[1,1,1,Nb]);
        
        for ii=1:max(N(:))
            I_hat_rep=I_hat_rep.*cfa_rep-I_tilde_rep;
            idx=repmat(shiftdim((N<=ii),-2),[L1,L2,1,1]);
            Delta(idx)=I_hat_rep(idx);
        end
        % Delta_tilde=imfilter(Delta,H,'replicate');
            Delta_tilde=zeros(size(Delta));
            for kk=1:Nb
                Delta_tilde=interp2D_mosaic(Delta(:,:,:,kk),'mask',cfa_rep(:,:,:,kk),'method','WB','nn',100);
            end
        I_hat=permute(sum((I_hat_rep-Delta_tilde).*repmat(cfa,[1,1,1,Nb]),3),[1,2,4,3]);
        
    elseif strcmpi(type,'ID')
        
        I_M=filterM(I_MSFA,cfa,period);
        I_tilde=I_MSFA.*cfa;
        % I_M_rep=I_M,[1,1,Nb]);
        Delta=I_tilde-I_M.*cfa;
        % Delta_hat=imfilter(Delta,H,'replicate');
        Delta_hat=interp2D_mosaic(Delta,'mask',cfa,'method','WB','nn',100);
        I_hat=I_M+Delta_hat;
        
    elseif strcmpi(type,'IID')
        
        Nbiter=7;
        I_hat=demosaic_algorithm_ID(I_MSFA,cfa,'ID',period);
        I_tilde=I_MSFA.*cfa;
        for ii=1:Nbiter
            I_M=sum(I_hat,3)/Nb;
            I_M_rep=repmat(I_M,[1,1,Nb]);
            Delta=I_tilde-I_M_rep.*cfa;
            % Delta_hat=imfilter(Delta,H,'replicate');
            Delta_hat=interp2D_mosaic(Delta,'mask',cfa,'method','WB','nn',100);
            I_hat=I_M_rep+Delta_hat;
        end
            
    end
end

