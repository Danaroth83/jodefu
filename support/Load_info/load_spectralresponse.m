%%  LOAD SENSOR SPECTRAL RESPONSES
%
% Author:
% Daniele Picone
%
% Description:
% Loads the centroid and the dispersion of an assigned spectral response
%
% Usage:
% [central,dispersion,wavelength_nm,response] = load_spectralresponse(sensor,type,Bands)
%
% Input:
% sensor: String describing the sensor
% type: Type of sensor data (eg: 'PAN','MS','HS')
% Bands: Subset of bands taken from the sensor (length: Nb)
%
% Output:
% central: central wavelengths of the spectral responses (length: Nb)
% dispersion: Bandwidth of each spectral response (length: Nb)
% wavelength_nm: Wavelenght in nanometer of the spectral response (length: Nw)
% response: Spectral response of the sensor (sizes: Nb x Nw)


function [ central,dispersion,wavelength_nm,response ] = load_spectralresponse( sensor,type,Bands )

if nargin<=1, type='MS'; end
if nargin<=2, Bands=[]; end

if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'HYP') && strcmpi(type,'PAN'), sensor='ALI'; end
if strcmpi(sensor,'AVIRIS_Mof'), sensor='AVIRIS176'; end
if any(strcmpi(sensor,{'DE2','Pleiades','Cave','U260','CZ','Nuance','FNA','Pulnix','ROSIS','CHR','Varispec','none'}))
    [central,dispersion]=load_wavelength(sensor,type,Bands);
    wavelength_nm=(floor(min(central-dispersion/2)):1:ceil(max(central+dispersion/2))).';
    response=zeros(length(central),length(wavelength_nm));
    for ii=1:length(central)
        response(ii,wavelength_nm>=central(ii)-dispersion(ii)/2 & wavelength_nm<=central(ii)+dispersion(ii)/2)=1;
    end
    return;
end


current_folder=fileparts(mfilename('fullpath'));
spectrum_file=fullfile(current_folder,'..','..','data','relative_spectral_responses',[sensor,'_Spectral_Responses.mat']);
load(spectrum_file,'Spectral_Responses_Matrix','wavelength_nm');

if ~strcmp(type,'PAN')
    if nargin<=2
        if ~any(strcmpi(sensor,{'HYP','AVIRIS','AVIRIS176'}))
            Bands=1:size(Spectral_Responses_Matrix,1)-1;
        else
            Bands=1:size(Spectral_Responses_Matrix,1);
        end
    end
    Nbands=length(Bands);
    Nwavelengths=length(wavelength_nm);
    response=Spectral_Responses_Matrix(Bands,:);
    summation=sum(response,2);
    wavelength_mat=repmat(wavelength_nm.',[Nbands,1]);
    central=(sum(wavelength_mat.*response,2)./summation);
    dispersion=sqrt(sum(response.*(wavelength_mat-repmat(central,[1,Nwavelengths])).^2,2)./summation);
    dispersion=dispersion';
    central=central';
else
    response=Spectral_Responses_Matrix(end,:);
    summation=sum(response,2);
    central=sum(wavelength_nm.'.*response)/summation;
    dispersion=sqrt(sum(response.*(wavelength_nm.'-repmat(central,[1,length(wavelength_nm)])).^2)/summation);
end

end

