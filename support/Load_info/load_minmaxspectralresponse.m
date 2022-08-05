function [ min_sr,max_sr ] = load_minmaxspectralresponse( sensor,type,Bands )
%LOAD_MINMAXSPECTRALRESPONSE
%   Loads the edges of a spectral response


if nargin<=1, type='MS'; end
if nargin<=2, Bands=[]; end

if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'IKO'), sensor='IKONOS'; end
if strcmpi(sensor,'HYP') && strcmpi(type,'PAN'), sensor='ALI'; end
if strncmpi(sensor,'WV2',3), sensor='WV2'; end
if strncmpi(sensor,'AVIRIS',6), sensor='AVIRIS'; end
% if nargin<=2, im_tag=[]; end
% if strncmpi(im_tag,'Beijing',7) && strcmpi(sensor,'WV3'), sensor='WV34bands'; end
if strcmpi(sensor,'none') && strcmpi(type,'PAN')
    min_sr=400; max_sr=740; disp('A simulated PAN sensor was employed (Wavelenghts: 400-740 nm)');
    return;
end
if any(strcmpi(sensor,{'DE2','U260','Nuance','ROSIS','Pleiades','CHR'}))
    [central,dispersion]=load_wavelength(sensor,'MS',Bands);
    min_sr=central-dispersion/2;
    max_sr=central+dispersion/2;
    return;
end


stopband_amp=0.5;

current_folder=fileparts(mfilename('fullpath'));
spectrum_file=fullfile(current_folder,'..','..','data','relative_spectral_responses',[sensor,'_Spectral_Responses.mat']);
load(spectrum_file,'Spectral_Responses_Matrix','wavelength_nm');

if ~strcmp(type,'PAN')
    if isempty(Bands)
        if ~any(strcmpi(sensor,{'HYP','AVIRIS'}))
            Bands=1:size(Spectral_Responses_Matrix,1)-1;
        else
            Bands=1:size(Spectral_Responses_Matrix,1);
        end
    end
    response=Spectral_Responses_Matrix(Bands,:);
else
    Bands=1;
    response=Spectral_Responses_Matrix(end,:);
end

Nbands=length(Bands);
Nwavelengths=length(wavelength_nm);
maxamplitude=max(response,[],2);
higher=response>stopband_amp*repmat(maxamplitude,[1,Nwavelengths]);
min_sr=zeros(1,Nbands);
max_sr=zeros(1,Nbands);
for ii=1:Nbands
    min_sr(ii)=wavelength_nm(find(higher(ii,:),1));
    max_sr(ii)=wavelength_nm(find(higher(ii,:),1,'last'));
end
    
end