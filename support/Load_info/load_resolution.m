function [ spatial_res,radiometric_res ] = load_resolution( sensor, im_tag, type )
%LOAD_RESOLUTION Loads the spatial resolution in meters of specific satellite sensors
if nargin<=2 || isempty(type), type='MS'; end
if nargin<=1, im_tag=[]; end
if strcmpi(sensor,'IKO'), sensor='IKONOS'; end
if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'CHRIS'), sensor='CHR'; end
if strcmpi(sensor,'GE1'), sensor='GeoEye1'; end

if isequal(type,1) || strcmp(type,'PAN')
    switch sensor
        case 'QB'
            spatial_res=0.6;
        case 'IKONOS'
            spatial_res=0.8;
        case 'WV1'
            spatial_res=0.5;
        case {'WV2','WV2_test','WV24bands'}
            spatial_res=0.4;
        case {'WV3','WV34bands'}
            if strncmpi(im_tag,'Adelaide',8) || strncmpi(im_tag,'RdJ',3) || strncmpi(im_tag,'Tripoli',7) || strncmpi(im_tag,'Janeiro',7)
                spatial_res=0.3;
            elseif strncmpi(im_tag,'Beijing',7)
                spatial_res=0.4;
            elseif strncmpi(im_tag,'Sydney',6)
                spatial_res=0.5;
            else
                spatial_res=0.4;
            end
        case 'GeoEye1'
            spatial_res=0.5;
        case 'DE2'
            spatial_res=1;
        case {'ALI','HYP'}
            spatial_res=10;
        case {'AVIRIS','AVIRIS176','AVIRIS220','AVIRIS_Mof'} % Simulated Panchromatic
            spatial_res=6;
        case 'Pleiades'
            spatial_res=0.7;
        case {'CHR','CHRIS'}
            spatial_res=4.5;
        otherwise
            spatial_res=NaN;
    end
else
    switch sensor
        case 'QB'
            spatial_res=2.4;
        case 'IKONOS'
            spatial_res=3.2;
        case 'WV1'
            spatial_res=2;
        case {'WV2','WV24bands'}
            spatial_res=1.6;
        case {'WV3','WV34bands'}
            if strncmpi(im_tag,'Adelaide',8) || strncmpi(im_tag,'RdJ',3) || strncmpi(im_tag,'Tripoli',7)|| strncmpi(im_tag,'Janeiro',7)
                spatial_res=1.2;
            elseif strncmpi(im_tag,'Beijing',7)
                spatial_res=1.6;
            elseif strncmpi(im_tag,'Sydney',6)
                spatial_res=2;
            else
                spatial_res=1.6;
            end
        case 'GeoEye1'
            spatial_res=2;
        case 'DE2'
            spatial_res=4;
        case 'ALI'
            spatial_res=30;
        case 'HYP'
            spatial_res=30;
        case {'AVIRIS','AVIRIS176','AVIRIS220','AVIRIS_Mof'}
            spatial_res=30;
        case 'Pleiades'
            spatial_res=2.8;
        case 'ROSIS'
            spatial_res=1.3;
         case {'CHR','CHRIS'}
            spatial_res=18;
        otherwise
            spatial_res=NaN;
    end
end

if nargout>=2
    switch sensor
        case {'WV1','WV2','WV24bands','WV3','WV34bands','IKONOS','QB','GeoEye1','DE2','Pleiades'}
            radiometric_res=11;
        case {'ALI','HYP','ROSIS'}
            radiometric_res=15;
         case {'CHR','CHRIS'}
            radiometric_res=16;
        case {'AVIRIS','AVIRIS_Mof','AVIRIS176','AVIRIS220'}
            radiometric_res=13;
        case {'CAVE','U260'}
            radiometric_res=16;
            if strncmpi(im_tag,'CAVE32',6) || strncmpi(im_tag,'watercolors',6)
                radiometric_res=8;
            end
        case {'CS','Nuance'}
            radiometric_res=0;
        case {'FNA','Pulnix'}
            radiometric_res=0;
        case {'TT','Varispec'}
            radiometric_res=0;
        otherwise
            radiometric_res=NaN;
    end
end
end

