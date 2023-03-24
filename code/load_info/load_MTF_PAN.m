function  MTF_PAN = load_MTF_PAN( sensor_PAN )
%LOAD_MTF Loads the Nyquist gain of the MS sensor MTF
switch sensor_PAN
    case 'IKONOS'
        MTF_PAN=0.17;
    case 'GeoEye1'
        MTF_PAN=0.16;
    case 'QB'
        MTF_PAN=0.15;
    case 'WV2'
        MTF_PAN=0.11;
    case {'WV3','WV34bands'}
        MTF_PAN=0.14;
    case 'ALI'
        MTF_PAN=0.13;
    otherwise
        MTF_PAN=0.15;
end