function [ wavelength,bandwidth ] = load_wavelength( sensor,type,Bands )
%LOAD_WAVELENGTH Summary of this function goes here
%   Detailed explanation goes here
if nargin<=1, type='PAN'; end
if strcmpi(sensor,'Hyperion'); sensor='HYP'; end
if strcmpi(sensor,'IKO'); sensor='IKONOS'; end
if strcmpi(sensor,'CHRIS'); sensor='CHR'; end
if strncmpi(sensor,'WV2',3), sensor='WV2'; end
if strncmpi(sensor,'AVIRIS',6), sensor='AVIRIS'; end

current_folder=fileparts(mfilename('fullpath'));
data_folder=fullfile(current_folder,'..','..','data','relative_spectral_responses');

switch sensor
    case 'HYP'
        if strcmpi(type,'PAN')
            Cen_wl=585;
            bw=210;
        else
            Cen_wl=[355.5900,365.7600,375.9400,386.1100,396.2900,406.4600,416.6400,426.8200,436.9900,447.1700,457.3400,467.5200,477.6900,487.8700,498.0400,508.2200,518.3900,528.5700,538.7400,548.9200,559.0900,569.2700,579.4500,589.6200,599.8000,609.9700,620.1500,630.3200,640.5000,650.6700,660.8500,671.0200,681.2000,691.3700,701.5500,711.7200,721.9000,732.0700,742.2500,752.4300,762.6000,772.7800,782.9500,793.1300,803.3000,813.4800,823.6500,833.8300,844.0000,851.9200,854.1800,862.0100,864.3500,872.1000,874.5300,882.1900,884.7000,892.2800,894.8800,902.3600,905.0500,912.4500,915.2300,922.5400,925.4100,932.6400,935.5800,942.7300,945.7600,952.8200,955.9300,962.9100,966.1100,972.9900,976.2800,983.0800,986.4600,993.1700,996.6300,1003.3000,1006.8100,1013.3000,1016.9800,1023.4000,1027.1600,1033.4900,1037.3300,1043.5900,1047.5100,1053.6900,1057.6800,1063.7900,1073.8900,1083.9900,1094.0900,1104.1900,1114.1900,1124.2800,1134.3800,1144.4800,1154.5800,1164.6800,1174.7700,1184.8700,1194.9700,1205.0700,1215.1700,1225.1700,1235.2700,1245.3600,1255.4600,1265.5600,1275.6600,1285.7600,1295.8600,1305.9600,1316.0500,1326.0500,1336.1500,1346.2500,1356.3500,1366.4500,1376.5500,1386.6500,1396.7400,1406.8400,1416.9400,1426.9400,1437.0400,1447.1400,1457.2300,1467.3300,1477.4300,1487.5300,1497.6300,1507.7300,1517.8300,1527.9200,1537.9200,1548.0200,1558.1200,1568.2200,1578.3200,1588.4200,1598.5100,1608.6100,1618.7100,1628.8100,1638.8100,1648.9000,1659.0000,1669.1000,1679.2000,1689.3000,1699.4000,1709.5000,1719.6000,1729.7000,1739.7000,1749.7900,1759.8900,1769.9900,1780.0900,1790.1900,1800.2900,1810.3800,1820.4800,1830.5800,1840.5800,1850.6800,1860.7800,1870.8700,1880.9800,1891.0700,1901.1700,1911.2700,1921.3700,1931.4700,1941.5700,1951.5700,1961.6600,1971.7600,1981.8600,1991.9600,2002.0600,2012.1500,2022.2500,2032.3500,2042.4500,2052.4500,2062.5500,2072.6500,2082.7500,2092.8400,2102.9400,2113.0400,2123.1400,2133.2400,2143.3400,2153.3400,2163.4300,2173.5300,2183.6300,2193.7300,2203.8300,2213.9300,2224.0300,2234.1200,2244.2200,2254.2200,2264.3200,2274.4200,2284.5200,2294.6100,2304.7100,2314.8100,2324.9100,2335.0100,2345.1100,2355.2100,2365.2000,2375.3000,2385.4000,2395.5000,2405.6000,2415.7000,2425.8000,2435.8900,2445.9900,2456.0900,2466.0900,2476.1900,2486.2900,2496.3900,2506.4800,2516.5900,2526.6800,2536.7800,2546.8800,2556.9800,2566.9800,2577.0800];
            bw=[11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3871,11.3784,11.3538,11.3133,11.2580,11.1907,11.1119,11.0245,10.9321,10.8368,10.7407,10.6482,10.5607,10.4823,10.4147,10.3595,10.3188,10.2942,10.2856,10.2980,10.3349,10.3909,10.4592,10.5322,10.6004,10.6562,10.6933,10.7058,10.7276,10.7907,10.8833,10.9938,11.1044,11.1980,11.2600,11.2824,11.2822,11.0457,11.2816,11.0457,11.2809,11.0457,11.2797,11.0457,11.2782,11.0457,11.2771,11.0457,11.2765,11.0457,11.2756,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0457,11.2754,11.0451,11.2754,11.0423,11.2754,11.0372,11.2754,11.0302,11.2754,11.0218,11.0122,11.0013,10.9871,10.9732,10.9572,10.9418,10.9248,10.9065,10.8884,10.8696,10.8513,10.8335,10.8154,10.7979,10.7822,10.7663,10.7520,10.7385,10.7270,10.7174,10.7091,10.7022,10.6970,10.6946,10.6937,10.6949,10.6996,10.7058,10.7163,10.7283,10.7437,10.7612,10.7807,10.8034,10.8267,10.8534,10.8818,10.9110,10.9422,10.9743,11.0074,11.0414,11.0759,11.1108,11.1461,11.1811,11.2156,11.2496,11.2826,11.3146,11.3460,11.3753,11.4037,11.4302,11.4538,11.4760,11.4958,11.5133,11.5286,11.5404,11.5505,11.5580,11.5621,11.5634,11.5617,11.5563,11.5477,11.5346,11.5193,11.5002,11.4789,11.4548,11.4279,11.3994,11.3688,11.3366,11.3036,11.2696,11.2363,11.2007,11.1666,11.1333,11.1018,11.0714,11.0424,11.0155,10.9912,10.9698,10.9508,10.9355,10.9230,10.9139,10.9083,10.9069,10.9057,10.9013,10.8951,10.8854,10.8740,10.8591,10.8429,10.8242,10.8039,10.7820,10.7592,10.7342,10.7092,10.6834,10.6572,10.6312,10.6052,10.5803,10.5560,10.5328,10.5101,10.4904,10.4722,10.4552,10.4408,10.4285,10.4197,10.4129,10.4088,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077,10.4077];
        end
    case 'ALI'
        if strcmpi(type,'PAN')
            Cen_wl=585;
            bw=210;
        else
            Cen_wl=[443,482.5,565,660,790,867.5,1250,1650,2215];
            bw=[20,65,80,60,30,45,100,200,270];
        end
    case 'QB'
        if strcmpi(type,'PAN')
            Cen_wl=729;
            bw=648;
        else
            Cen_wl=[487.5,545,650,816.5];
            bw=[115,154,120,203];
        end
    case 'IKONOS'
        if strcmpi(type,'PAN')
            Cen_wl=675;
            bw=450;
        else
            Cen_wl=[480.5,550.5,665,805];
            bw=[71,89,66,96];
        end
    case 'GeoEye1'
        if strcmpi(type,'PAN')
            Cen_wl=675;
            bw=350;
        else
            Cen_wl=[480,545,677.5,850];
            bw=[60,70,45,140];
        end
    case 'WV2'
        if strcmpi(type,'PAN')
            Cen_wl=625;
            bw=350;
        else
            Cen_wl=[427,478,546,608,659,724,831,908]; % Coastal blue, Blue, Green, Yellow, Red, Red edge, NIR1 and NIR2
            bw=[50,60,70,40,60,40,125,180];
        end
    case 'WV24bands'
        if strcmpi(type,'PAN')
            Cen_wl=625;
            bw=350;
        else
            Cen_wl=[478,546,659,831]; % Blue, Green, Red, NIR1
            bw=[60,70,60,125];
        end
    case 'WV3'
        if strcmpi(type,'PAN')
            Cen_wl=625;
            bw=450;
        else
            Cen_wl=[425,480,545,605,660,725,832.5,950,1210,1570,1660,1730,2165,2205,2260,2330];
            bw=[50,60,70,40,60,40,125,180,30,40,40,40,40,40,60,70];
        end
    case 'WV34bands'
        if strcmpi(type,'PAN')
            Cen_wl=625;
            bw=450;
        else
            Cen_wl=[480,545,660,832.5];
            bw=[60,70,60,125];
        end
    case 'Pleiades'
        if strcmpi(type,'PAN')
            Cen_wl=600;
            bw=360;
        else
            Cen_wl=[490,550,660,845];
            bw=[80,80,80,140];
        end
    case 'DE2'
        if strcmpi(type,'PAN')
            Cen_wl=675;
            bw=450;
        else
            Cen_wl=[465,545,660,825];
            bw=[90,70,120,130];
        end
    case {'CAVE','U260'}
        Cen_wl=400:10:700;
        bw=10*ones(1,length(Cen_wl));
    case {'CZ','Nuance'}
        Cen_wl=420:10:720;
        bw=10*ones(1,length(Cen_wl));
    case {'FNA','Pulnix'}
        Cen_wl=400:10:720;
        bw=10*ones(1,length(Cen_wl));
    case {'TT','VariSpec'}
        Cen_wl=420:10:720;
        bw=10*ones(1,length(Cen_wl));
    case 'AVIRIS'
        load(fullfile(data_folder,'AVIRIS_Spectral_Responses.mat'),'central','bandwidth');
        Cen_wl=central;
        bw=bandwidth;
    case 'AVIRIS220'
        load(fullfile(data_folder,'AVIRIS_Spectral_Responses.mat'),'central','bandwidth');
        Cen_wl=central([2:32,34:96,98:160,162:end]);
        bw=bandwidth([2:32,34:96,98:160,162:end]);
    case 'AVIRIS176'
        load(fullfile(data_folder,'AVIRIS_Spectral_Responses.mat'),'central','bandwidth');
        Cen_wl=central([2:32,34:96,98:160,162:end]);
        bw=bandwidth([2:32,34:96,98:160,162:end]);
        Cen_wl=Cen_wl(1:176);
        bw=bw(1:176);
        % warning('Central Bandwidths are likely not accurate');
    case 'ROSIS'
        Cen_wl=430:(850-430)/102:850;
        bw=15.2*ones(1,length(Cen_wl));
        % warning('Central Bandwidths are likely not accurate');
    case {'CHRIS','CHR'}
        if strcmpi(type,'PAN')
            Cen_wl=550;
            bw=300;
        else
            Cen_wl=415:(1050-415)/17:1050;
            bw=(1050-415)/17*ones(1,length(Cen_wl));
        end
        % warning('Central Bandwidths are likely not accurate');
    case 'none'
        if strcmpi(type,'PAN')
            Cen_wl=560; bw=360;
        else
            error('Unknown wavelengths characteristics');
        end
    otherwise
        error('Unknown wavelengths characteristics');
end

if nargin<=2 || isempty(Bands), Bands=length(Cen_wl); end

if strcmpi(type,'PAN')
    wavelength=Cen_wl;
    bandwidth=bw;
else
    wavelength=Cen_wl(Bands);
    bandwidth=bw(Bands);
end

end

