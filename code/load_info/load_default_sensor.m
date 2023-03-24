function [sensor_lq,sensor_hq]=load_default_sensor(place)

if strncmpi(place,'RdJ',3) || strncmpi(place,'Sydney',6)
    sensor_lq='HYP';
    sensor_hq='WV3';
elseif strncmpi(place,'Beijing',7)
    sensor_lq='WV34bands';
    sensor_hq='WV34bands';
elseif strncmpi(place,'SanFrancisco',12)
    sensor_lq='QB';
    sensor_hq='QB';
elseif strncmpi(place,'SaoPaulo',8)
    sensor_lq='HYP';
    sensor_hq='IKONOS';
elseif strncmpi(place,'Sofia',5) || strncmpi(place,'Sudbury',7) || strncmpi(place,'Fuji',4)
    sensor_lq='HYP';
    sensor_hq='ALI';    
elseif strncmpi(place,'WV2',3)
    sensor_lq='WV2_test';
    sensor_hq='WV2_test';
elseif strncmpi(place,'Washington',10)
    sensor_lq='WV24bands';
    sensor_hq='WV24bands';
elseif strncmpi(place,'Rome',4) ||strncmpi(place,'Rio',3) || strncmpi(place,'Stockholm',9)
    sensor_lq='WV2';
    sensor_hq='WV2';
elseif strncmpi(place,'Tripoli',7) || strncmpi(place,'Janeiro',7)
    sensor_lq='WV3';
    sensor_hq='WV3';
elseif strncmpi(place,'China',5) || strncmpi(place, 'Fields', 6)
    sensor_lq='IKONOS';
    sensor_hq='IKONOS';
elseif strncmpi(place,'Hobart',6)
    sensor_lq='GeoEye1';
    sensor_hq='GeoEye1';
elseif strncmpi(place,'Indianapolis',12)
    sensor_lq='QB';
    sensor_hq='QB';
elseif strncmpi(place,'Tls',3) || strncmpi(place,'Toulouse',8)
    sensor_lq='IKONOS';
    sensor_hq='IKONOS';
elseif strncmpi(place,'Chris',5)
    sensor_lq='CHR';
    sensor_hq='CHR';
elseif strncmpi(place,'HYP',3) || strncmpi(place,'Paris',5)
    sensor_lq='HYP';
    sensor_hq='ALI';
elseif strncmpi(place,'Moffett',7)
    sensor_lq='AVIRIS176';
    sensor_hq='none';
elseif strncmpi(place,'Vancouver',9)
    sensor_lq='DE2';
    sensor_hq='DE2';
elseif strncmpi(place,'Pleiades',9)
    sensor_lq='Pleiades';
    sensor_hq='Pleiades';
elseif strncmpi(place,'Pavia',5) || strncmpi(place,'ROSIS',5)
    sensor_lq='ROSIS';
    sensor_hq='none';
elseif strncmpi(place,'CZ',2)
    sensor_lq='Nuance';
    sensor_hq='none';
elseif strncmpi(place,'CAVE',4)
    sensor_lq='U260';
    sensor_hq='none';
elseif strncmpi(place,'FNA',3)
    sensor_lq='Pulnix';
    sensor_hq='none';
elseif strncmpi(place,'TT',2)
    sensor_lq='VariSpec';
    sensor_hq='none';    
else
    sensor_lq='none';
    sensor_hq='none';
end