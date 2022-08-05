function [place,cut,bandselection,sensor_lq,sensor_hq_PAN,sensor_hq_MS]=parse_im_tag(im_tag)

sensor_database={'WV3','WV34bands','WV2','WV1','IKONOS','QB',...
    'QuickBird','GE1','GeoEye1','AVIRIS','HYP','Hyperion','ALI','Chris',...
    'CHR','U260','Nuance','Pulnix','Varispec'};

idx=strfind(im_tag,'_');
if isempty(idx)
    place=im_tag;
    cut='dftcut';
    bandselection='dft';
    [sensor_lq,sensor_hq_PAN]=load_default_sensor(place);
    sensor_hq_MS=sensor_hq_PAN;
else
    place=im_tag(1:idx(1)-1);
    if length(im_tag)-idx(1)>=3 && strcmpi(im_tag(idx(1)+(1:3)),'cut')
        if numel(idx)>=2
            cut=im_tag(idx(1)+1:idx(2)-1);
        else
            cut=im_tag(idx(1)+1:end);
        end
        idx=idx(2:end);
    else
        cut='dftcut';
    end
    if numel(idx)==0
        bandselection='dft';
        [sensor_lq,sensor_hq_PAN]=load_default_sensor(place);
        sensor_hq_MS=sensor_hq_PAN;
    else
        if ~any(strcmpi(im_tag(idx(end)+1:end),sensor_database))
            bandselection=im_tag(idx(end)+1:end);
            im_tag=im_tag(1:idx(end)-1);
            idx=idx(1:end-1);
        else
            bandselection='dft';
        end
        if numel(idx)==0
            [sensor_lq,sensor_hq_PAN]=load_default_sensor(place);
            sensor_hq_MS=sensor_hq_PAN;
        elseif numel(idx)==1
            sensor_lq=load_default_sensor(place);
            sensor_hq_PAN=im_tag(idx(end)+1:end);
            sensor_hq_MS=sensor_hq_PAN;
        elseif numel(idx)==2
            sensor_lq=im_tag(idx(1)+1:idx(2)-1);
            sensor_hq_PAN=im_tag(idx(end)+1:end);
            sensor_hq_MS=sensor_hq_PAN;
        elseif numel(idx)==3
            sensor_lq=im_tag(idx(1)+1:idx(2)-1);
            sensor_hq_MS=im_tag(idx(2)+1:idx(3)-1);
            sensor_hq_PAN=im_tag(idx(end)+1:end);
        end
    end
end

if strncmpi(place,'Sofia',5), place='Sofia'; end
if strncmpi(place,'Sudbury',7), place='Sudbury'; end
if strncmpi(place,'Tls',3), place='Toulouse'; end
if strncmpi(place,'HYP',3), place='Paris'; end
if strncmpi(place,'ROSIS',5), place='Pavia'; end
if strncmpi(place,'WV2',5), place='Rome'; end
if strncmpi(place,'CHR',3), place='Chris'; end

if strcmpi(place,'Beijing')
    if strcmpi(sensor_lq,'WV3'), sensor_lq='WV34bands'; end
    if strcmpi(sensor_hq_MS,'WV3'), sensor_hq_MS='WV34bands'; end
    if strcmpi(sensor_hq_PAN,'WV3'), sensor_hq_PAN='WV34bands'; end
end
% if strcmpi(place,'Hobart'), place='Hobart1'; end
% if strcmpi(place,'Indianapolis'), place='Indianapolis1'; end

end