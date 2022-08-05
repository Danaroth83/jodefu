function [vcut,hcut,edge_cut_lq,edge_cut_hq,Qblocks_size]=load_cut(im_tag,cut_label,sensor_lq,sensor_hq,type_hq)
if nargin<=3, sensor_hq=sensor_lq; end
if nargin<=4, type_hq='PAN'; end
if strcmpi(sensor_lq,'IKO'), sensor_lq='IKONOS'; end
if strcmpi(sensor_hq,'IKO'), sensor_hq='IKONOS'; end
if strcmpi(sensor_lq,'Hyperion'), sensor_lq='HYP'; end
if strcmpi(sensor_hq,'Hyperion'), sensor_hq='HYP'; end
if strcmpi(sensor_lq,'CHRIS'), sensor_lq='CHR'; end
if strcmpi(sensor_hq,'CHRIS'), sensor_hq='CHR'; end
if strcmpi(sensor_lq,'GE1'), sensor_lq='GeoEye1'; end
if strcmpi(sensor_hq,'GE1'), sensor_hq='GeoEye1'; end
if strcmpi(sensor_lq,'Deimos-2'), sensor_lq='DE2'; end
if strcmpi(sensor_hq,'Deimos-2'), sensor_hq='DE2'; end

ratio=load_resolution(sensor_lq,im_tag,'MS')/load_resolution(sensor_hq,im_tag,type_hq);
edge_cut_lq=0;
edge_cut_hq=0;
Qblocks_size=32;

if strncmpi(im_tag,'Beijing',7)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=138;
        switch cut_label
            case 'cut2'
                vcut=20+(1:144); hcut=68+(1:96);      %Full image
            case 'cut3'
                vcut=8+(1:156); hcut=128+(1:36);
            case 'cut4'
                vcut=-74+(1:288); hcut=148+(1:126);
            case 'cut5'
                vcut=4+(1:156); hcut=4+(1:156);
            otherwise
                vcut=68+(1:96); hcut=116+(1:48);
        end
    elseif any(strcmpi(sensor_lq,{'WV3','WV34bands'}))
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        switch cut_label
            case {'cut2','cut128'}
                vcut=1972+(-63:64);  hcut=1714+(-63:64);
            case {'cut3','cut256'}
                vcut=1972+(-127:128);  hcut=1714+(-127:128);
            otherwise
                vcut=1972+(-255:256);  hcut=1714+(-255:256); %Bird's nest cut
                % vcut=1716+(1:512);  hcut=1458+(1:512); %Bird's nest cut
        end
        vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=138*ratio;
    end
elseif strncmpi(im_tag,'RdJ',3)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=69;
        switch cut_label
            case 'cut2'
                vcut=1+(1:108); hcut=2+(1:132);  % Full Image
            otherwise
                vcut=1+(1:108); hcut=12+(1:120);
        end
    elseif strcmpi(sensor_lq,'WV3')
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        vcut=44+(1:512); hcut=1424-256+(1:512);
        vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=69*ratio;
    end
elseif strncmpi(im_tag,'SaoPaulo',8)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=138;
        switch cut_label
            case 'cut2'
                vcut=40+(1:192); hcut=40+(1:192);
            otherwise
                vcut=40+(1:192); hcut=88+(1:144);
        end
    elseif strcmpi(sensor_lq,'IKONOS')
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        switch cut_label
            case 'cut2'
                vcut=800+(1:256); hcut=800+(1:256);
            otherwise
                vcut=800+(1:512); hcut=800+(1:512);
        end
        vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=138*ratio;
    end
elseif strncmpi(im_tag,'SanFrancisco',12)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=138;
        switch cut_label
            case 'cut2'
                vcut=29:160; hcut=29:160;
                % vcut=25:160; hcut=25:160;
            case 'cut3'
                vcut=2+(1:294); hcut=22+(1:132);
            otherwise
                vcut=88+(1:72); hcut=40+(1:120);
        end
    elseif strcmpi(sensor_lq,'QB')
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        switch cut_label
            case 'cut2'
                vcut=1000+(1:256); hcut=1504+(1:256);
            otherwise
                vcut=872+(1:512); hcut=1376+(1:512);
        end
        vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=138*ratio;
    end
elseif strncmpi(im_tag,'Sydney',6)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=138;
        vcut=1:132; hcut=1:132;
    elseif strcmpi(sensor_lq,'WV3')
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        vcut=500+(1:512); hcut=1100+(1:512);
        vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=138*ratio;
    end
elseif strncmpi(im_tag,'Stockholm',9)
    edge_cut_lq=92;
    edge_cut_hq=92*ratio;
    if any(strcmpi(cut_label,{'cut3','cut256'}))
        vcut=1008+(-127:128); hcut=880+(-127:128);
    else
        % vcut=1008+(-255:256); hcut=880+(-255:256);
        vcut=1136+(-255:256); hcut=752+(-255:256);
    end
    vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
elseif strncmpi(im_tag,'Washington',10)
    edge_cut_lq=92;
    edge_cut_hq=92*ratio;
    if any(strcmpi(cut_label,{'cut3','cut256'}))
        vcut=1503+(1:256); hcut=1328+(1:256);
    else
        vcut=1375+(1:512); hcut=1200+(1:512);
    end
    vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
elseif strncmpi(im_tag,'Tripoli',7)
    edge_cut_lq=92;
    edge_cut_hq=92*ratio;
    if any(strcmpi(cut_label,{'cut2','cut128'}))
        vcut=1256+(1:128); hcut=881+(1:128);
    elseif any(strcmpi(cut_label,{'cut3','cut256'}))
        % vcut=1128+(1:256); hcut=753+(1:256);
        vcut=1256+(-255:0); hcut=880+(-127:128);
    else
        vcut=1256+(-255:256); hcut=880+(-255:256);
    end
    vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
elseif strncmpi(im_tag,'Janeiro',7)
    edge_cut_lq=92;
    edge_cut_hq=92*ratio;
    if any(strcmpi(cut_label,{'cut2','cut128'}))
        vcut=1256+(-63:64); hcut=2256+(-63:64);
    elseif any(strcmpi(cut_label,{'cut3','cut256'}))
        % vcut=1256+(-127:128); hcut=2256+(-127:128);
        vcut=1256+48+(-255:0); hcut=2256-64-32+(-255:0);
    else
        vcut=1256+(-255:256); hcut=2256+(-255:256);
    end
    vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
elseif strncmpi(im_tag,'Vancouver',9)
    edge_cut_lq=92;
    edge_cut_hq=92*ratio;
    switch cut_label
        case 'cut2'
            vcut=511+(1:500); hcut=299+(1:500);
        otherwise
            vcut=650+(1:512); hcut=263+(1:512);
    end
    vcut=vcut-edge_cut_lq; hcut=hcut-edge_cut_lq;
elseif strncmpi(im_tag,'Sofia',5) || strncmpi(im_tag,'Sudbury',7)
    Qblocks_size=30;
    edge_cut_lq=43;
    edge_cut_hq=43*ratio;
    switch cut_label
        case 'cut2'
            % vcut=97:288; hcut=1:192;   % Square; RR power of 2; black corner
            vcut=1:384; hcut=49:144;   % Rectangular; RR power of 2; no black corner
            % vcut=97:288; hcut=19:174;  % Rectangular; RR not power of 2; no black corner
            % vcut=113:272; hcut=17:176; % Square; RR not power of 2; no black corner
        otherwise
            if strncmpi(im_tag,'Sofia1GST',9) || strncmpi(im_tag,'Sudbury1GST',11)
                vcut=97:288; hcut=22:177;
            else
                vcut=97:288; hcut=19:174;
            end
    end
elseif strncmpi(im_tag,'Fuji',4)
    Qblocks_size=30;
    edge_cut_lq=43;
    edge_cut_hq=43*ratio;
    vcut=1:105; hcut=1:105;
elseif strncmpi(im_tag,'China',5)
    switch cut_label
        case 'cut2'
            vcut=48+(1:128); hcut=48+(1:128);
        otherwise
            vcut=1:300; hcut=1:300;
    end
elseif strncmpi(im_tag,'Hobart',6)
    switch cut_label
        case 'cut2'
            vcut=1:128; hcut=1:128;
        otherwise
            vcut=1:512; hcut=1:512;
    end
elseif strncmpi(im_tag,'Rio',3)
    vcut=1:64; hcut=1:64;
elseif strncmpi(im_tag,'Rome',4)
    vcut=1:300; hcut=1:300;
elseif strncmpi(im_tag,'Indianapolis',12)
    switch cut_label
        case 'cut2'
            vcut=1:100; hcut=1:100;
        otherwise
            vcut=1:128; hcut=1:128;
    end
elseif strncmpi(im_tag,'Tls',3) || strncmpi(im_tag,'Toulouse',8)
    vcut=1:128; hcut=1:128;
elseif strncmpi(im_tag,'CHR',3)
    vcut=1:100; hcut=1:100;
elseif strncmpi(im_tag,'HYP',3) || strncmpi(im_tag,'Paris',5)
    Qblocks_size=30;
    switch cut_label
        case 'cut2'
            vcut=31:130; hcut=181:280;
        case 'cut3'
             vcut=11:160; hcut=141:290;
        case 'cut4'
            vcut=1:361; hcut=1:348;
        otherwise
            vcut=24+(1:108); hcut=158+(1:108);
    end
elseif strncmpi(im_tag,'Pavia',5) || strncmpi(im_tag,'ROSIS',5)
    vcut=50+(1:512); hcut=20+(1:256);
elseif strncmpi(im_tag,'Moffett',7)
    Qblocks_size=30;
    switch cut_label
        case 'cut1'
            % vcut=7+(1:64); hcut=2+(1:32);
            vcut=35+(1:320); hcut=10+(1:160);
        case 'cut2'
            % vcut=1:78; hcut=1:36; % Almost full image
            vcut=1:390; hcut=1:180;
        otherwise
            % vcut=1:79; hcut=1:37;
            vcut=1:395; hcut=1:185;
    end
elseif strncmpi(im_tag,'Pleiades',8)
    switch cut_label
        case 'cut2'
            vcut=256+(1:512); hcut=256+(1:512);
        case 'cut3'
            vcut=1:256; hcut=1:256;
        otherwise
            vcut=1:1024; hcut=1:1024;
    end
elseif strncmpi(im_tag,'CZ',2)
    switch cut_label
        case {'cut1','cut512'}
            vcut=520+(-255:256); hcut=696+(-255:256);
        case {'cut2','cut256'}
            vcut=520+(-127:128); hcut=696+(-127:128);
        case {'cut3','cut128'}
            vcut=520+(-63:64); hcut=696+(-63:64);
        otherwise
            vcut=1:1040; hcut=1:1392;
    end
elseif strncmpi(im_tag,'CAVE',4)
    switch cut_label
        case {'cut2','cut256'}
            vcut=256+(-127:128); hcut=256+(-127:128);
        case {'cut3','cut128'}
            vcut=256+(-63:64); hcut=256+(-63:64);
        otherwise
            vcut=1:512; hcut=1:512;
    end
elseif strncmpi(im_tag,'TT',2)
    switch cut_label
        case {'cut2','cut256'}
            vcut=250+(-127:128); hcut=250+(-127:128);
        case {'cut3','cut128'}
            vcut=250+(-63:64); hcut=250+(-63:64);
        otherwise
            vcut=[]; hcut=[];
    end    
elseif strncmpi(im_tag,'FNA',3)
    if strncmpi(im_tag,'FNAc',4)
        switch cut_label
            case {'cut1','cut512'}
                vcut=510+(-255:256); hcut=670+(-255:256);
            case {'cut2','cut256'}
                vcut=510+(-127:128); hcut=670+(-127:128);
            case {'cut3','cut128'}
                vcut=510+(-63:64); hcut=670+(-63:64);
            otherwise
                vcut=[]; hcut=[];
        end
    elseif strncmpi(im_tag,'FNAh',4)
        switch cut_label
            case {'cut1','cut512'}
                vcut=508+(-255:256); hcut=634+(-255:256);
            case {'cut2','cut256'}
                vcut=508+(-127:128); hcut=634+(-127:128);
            case {'cut3','cut128'}
                vcut=508+(-63:64); hcut=634+(-63:64);
            otherwise
                vcut=[]; hcut=[];
        end
     elseif strncmpi(im_tag,'FNA',3)
        switch cut_label
            case {'cut1','cut512'}
                vcut=341+(-255:256); hcut=209+(-255:256);
            case {'cut2','cut256'}
                vcut=341+(-127:128); hcut=209+(-127:128);
            case {'cut3','cut128'}
                vcut=341+(-63:64); hcut=209+(-63:64);
            otherwise
                vcut=[]; hcut=[];
        end
    end
end
    
end