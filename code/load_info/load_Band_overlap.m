function Band_overlap=load_Band_overlap(sensor_lq,sensor_hq,type_lq,type_hq)
if nargin<=2, type_lq='MS'; end
if nargin<=3, type_hq='PAN'; end

% cd Assignment
% [cent_lq,disp_lq]=load_spectralresponse(sensor_lq,type_lq,im_tag);
% [cent_hq,disp_hq]=load_spectralresponse(sensor_hq,type_hq,im_tag);
% cd ..
% dispersionratio=1.5;
% min_lq=cent_lq-dispersionratio*disp_lq;
% max_lq=cent_lq+dispersionratio*disp_lq;
% min_hq=cent_hq-dispersionratio*disp_hq;
% max_hq=cent_hq+dispersionratio*disp_hq;

[min_lq,max_lq]=load_minmaxspectralresponse(sensor_lq,type_lq);
[min_hq,max_hq]=load_minmaxspectralresponse(sensor_hq,type_hq);


Nhq=numel(min_hq);
Band_overlap=cell(1,Nhq);
for ii=1:Nhq
    min_hq_temp=min_hq(ii);
    max_hq_temp=max_hq(ii);
    Band_overlap{ii}=find((min_lq>=min_hq_temp & min_lq<=max_hq_temp) | (max_lq>=min_hq_temp & max_lq<=max_hq_temp));
end
if strcmpi(type_hq,'PAN'), Band_overlap=Band_overlap{1}; end

end