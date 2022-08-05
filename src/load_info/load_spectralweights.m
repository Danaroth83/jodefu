function spectralweights=load_spectralweights(sensor_MS,sensor_PAN,Bands,type)

if nargin<=4, type='MS'; end

[~,~,wavelength_nm,Response_MS]=load_spectralresponse(sensor_MS,type,Bands );
[~,~,wavelength_nm_PAN,Response_PAN]=load_spectralresponse( sensor_PAN,'PAN');
Response_MS=abs(Response_MS./repmat(max(Response_MS,[],2),[1,size(Response_MS,2)]));
Response_PAN=abs(Response_PAN./max(Response_PAN(:)));

wavelength_nm_min_idx=find(wavelength_nm>=wavelength_nm_PAN(1),1,'first');
wavelength_nm_max_idx=find(wavelength_nm<=wavelength_nm_PAN(end),1,'last');
Response_MS=Response_MS(:,wavelength_nm_min_idx:wavelength_nm_max_idx);

Response_PAN_new=zeros(size(Response_PAN,1),length(wavelength_nm_min_idx:wavelength_nm_max_idx));
for ii=1:size(Response_PAN,1)
    Response_PAN_new(ii,:) = interp1(wavelength_nm_PAN,Response_PAN(ii,:),wavelength_nm(wavelength_nm_min_idx:wavelength_nm_max_idx),'pchip');
end

% [~,idx_maxwv]=max(Response_MS,[],2);
% div_Response_PAN=shiftdim(Response_PAN_new(:,idx_maxwv),1);
% div_Response_PAN(div_Response_PAN<0.1)=1;
% Response_PAN_norm=Response_PAN_new./div_Response_PAN;

Response_PAN_norm=Response_PAN_new;

spectralweights=sum(Response_MS.*Response_PAN_norm,2);
spectralweights=spectralweights/sum(spectralweights(:));
