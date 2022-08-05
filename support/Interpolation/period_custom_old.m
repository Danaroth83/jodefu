function period=period_custom(mask)
    
    current_folder = fileparts(mfilename('fullpath'));
    addpath(fullfile(current_folder,'..','Validation'));
    max_test=17;

    max_test=min(max_test,[size(mask,1),size(mask,2)]);
    
    for ii=1:max_test(1)
        if SCC(mask,circshift(mask,ii,1),'none')>0.9, break; end
    end
    for jj=1:max_test(2)
        if SCC(mask,circshift(mask,jj,2),'none')>0.9, break; end
    end
    if ii==max_test(1), ii=size(mask,1); end
    if jj==max_test(2), jj=size(mask,2); end
    period=[ii,jj];
end