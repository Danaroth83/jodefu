% BINARY MASK PERIOD CHECKER
% Note: It returns an empty output if the mask is not periodic

function period=period_custom(mask)

    max_period=16;

    period=[0,0];

    for ii=1:min(max_period,size(mask,1))
        if isequal(mask(1:end-ii,:,:),mask(ii+1:end,:,:)), period(1)=ii; break; end
    end

    for ii=1:min(max_period,size(mask,2))
        if isequal(mask(1:end-ii,:,:),mask(ii+1:end,:,:)), period(2)=ii; break; end
    end

    if any(period==0), period=[size(mask,1),size(mask,2)]; end
end