%% FIND ZERO PADDING

% Description:
% Given the maximum shift of the mask, finds the most efficient zero
% padding to the original image to perform all the mosaicking operations
%
% Usage:
% padding=padding=find_zeropad(mask,shift)
%
% Input:
% mask:  The mosaicking mask, whose sizes are L1 x L2 x Nb
%        (image columns x rows x band)
% shift: The requested shift of the mask (sizes: 2 x Nb), where the first
%        column defines the shift in the horizontal direction and the
%        second in the vertical direction
%
% Output:
% padding: A 4 he necessary padding for the image, whose order is
%          [left, up, right, down]


function padding=find_zeropad(mask,shift)

    padding=[-min(min(shift(:,1)),0),-min(min(shift(:,2)),0),...
             max(max(shift(:,1)),0),max(max(shift(:,2)),0)];

    for ii=padding(1):-1:1
        if all(mask(1:ii,:,:)==0), padding(1)=padding(1)-ii; break; end
    end
    for ii=padding(2):-1:1
        if all(mask(:,1:ii,:)==0), padding(2)=padding(2)-ii; break; end
    end
    for ii=padding(3):-1:1
        if all(mask(end-ii+1:end,:,:)==0), padding(3)=padding(3)-ii; break; end
    end
    for ii=padding(4):-1:1
        if all(mask(:,end-ii+1:end,:)==0), padding(4)=padding(4)-ii; break; end
    end

end