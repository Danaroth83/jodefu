%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           ATWT decomposition in log2(ratio) levels, based on the Starck 
%           and Murtagh (S&M) Filter [Strang96]
%          
% Interface:
%           [Image_LP,Image_HP] = ATWT(Orig_Image,ratio)
%
% Inputs:
%           Orig_Image: Original image
%           ratio:      Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           Image_LP:   Low Pass part of the original image 
%           Image_HP:   High Pass part of the original image 
% 
% References:
%           [Strang96]      G. Strang and T. Nguyen, Wavelets and Filter Banks, 2nd ed. Wellesley,
%                           MA: Wellesley Cambridge Press, 1996.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Image_LP,Image_HP] = ATWT_dec(Orig_Image,ratio)


[Height,Width,Bands]=size(Orig_Image);

h=[1 4 6 4 1 ]/16;
g=[0 0 1 0 0 ]-h;
htilde=[ 1 4 6 4 1]/16;
gtilde=[ 0 0 1 0 0 ]+htilde;
h=sqrt(2)*h;
g=sqrt(2)*g;
htilde=sqrt(2)*htilde;
gtilde=sqrt(2)*gtilde;
WF={h,g,htilde,gtilde};
Levels = ceil(log2(ratio));
for i=1:Bands    
    WT = ndwt2(Orig_Image(:,:,i),Levels,WF);
    
    for ii = 2 : numel(WT.dec), WT.dec{ii} = zeros(size(WT.dec{ii})); end
    
    Image_LP(:,:,i) = indwt2(WT,'c');
end

Image_HP = Orig_Image - Image_LP;

