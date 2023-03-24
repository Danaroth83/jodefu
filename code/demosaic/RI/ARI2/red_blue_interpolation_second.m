function [red, blue] = red_blue_interpolation_second(green, red, blue, mask, eps, DynamicRange)

% inverse mask
imr = zeros(size(mask(:,:,1)));
img = zeros(size(mask(:,:,1)));
imb = zeros(size(mask(:,:,1)));
imr(find(mask(:,:,1)==0)) = 1;
img(find(mask(:,:,2)==0)) = 1;
imb(find(mask(:,:,3)==0)) = 1;
imask = mask;
imask(:,:,1) = imr;
imask(:,:,2) = img;
imask(:,:,3) = imb;

%% Iterpolate R and B at G pixels
% Step (i): iterative directional interpolation
% initial linear interpolation 
F1 = [1/2,0,1/2];
Guider1 = red+imfilter(red,F1,'replicate').*mask(:,:,2);
Guideg1 = green;
Guideb1 = blue+imfilter(blue,F1,'replicate').*mask(:,:,2);
F2 = F1';
Guider2 = red+imfilter(red,F2,'replicate').*mask(:,:,2);
Guideg2 = green;
Guideb2 = blue+imfilter(blue,F2,'replicate').*mask(:,:,2);

% initial guided filter window size for RI
h = 2;
v = 2;
% initial guided filter window size for MLRI
h2 = 2;
v2 = 0;
% maximum iteration number
itnum = 2;

% initialization of iteration criteria
  RI_w2R1 = ones(size(mask(:,:,1)))*1e32;
  RI_w2R2 = ones(size(mask(:,:,1)))*1e32;
MLRI_w2R1 = ones(size(mask(:,:,1)))*1e32;
MLRI_w2R2 = ones(size(mask(:,:,1)))*1e32;
  RI_w2B1 = ones(size(mask(:,:,1)))*1e32;
  RI_w2B2 = ones(size(mask(:,:,1)))*1e32;
MLRI_w2B1 = ones(size(mask(:,:,1)))*1e32;
MLRI_w2B2 = ones(size(mask(:,:,1)))*1e32;

% initial guide image for RI/MLRI
  RI_Guideg1 = Guideg1;
  RI_Guider1 = Guider1;
  RI_Guideb1 = Guideb1;
  RI_Guideg2 = Guideg2;
  RI_Guider2 = Guider2;
  RI_Guideb2 = Guideb2;
MLRI_Guideg1 = Guideg1;
MLRI_Guider1 = Guider1;
MLRI_Guideb1 = Guideb1;
MLRI_Guideg2 = Guideg2;
MLRI_Guider2 = Guider2;
MLRI_Guideb2 = Guideb2;

% initialization of interpolated R and B values
  RI_R1 = Guider1;
  RI_R2 = Guider2;
MLRI_R1 = Guider1;
MLRI_R2 = Guider2;
  RI_B1 = Guideb1;
  RI_B2 = Guideb2;
MLRI_B1 = Guideb1;
MLRI_B2 = Guideb2;

%% Iterative horizontal and vertical interpolation
for ittime = 1:itnum
	% generate horizontal and vertical tentative estimate by RI
    M = ones(size(mask(:,:,1))); % mask
	RI_tentativeR1 = guidedfilter(RI_Guideg1, RI_Guider1, M, h, v, eps);
	RI_tentativeB1 = guidedfilter(RI_Guideg1, RI_Guideb1, M, h, v, eps);
	RI_tentativeR2 = guidedfilter(RI_Guideg2, RI_Guider2, M, v, h, eps);
	RI_tentativeB2 = guidedfilter(RI_Guideg2, RI_Guideb2, M, v, h, eps);
	% generate horizontal tentative estimate by MLRI
	F1 = [-1,0,2,0,-1];
	difR = imfilter(MLRI_Guider1, F1,'replicate');
	difG = imfilter(MLRI_Guideg1, F1,'replicate');
	difB = imfilter(MLRI_Guideb1, F1,'replicate');
	MLRI_tentativeR1 = guidedfilter_MLRI(MLRI_Guideg1, MLRI_Guider1, M, difG, difR, imask(:,:,2), h2, v2, eps);
	MLRI_tentativeB1 = guidedfilter_MLRI(MLRI_Guideg1, MLRI_Guideb1, M, difG, difB, imask(:,:,2), h2, v2, eps);
	% generate vertical tentative estimate by MLRI
	F2 = [-1;0;2;0;-1];
	difR = imfilter(MLRI_Guider2, F2,'replicate');
	difG = imfilter(MLRI_Guideg2, F2,'replicate');
	difB = imfilter(MLRI_Guideb2, F2,'replicate');
	MLRI_tentativeR2 = guidedfilter_MLRI(MLRI_Guideg2, MLRI_Guider2, M, difG, difR, imask(:,:,2), v2, h2, eps);
	MLRI_tentativeB2 = guidedfilter_MLRI(MLRI_Guideg2, MLRI_Guideb2, M, difG, difB, imask(:,:,2), v2, h2, eps);
    
	% calculate residuals of RI and MLRI
	  RI_residualR1 = (red  -   RI_tentativeR1) .* imask(:,:,2);
	  RI_residualB1 = (blue -   RI_tentativeB1) .* imask(:,:,2);
	  RI_residualR2 = (red  -   RI_tentativeR2) .* imask(:,:,2);
	  RI_residualB2 = (blue -   RI_tentativeB2) .* imask(:,:,2);
	MLRI_residualR1 = (red  - MLRI_tentativeR1) .* imask(:,:,2);
	MLRI_residualB1 = (blue - MLRI_tentativeB1) .* imask(:,:,2);
	MLRI_residualR2 = (red  - MLRI_tentativeR2) .* imask(:,:,2);
	MLRI_residualB2 = (blue - MLRI_tentativeB2) .* imask(:,:,2);

	% horizontal and vertical linear interpolation of residuals
	K1 = [1/2,0,1/2];
	  RI_residualR1 = imfilter(  RI_residualR1, K1, 'replicate');
	  RI_residualB1 = imfilter(  RI_residualB1, K1, 'replicate');
	MLRI_residualR1 = imfilter(MLRI_residualR1, K1, 'replicate');
	MLRI_residualB1 = imfilter(MLRI_residualB1, K1, 'replicate');
	K2 = K1';
	  RI_residualR2 = imfilter(  RI_residualR2, K2, 'replicate');
	  RI_residualB2 = imfilter(  RI_residualB2, K2, 'replicate');
	MLRI_residualR2 = imfilter(MLRI_residualR2, K2, 'replicate');
	MLRI_residualB2 = imfilter(MLRI_residualB2, K2, 'replicate');
    
	% add tentative estimate
	  RI_R1 = (  RI_tentativeR1 +   RI_residualR1) .* mask(:,:,2);
	  RI_B1 = (  RI_tentativeB1 +   RI_residualB1) .* mask(:,:,2);
	  RI_R2 = (  RI_tentativeR2 +   RI_residualR2) .* mask(:,:,2);
	  RI_B2 = (  RI_tentativeB2 +   RI_residualB2) .* mask(:,:,2);
	MLRI_R1 = (MLRI_tentativeR1 + MLRI_residualR1) .* mask(:,:,2);
	MLRI_B1 = (MLRI_tentativeB1 + MLRI_residualB1) .* mask(:,:,2);
	MLRI_R2 = (MLRI_tentativeR2 + MLRI_residualR2) .* mask(:,:,2);
	MLRI_B2 = (MLRI_tentativeB2 + MLRI_residualB2) .* mask(:,:,2);

   %% Step(ii): adaptive selection of iteration at each pixel
    % calculate iteration criteria
	  RI_criR1 = (  RI_Guider1 -   RI_tentativeR1) .* M;
	  RI_criB1 = (  RI_Guideb1 -   RI_tentativeB1) .* M;
	  RI_criR2 = (  RI_Guider2 -   RI_tentativeR2) .* M;
	  RI_criB2 = (  RI_Guideb2 -   RI_tentativeB2) .* M;
	MLRI_criR1 = (MLRI_Guider1 - MLRI_tentativeR1) .* M;
	MLRI_criB1 = (MLRI_Guideb1 - MLRI_tentativeB1) .* M;
	MLRI_criR2 = (MLRI_Guider2 - MLRI_tentativeR2) .* M;
	MLRI_criB2 = (MLRI_Guideb2 - MLRI_tentativeB2) .* M;

    % calculate gradient of iteration criteria
	F1 = [-1,0,1];
	  RI_difcriR1 = abs(imfilter(  RI_criR1, F1, 'replicate'));
	  RI_difcriB1 = abs(imfilter(  RI_criB1, F1, 'replicate'));
	MLRI_difcriR1 = abs(imfilter(MLRI_criR1, F1, 'replicate'));
	MLRI_difcriB1 = abs(imfilter(MLRI_criB1, F1, 'replicate'));
	F2 = F1';
	  RI_difcriR2 = abs(imfilter(  RI_criR2, F2, 'replicate'));
	  RI_difcriB2 = abs(imfilter(  RI_criB2, F2, 'replicate'));
	MLRI_difcriR2 = abs(imfilter(MLRI_criR2, F2, 'replicate'));
	MLRI_difcriB2 = abs(imfilter(MLRI_criB2, F2, 'replicate'));
    
	% absolute value of iteration criteria
	  RI_criR1 = abs(  RI_criR1);
	  RI_criB1 = abs(  RI_criB1);
	  RI_criR2 = abs(  RI_criR2);
	  RI_criB2 = abs(  RI_criB2);
	MLRI_criR1 = abs(MLRI_criR1);
	MLRI_criB1 = abs(MLRI_criB1);
	MLRI_criR2 = abs(MLRI_criR2);
	MLRI_criB2 = abs(MLRI_criB2);
    
	% directional map of iteration criteria
	  RI_criR1 =   RI_criR1 +   RI_criB1;
	  RI_criB1 =   RI_criB1 +   RI_criR1;
	  RI_criR2 =   RI_criR2 +   RI_criB2;
	  RI_criB2 =   RI_criB2 +   RI_criR2;
	MLRI_criR1 = MLRI_criR1 + MLRI_criB1;
	MLRI_criB1 = MLRI_criB1 + MLRI_criR1;
	MLRI_criR2 = MLRI_criR2 + MLRI_criB2;
	MLRI_criB2 = MLRI_criB2 + MLRI_criR2;
    
	% directional gradient map of iteration criteri
	  RI_difcriR1 =   RI_difcriR1 +   RI_difcriB1;
	  RI_difcriB1 =   RI_difcriB1 +   RI_difcriR1;
	  RI_difcriR2 =   RI_difcriR2 +   RI_difcriB2;
	  RI_difcriB2 =   RI_difcriB2 +   RI_difcriR2;
	MLRI_difcriR1 = MLRI_difcriR1 + MLRI_difcriB1;
	MLRI_difcriB1 = MLRI_difcriB1 + MLRI_difcriR1;
	MLRI_difcriR2 = MLRI_difcriR2 + MLRI_difcriB2;
	MLRI_difcriB2 = MLRI_difcriB2 + MLRI_difcriR2;
    
    % smoothing of iteration criteria
	sigma = 2;
	F1 = fspecial('gaussian', [5,5], sigma);
	  RI_criR1 = imfilter(  RI_criR1, F1, 'replicate');
	MLRI_criR1 = imfilter(MLRI_criR1, F1, 'replicate');
	  RI_criB1 = imfilter(  RI_criB1, F1, 'replicate');
	MLRI_criB1 = imfilter(MLRI_criB1, F1, 'replicate');
	  RI_difcriR1 = imfilter(  RI_difcriR1, F1, 'replicate');
	MLRI_difcriR1 = imfilter(MLRI_difcriR1, F1, 'replicate');
	  RI_difcriB1 = imfilter(  RI_difcriB1, F1, 'replicate');
	MLRI_difcriB1 = imfilter(MLRI_difcriB1, F1, 'replicate');
	F2 = fspecial('gaussian', [5,5],sigma);
	  RI_criR2 = imfilter(  RI_criR2, F2, 'replicate');
	MLRI_criR2 = imfilter(MLRI_criR2, F2, 'replicate');
	  RI_criB2 = imfilter(  RI_criB2, F2, 'replicate');
	MLRI_criB2 = imfilter(MLRI_criB2, F2, 'replicate');
	  RI_difcriR2 = imfilter(  RI_difcriR2, F2, 'replicate');
	MLRI_difcriR2 = imfilter(MLRI_difcriR2, F2, 'replicate');
	  RI_difcriB2 = imfilter(  RI_difcriB2, F2, 'replicate');
	MLRI_difcriB2 = imfilter(MLRI_difcriB2, F2, 'replicate');

    % calcualte iteration criteria
	  RI_wR1 = (  RI_criR1.^2).*(  RI_difcriR1);
	  RI_wR2 = (  RI_criR2.^2).*(  RI_difcriR2);
	MLRI_wR1 = (MLRI_criR1.^2).*(MLRI_difcriR1);
	MLRI_wR2 = (MLRI_criR2.^2).*(MLRI_difcriR2);
	  RI_wB1 = (  RI_criB1.^2).*(  RI_difcriB1);
	  RI_wB2 = (  RI_criB2.^2).*(  RI_difcriB2);
	MLRI_wB1 = (MLRI_criB1.^2).*(MLRI_difcriB1);
	MLRI_wB2 = (MLRI_criB2.^2).*(MLRI_difcriB2);

    % find smaller criteria pixels
	  RI_piR1 = find(   RI_wR1 <   RI_w2R1);
	  RI_piR2 = find(   RI_wR2 <   RI_w2R2);
	MLRI_piR1 = find( MLRI_wR1 < MLRI_w2R1);
	MLRI_piR2 = find( MLRI_wR2 < MLRI_w2R2);
	  RI_piB1 = find(   RI_wB1 <   RI_w2B1);
	  RI_piB2 = find(   RI_wB2 <   RI_w2B2);
	MLRI_piB1 = find( MLRI_wB1 < MLRI_w2B1);
	MLRI_piB2 = find( MLRI_wB2 < MLRI_w2B2);

    % guide updating
	  RI_Guider1 = red  +   RI_R1;
	  RI_Guideb1 = blue +   RI_B1;
	  RI_Guider2 = red  +   RI_R2;
	  RI_Guideb2 = blue +   RI_B2;
	MLRI_Guider1 = red  + MLRI_R1;
	MLRI_Guideb1 = blue + MLRI_B1;
	MLRI_Guider2 = red  + MLRI_R2;
	MLRI_Guideb2 = blue + MLRI_B2;
	
    % select smallest iteration criteria at each pixel
	  RI_R1(  RI_piR1) =   RI_Guider1(  RI_piR1);
	MLRI_R1(MLRI_piR1) = MLRI_Guider1(MLRI_piR1);
	  RI_R2(  RI_piR2) =   RI_Guider2(  RI_piR2);
	MLRI_R2(MLRI_piR2) = MLRI_Guider2(MLRI_piR2);
	  RI_B1(  RI_piB1) =   RI_Guideb1(  RI_piB1);
	MLRI_B1(MLRI_piB1) = MLRI_Guideb1(MLRI_piB1);
	  RI_B2(  RI_piB2) =   RI_Guideb2(  RI_piB2);
	MLRI_B2(MLRI_piB2) = MLRI_Guideb2(MLRI_piB2);

    % update minimum iteration criteria
	  RI_w2R1(  RI_piR1) =   RI_wR1(  RI_piR1);
	  RI_w2R2(  RI_piR2) =   RI_wR2(  RI_piR2);
	  RI_w2B1(  RI_piB1) =   RI_wB1(  RI_piB1);
	  RI_w2B2(  RI_piB2) =   RI_wB2(  RI_piB2);
	MLRI_w2R1(MLRI_piR1) = MLRI_wR1(MLRI_piR1);
	MLRI_w2R2(MLRI_piR2) = MLRI_wR2(MLRI_piR2);
	MLRI_w2B1(MLRI_piB1) = MLRI_wB1(MLRI_piB1);
	MLRI_w2B2(MLRI_piB2) = MLRI_wB2(MLRI_piB2);

	% guided filter window size update
	h = h+1;
	v = v+1;
	h2 = h2+1;
	v2 = v2+1;
end

%% Step(iii): adaptive combining
% combining weight
  RI_w2R1 = 1./(  RI_w2R1+1e-10);
  RI_w2R2 = 1./(  RI_w2R2+1e-10);
MLRI_w2R1 = 1./(MLRI_w2R1+1e-10);
MLRI_w2R2 = 1./(MLRI_w2R2+1e-10);
  RI_w2B1 = 1./(  RI_w2B1+1e-10);
  RI_w2B2 = 1./(  RI_w2B2+1e-10);
MLRI_w2B1 = 1./(MLRI_w2B1+1e-10);
MLRI_w2B2 = 1./(MLRI_w2B2+1e-10);

wR = RI_w2R1+RI_w2R2+MLRI_w2R1+MLRI_w2R2;
wB = RI_w2B1+RI_w2B2+MLRI_w2B1+MLRI_w2B2;

% combining
red2  = (RI_w2R1.*RI_R1+RI_w2R2.*RI_R2+MLRI_w2R1.*MLRI_R1+MLRI_w2R2.*MLRI_R2)./(wR+1e-32);
blue2 = (RI_w2B1.*RI_B1+RI_w2B2.*RI_B2+MLRI_w2B1.*MLRI_B1+MLRI_w2B2.*MLRI_B2)./(wB+1e-32);

% output of the second step
red  =  red + red2.*mask(:,:,2);
blue = blue + blue2.*mask(:,:,2);

red  = clip( red, 0, DynamicRange);
blue = clip(blue, 0, DynamicRange);

end
