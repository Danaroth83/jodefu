function green = green_interpolation(mosaic, mask, pattern, sigma, DynamicRange)

% mask
imask = (mask == 0);

% raw CFA data
rawq = sum(mosaic, 3);

%%% Calculate Horizontal and Vertical Color Differences %%%
% mask
maskGr = zeros(size(rawq));
maskGb = zeros(size(rawq));

if strcmp(pattern,'grbg')
    maskGr(1:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
    maskGb(2:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
elseif strcmp(pattern,'rggb')
    maskGr(1:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
    maskGb(2:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
elseif strcmp(pattern,'gbrg')
    maskGb(1:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
    maskGr(2:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
elseif strcmp(pattern,'bggr')
    maskGb(1:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
    maskGr(2:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
end

% guide image
Kh = [1/2,0,1/2];
Kv = Kh';
rawh = imfilter(rawq, Kh, 'replicate');
rawv = imfilter(rawq, Kv, 'replicate');

Guidegh = mosaic(:,:,2) + rawh.*mask(:,:,1) + rawh.*mask(:,:,3);
Guiderh = mosaic(:,:,1) + rawh.*maskGr;
Guidebh = mosaic(:,:,3) + rawh.*maskGb;
Guidegv = mosaic(:,:,2) + rawv.*mask(:,:,1) + rawv.*mask(:,:,3);
Guiderv = mosaic(:,:,1) + rawv.*maskGb;
Guidebv = mosaic(:,:,3) + rawv.*maskGr;

% tentative image
h = 5;
v = 0;
eps = 0;
tentativeRh = guidedfilter(Guidegh, mosaic(:,:,1), mask(:,:,1), h, v, eps);
tentativeGrh = guidedfilter(Guiderh, mosaic(:,:,2).*maskGr, maskGr, h, v, eps);
tentativeGbh = guidedfilter(Guidebh, mosaic(:,:,2).*maskGb, maskGb, h, v, eps);
tentativeBh = guidedfilter(Guidegh, mosaic(:,:,3), mask(:,:,3), h, v, eps);
tentativeRv = guidedfilter(Guidegv, mosaic(:,:,1), mask(:,:,1), v, h, eps);
tentativeGrv = guidedfilter(Guiderv, mosaic(:,:,2).*maskGb, maskGb, v, h, eps);
tentativeGbv = guidedfilter(Guidebv, mosaic(:,:,2).*maskGr, maskGr, v, h, eps);
tentativeBv = guidedfilter(Guidegv, mosaic(:,:,3), mask(:,:,3), v, h, eps);

% residual
residualGrh = ( mosaic(:,:,2) - tentativeGrh ) .* maskGr;
residualGbh = ( mosaic(:,:,2) - tentativeGbh ) .* maskGb;
residualRh = ( mosaic(:,:,1) - tentativeRh ) .* mask(:,:,1);
residualBh = ( mosaic(:,:,3) - tentativeBh ) .* mask(:,:,3);
residualGrv = ( mosaic(:,:,2) - tentativeGrv ) .* maskGb;
residualGbv = ( mosaic(:,:,2) - tentativeGbv ) .* maskGr;
residualRv = ( mosaic(:,:,1) - tentativeRv ) .* mask(:,:,1);
residualBv = ( mosaic(:,:,3) - tentativeBv ) .* mask(:,:,3);

% residual interpolation
Kh = [1/2,0,1/2];
residualGrh = imfilter(residualGrh, Kh, 'replicate');
residualGbh = imfilter(residualGbh, Kh, 'replicate');
residualRh = imfilter(residualRh, Kh, 'replicate');
residualBh = imfilter(residualBh, Kh, 'replicate');

Kv = Kh';
residualGrv = imfilter(residualGrv, Kv, 'replicate');
residualGbv = imfilter(residualGbv, Kv, 'replicate');
residualRv = imfilter(residualRv, Kv, 'replicate');
residualBv = imfilter(residualBv, Kv, 'replicate');

% add tentative image
Grh = ( tentativeGrh + residualGrh ) .* mask(:,:,1);
Gbh = ( tentativeGbh + residualGbh ) .* mask(:,:,3);
Rh = ( tentativeRh + residualRh ) .* maskGr;
Bh = ( tentativeBh + residualBh ) .* maskGb;
Grv = ( tentativeGrv + residualGrv ) .* mask(:,:,1);
Gbv = ( tentativeGbv + residualGbv ) .* mask(:,:,3);
Rv = ( tentativeRv + residualRv ) .* maskGb;
Bv = ( tentativeBv + residualBv ) .* maskGr;

% vertical and horizontal color difference 
difh = mosaic(:,:,2) + Grh + Gbh - mosaic(:,:,1) - mosaic(:,:,3) - Rh -Bh;
difv = mosaic(:,:,2) + Grv + Gbv - mosaic(:,:,1) - mosaic(:,:,3) - Rv -Bv;

%%% Combine Vertical and Horizontal Color Differences %%%
% color difference gradient
Kh = [1,0,-1];
Kv = Kh';
difh2 = abs(imfilter(difh, Kh, 'replicate'));
difv2 = abs(imfilter(difv, Kv, 'replicate'));

% directional weight
K = ones(5,5);
wh = imfilter(difh2, K, 'replicate');
wv = imfilter(difv2, K, 'replicate');
Kw = [1,0,0,0,0]; 
Ke = [0,0,0,0,1];
Ks = Ke'; 
Kn = Kw';
Ww = imfilter(wh, Kw, 'replicate');
We = imfilter(wh, Ke, 'replicate');
Wn = imfilter(wv, Kn, 'replicate');
Ws = imfilter(wv, Ks, 'replicate');
Ww = 1 ./ (Ww.*Ww + 1E-32);
We = 1 ./ (We.*We + 1E-32);
Ws = 1 ./ (Ws.*Ws + 1E-32);
Wn = 1 ./ (Wn.*Wn + 1E-32);

% combine directional color differences
h = fspecial('gaussian', [1,9], sigma);
Ke = [0,0,0,0,1,1,1,1,1] .* h; 
Kw = [1,1,1,1,1,0,0,0,0] .* h;
Ke = Ke / sum(Ke, 2);
Kw = Kw / sum(Kw, 2);
Ks = Ke'; 
Kn = Kw';
difn = imfilter(difv, Kn, 'replicate');
difs = imfilter(difv, Ks, 'replicate');
difw = imfilter(difh, Kw, 'replicate');
dife = imfilter(difh, Ke, 'replicate');
Wt = Ww + We + Wn + Ws;
dif = (Wn.*difn + Ws.*difs + Ww.*difw + We.*dife) ./ Wt ;

%%% Calculate Green by adding bayer raw data %%% 
green = dif + rawq;
green = green.*imask(:,:,2) + rawq.*mask(:,:,2);

% clip to 0-255
green = clip(green, 0, DynamicRange);

end
