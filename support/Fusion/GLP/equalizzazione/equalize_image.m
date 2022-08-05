function v1 = equalize_image (v1, v2)
% equalizza v1 rispetto a v2, tre diversi metodi

% equalizzazione eq_cdf
v1 = eq_cdf (v1, v2, 1000);

% equalizzazione histeq
% hgram = imhist (v2);
% v1 = histeq (v1, hgram);

% equalizzazione media-varianza
% v1 = (v1 - mean(v1(:))) / std(v1(:)) * std(v2(:)) + mean(v2(:));

end

