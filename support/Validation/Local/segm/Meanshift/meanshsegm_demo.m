% MEANSHSEGM_DEMO Demo showing the usage of meanshsegm 
% CMP Vision algorithms http://visionbook.felk.cvut.cz
%
% Example
%
% Mean shift segmentation is applied to an RGB color image (264x
% 512 pixels) which takes several minutes.  It is possible to use
% a different color space or a grayscale image. Small regions can be
% eliminated by post-processing using remsmall.  

ImageDir='images/';%directory containing the images
addpath('..') ;
cmpviapath('..') ;
if (exist('output_images')~=7)
  mkdir('output_images');
end




if 1
img=imresize(imread([ImageDir 'spiderman.jpg']),0.1) ;
l=meanshsegm(img,10,30) ;

figure(1) ;
imagesc(img); % title('input image');
axis image ; axis off ; 
exportfig(gcf,'output_images/meanshsegm_input.eps') ;

figure(2) ; 
imagesc(label2rgb(l-1,'jet','w','shuffle')) ; 
axis image ; axis off ; 
exportfig(gcf,'output_images/meanshsegm_output.eps') ;
end

