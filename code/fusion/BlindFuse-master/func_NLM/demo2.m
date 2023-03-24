clear
clc
clf
colormap(gray)

ima=100*ones(100);
ima(50:100,:)=50;
ima(:,50:100)=2*ima(:,50:100);
fs=fspecial('average');
ima=imfilter(ima,fs,'symmetric');

% % create example image
% ima=100*ones(50);
% ima(25:50,:)=50;
% ima(:,25:50)=2*ima(:,25:50);
% fs=fspecial('average');
% ima=imfilter(ima,fs,'symmetric');

% add some noise
sigma=10;
rima=ima+sigma*randn(size(ima));

% show it
imagesc(rima)
drawnow

% denoise it
[fima,patches]=NLM_personal(rima,11,1,2500,sigma);
fima2=NLmeansfilter(rima,5,2,sigma);

% show results
clf
% subplot(2,2,1),imagesc(ima),title('original');
% subplot(2,2,2),imagesc(rima),title('noisy');
% subplot(2,2,3),imagesc(fima),title('filtered');
% subplot(2,2,4),imagesc(rima-fima),title('residuals');

imagesc(fima);
imagesc(fima2);


a=exchangepatches(fima,[50,50],patches,5,1);
