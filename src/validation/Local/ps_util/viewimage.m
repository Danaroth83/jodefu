function immagine = viewimage(immagine,tol1,tol2,tol3)
% viewimage(immagine,tol1)
% visualizzazione [3-2-1] di immagine a 3 bande con stretching lineare con saturazione
% fissata da tol1, esempio
% [0.01 0.99], cioe' 2% (default), uguale per le tre bande 
% oppure 
% viewimage(immagine,tol1,tol2,tol3)
% con tol1, tol2, tol3 diversi per le tre bande 
% con le tre coppie di valori in percentuale oppure in livello assoluto 
iptsetpref('ImshowBorder', 'tight')
immagine = double(immagine);
L=size(immagine,3);
if (L<3)
    immagine=immagine(:,:,[1 1 1]);
end

if nargin == 1
    tol1 = [0.01 0.99];
end
if nargin <= 2
    tol = [tol1;tol1;tol1];
    immagine = linstretch(immagine,tol);
%     figure,
    imshow(immagine(:,:,3:-1:1),[])
elseif nargin == 4
    if sum(tol1(2)+tol2(2)+tol3(2)) <= 3
        tol = [tol1;tol2;tol3];
        immagine = linstretch(immagine,tol);
%         figure,
        imshow(immagine(:,:,3:-1:1),[])
    else
        tol = [tol1;tol2;tol3];
        [N M L] = size(immagine);
        NM = N*M;
        for i=1:3
            b = reshape(double(uint16(immagine(:,:,i))),NM,1);
            b(b<tol(i,1))=tol(i,1);
            b(b>tol(i,2))=tol(i,2);
            b = (b-tol(i,1))/(tol(i,2)-tol(i,1));
            immagine(:,:,i) = reshape(b,N,M);
        end
%         figure,
        imshow(immagine(:,:,3:-1:1),[])
    end
end
iptsetpref('ImshowBorder', 'loose')

function immagine=linstretch(immagine,tol)
    [N M L] = size(immagine);
    NM = N*M;
    for i=1:3
        b = reshape(double(uint16(immagine(:,:,i))),NM,1);
        [hb,levelb] = hist(b,max(b)-min(b));
        chb = cumsum(hb);
        t(1)=ceil(levelb(find(chb>NM*tol(i,1), 1 )));
        t(2)=ceil(levelb(find(chb<NM*tol(i,2), 1, 'last' )));
        b(b<t(1))=t(1);
        b(b>t(2))=t(2);
        b = (b-t(1))/(t(2)-t(1));
        immagine(:,:,i) = reshape(b,N,M);
    end

    