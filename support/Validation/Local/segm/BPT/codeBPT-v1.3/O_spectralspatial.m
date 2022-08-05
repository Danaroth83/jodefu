function dist = O_spectralspatial(Mi,Mj)

% Computes the spectral-spatial similarity between two sets of endmembers
% and their mean abundances

% average abundances
Ai = Mi.A;
[r,c,p] = size(Ai);
Ai = reshape(Ai,r*c,p);
Ai_avg = mean(Ai);
Aj = Mj.A;
[r,c,p] = size(Aj);
Aj = reshape(Aj,r*c,p);
Aj_avg = mean(Aj);

[dist,~] = disimilitud_ab(Mi.E',Mj.E',Ai_avg,Aj_avg,'sam');