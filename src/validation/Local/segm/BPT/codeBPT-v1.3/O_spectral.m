function dist = O_spectral(Mi,Mj)

% Computes the spectral-spatial similarity between two sets of endmembers

[dist,~] = disimilitud(Mi.E,Mj.E,'sam','grana');