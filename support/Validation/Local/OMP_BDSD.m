function [a, indx,D1] = OMP_BDSD (D, y, delta, num_atoms)
% OMP new adpted for BDSD
% D = dictionary (matrix)
% y = observation (column vector)
% delta = maximum percentage of the relative admissible value
% keyboard
[L_atom, N_atom] = size (D);
a = zeros (num_atoms, 1);
res = y;
indx = num_atoms*ones (num_atoms, 1);
curr_delta = 100*norm(res)/norm(y);
% S = sum(D.*D).^0.5;
% D1=D./repmat(S,size(D,1),1);
D1=D;
j = 0;
while  j < num_atoms && j<N_atom && curr_delta > delta
    j = j+1;
    proj = D1' * res;
    [tmp, imax] = max(abs(proj));
    imax = imax(1);
    indx(j) = imax;
    Da = zeros (L_atom, 1);
    Di = D1(:,indx(1:j));
    yi = y;
    a(1:j) = pinv(Di) * yi;
    Da= Di * a(1:j);
    res = y - Da;
    curr_delta = 100*norm(res)/norm(y);
end



