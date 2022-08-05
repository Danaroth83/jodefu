function v = cell_averaging(u,z)
%CELL_AVERAGING average the 2D signal u over square zxz sized cells

[ny,nx] = size(u);
nx0 = nx/z; 
ny0 = ny/z; 
if ((nx0 ~=floor(nx0)) || (ny0 ~= floor(ny0)))
    error('the dimensions of u must be proportional to z'); 
end
idx = 1+z*(0:(nx0-1));
idy = 1+z*(0:(ny0-1));
v = zeros(ny0,nx0);
for dx = 0:(z-1)
    for dy = 0:(z-1)
        v = v + u(dy+idy,dx+idx);
    end
end
v = v/z^2;
end

