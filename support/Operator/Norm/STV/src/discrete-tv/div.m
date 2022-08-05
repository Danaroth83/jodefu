function d = div(px,py)
% DIV 
% compute d = discrete divergence (= opposite of the adjoint of the 'grad'
% operator) of the 2D field vector (px,py)

[ny,nx] = size(px); 
if (size(py,1) ~= ny) || (size(py,2) ~= nx)
    error('The two components of the input 2D vector field must have the same dimensions'); 
end
div_x = px - px(:,[1,1:(nx-1)]); 
div_x(:,1) = px(:,1); 
div_x(:,nx) = -px(:,nx-1); 
div_y = py - py([1,1:(ny-1)],:); 
div_y(1,:) = py(1,:); 
div_y(ny,:) = -py(ny-1,:); 
d = div_x + div_y; 
end

