%     Matlab Radial Basis Function Toolkit (MRBFT)
% 
%     Copyright (c) 2020 Daniele Picone
% 
%     Based on the original MRBFT Toolbox     
%
%     GNU General Public License ("GPL") copyright permissions statement:
%     **************************************************************************
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Inverse Multiquadratic

classdef imqx < rbfx
    methods
        function obj = imqx()  % constructor
          obj@rbfx();  % call constructor of the superclass 
        end
       
        function v = rbf(obj,r,s), v = 1./(1 + (s.*r).^2 ).^(1/2); end
       
        function d = D1(obj,r,s,x), d = -(x.*s.^2)./(1 + (s.*r).^2 ).^(3/2); end
       
        function d = D2(obj,r,s,x)
            d = s.^2.*(- r.^2.*s.^2 + 3*x.^2.*s.^2 -1)./(1 + (s.*r).^2 ).^(5/2);
            % d = 3*s.^4.*(-r.^2.*s.^2 + 4*s.^2*x.^2 - 1)./(r.^2*s.^2 + 1).^3;
        end
        
        function d = D3(obj, r, s, x) 
            d = s.^4.*(9*s.^2.*r.^2.*x-15*s.^2.*x.^3+9*x)./(1.0 + (s.*r).^2 ).^(7/2);
        end
        
        function d = D4(obj, r, s, x)
            d = (105*s.^8.*x.^4)./(s.^2.*r.^2+1).^(9/2)-(90*s.^6.*x.^2)./(s.^2.*r.^2 + 1).^(7/2) + (9*s.^4)./(s.^2.*r.^2 + 1).^(5/2);
        end
        
        function d = G(obj, r, s, x, y)   % Gradient
           d = -(s.^2.*(x + y))./(1 + (s.*r).^2 ).^(3/2);
        end
        
        function d = L(obj, r, s)         % Laplacian
           d = s.^2.*(s.^2.*r.^2-2)./(1 + (s.*r).^2 ).^(5/2);
           % d = 4*s.^2.*(r.^2.*s.^2 - 1)./(1 + (s.*r).^2 ).^3;
           % d = 3*s.^4.*(-r.^2.*s.^2 + 4*s.^2*x.^2 - 1)./(r.^2*s.^2 + 1).^3;
        end
        
         % x and y not used but required by the abstract function definition in the superclass
        function d = B(obj, r, s, ~, ~)   % Biharmonic operator   
            d= 3*s.^4.*(3*s.^4.*r.^4-24*s.^2.*r.^2+8)./(1 + (s.*r).^2 ).^(9/2);
            %d = 64*( s.^4 - 4*r.^2.*s.^6 + r.^4.*s.^8  )./(1 + (s.*r).^2 ).^5;
            % d= 64*(s.^4 -4*r.^2.*s.^6 + r.^4.*s.^8)
        end

     
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x ) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y)
           d= -x.*s.^2.*(12*y.^2.*s.^4-3*x.^2.*s.^4-3*s.^2)./(1 + (s.*r).^2 ).^(7/2);
       end
       
       function d = D22(obj, r, s, x, y)
           d=-3*s.^4.*(4*s.^4.*x.^4 - 27*s.^4.*x.^2.*y.^2 + 4*s.^4.*y.^4 + 3*s.^2.*x.^2 + 3*s.^2.*y.^2 - 1)./(1 + (s.*r).^2 ).^(9/2);
           % d = -(12*s.^6.*y.^2)./(s.^2.*r.^2 + 1).^(5/2) + (2*s.^4)./(s.^2.*r.^2 + 1).^(3/2) + (s.^4.*y.^2 + s.^2).* ((15*s.^4.*y.^2)./(s.^2.*r.^2 + 1).^(7/2) - (3*s.^2)./(s.^2.*r.^2 + 1).^(5/2));
       end
      
   end % methods
end  % class
