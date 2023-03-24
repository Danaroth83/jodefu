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


% Multiquadratic

classdef mqx < rbfx
    methods
        function obj = mqx()  % constructor
          obj@rbfx();  % call constructor of the superclass 
        end
       
        function v = rbf(obj,r,s)
            v = (1 + (s.*r).^2 ).^(1/2);
        end
       
        function d = D1(obj,r,s,x) 
            d= (x.*s.^2)./( 1 + (s.*r).^2).^(1/2);
        end
       
        function d = D2(obj,r,s,x)
            d = s.^2.*(1+ s.^2.*r.^2-s.^2.*x.^2)./( 1 + (s.*r).^2).^(3/2);
        end
        
        function d = D3(obj, r, s, x) 
            d = 3*x.*s.^4.*(1+ s.^2.*r.^2-s.^2.*x.^2)./( 1 + (s.*r).^2).^(5/2);
        end
        
        function d = D4(obj, r, s, x)
            d=3*s.^4.*(s.^4.*(r.^2-x.^2).*(5*x.^2 - r.^2) + 2*s.^2.*(3*x.^2 - r.^2) - 1)./( 1 + (s.*r).^2).^(7/2);
        end
        
        function d = G(obj, r, s, x, y)   % Gradient
           d= ((x+y).*s.^2)./( 1 + (s.*r).^2).^(1/2);
        end
        
        function d = L(obj, r, s)         % Laplacian
           d= s.^2.*(s.^2.*r.^2 + 2)./( 1 + (s.*r).^2).^(3/2);
        end
        
         % x and y not used but required by the abstract function definition in the superclass
        function d = B(obj, r, s, ~, ~)   % Biharmonic operator   
            d= s.^4.*(s.^4.*r.^2 + 8*s.^2.*r.^2 - 8)./( 1 + (s.*r).^2).^(7/2);
            % d= 3*s.^4.*(3*s.^4.*r.^4-24*s.^2.*r.^2+8)./(1 + (s.*r).^2 ).^(9/2);
        end

     
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x ) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y)
           d= s.^4.*x.*(-s.^2.*x^2 + 2*s.^2.*y.^2 - 1)./(1 + (s.*r).^2 ).^(5/2);
       end
       
       function d = D22(obj, r, s, x, y)
           d=s.^4.*(2*s.^4.*x.^4 - 11*s.^4.*x.^2.*y^2 + 2*s.^4.*y.^4 + s.^2.*x.^2 + s.^2.*y.^2 - 1)./(1 + (s.*r).^2 ).^(7/2);
           % d = -(12*s.^6.*y.^2)./(s.^2.*r.^2 + 1).^(5/2) + (2*s.^4)./(s.^2.*r.^2 + 1).^(3/2) + (s.^4.*y.^2 + s.^2).* ((15*s.^4.*y.^2)./(s.^2.*r.^2 + 1).^(7/2) - (3*s.^2)./(s.^2.*r.^2 + 1).^(5/2));
           %d = -8*s.^4.*(-1 + 4*r.^2.*s.^2 + 5*(x.^4 + y.^4).*s.^4 - 38*x.^2.*y.^2.*s.^4  )./(1 + (s.*r).^2 ).^5;
       end
      
   end % methods
end  % class
