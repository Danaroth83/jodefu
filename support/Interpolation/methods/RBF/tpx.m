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


% Thin Plate

classdef tpx < rbfx
    methods
        function obj = tpx()  % constructor
          obj@rbfx();  % call constructor of the superclass 
        end
       
        function v = rbf(obj,r,s)
            v = s.^2.*r.*log((s.*r).^r);
            % v = (s.*r).^2.*log(s.*r);
        end
       
        function d = D1(obj,r,s,x) 
            d = s.^2.*x.*(2*log((s.*r).^r)./r+1);
            % d = s.^2.*x.*(2*log(s.*r)+1);
        end
       
        function d = D2(obj,r,s,x)
            d = s.^2.*(2*r.*log((s.*r).^r) + 2*x.^2 + r.^2)./r.^2;
            % d = s.^2.*(2*r.^2.*log(s.*r) + 2*x.^2 + r.^2)./r.^2;
        end
        
        function d = D3(obj,r,s,x) 
            d = 2*s.^2.*x.*(3.*r.^2-2*x.^2)./r.^4;
        end
        
        function d = D4(obj,r,s,x)
            d = -2*s.^2.*(4*x.^4-3*r.^4)./r.^6;
        end
        
        function d = G(obj, r, s, x, y)   % Gradient
            d = s.^2.*(x+y).*(2*log((s.*r).^r)./r+1);
            % d = s.^2.*x.*(2*log(s.*r)+1);
        end
        
        function d = L(obj, r, s)         % Laplacian
            d = 4*s.^2.*(log((s.*r).^r)./r+1);
            % d = 4*s.^2.*(log(s.*r)+1);
        end
        
         % x and y not used but required by the abstract function definition in the superclass
        function d = B(obj, r, ~, ~, ~)   % Biharmonic operator   
            d= zeros(size(r));
            % d= 3*s.^4.*(3*s.^4.*r.^4-24*s.^2.*r.^2+8)./(1 + (s.*r).^2 ).^(9/2);
        end

     
% D12       
% mixed partial derivative
%      D_{xyy}  d = D12( r, s, x, y ) 
% or   D_{yxx}  d = D12( r, s, y, x ) 
% depending on the order of the x and y arguments
       
       function d = D12(obj, r, s, x, y)
           d= 2*s.^2.*x.*(x.^2-y.^2)./r.^4;
       end
       
       function d = D22(obj, r, s, x, y)
           d=-2*s.^2.*(r.^4 - 8*x.^2.*y.^2)./r.^6;
       end
      
   end % methods
end  % class
