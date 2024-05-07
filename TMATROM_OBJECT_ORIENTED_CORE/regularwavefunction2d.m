% Regular circular wavefunction
%
%  u = regularwavefunction2d(n,k,x) returns a regular wavefunction object u
%  with order n, wavenumber k and origin x.
%
% Also:
%
%   f = u.evaluate(z) returns the values f of the wavefunction at points z.
%
%   f = u.evaluate(z,mask) returns the values f of the wavefunction at
%   points z for which mask==1 and NaN elsewhere.
%
%   [dx,dy] = u.evaluateGradient(z) returns dx and dy the partial 
%   derivatives of the wavefunction in the x and y directions respectively
%   at the points z.
%
%   [dx,dy] = u.evaluateGradient(z,mask) returns dx and dy the partial 
%   derivatives of the wavefunction in the x and y directions respectively
%   at the points z for which mask==1 and NaN elsewhere.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: wavefunction2d, radiatingwavefunction2d, incident.
%
% Stuart C. Hawkins - 9 January 2023

% Copyright 2014, 2015, 2016, 2017, 2018, 2022, 2023, 2024 Stuart C. Hawkins and M. Ganesh.
% 	
% This file is part of TMATROM.
% 
% TMATROM is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATROM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TMATROM.  If not, see <http://www.gnu.org/licenses/>.


classdef regularwavefunction2d < wavefunction2d
    
    properties
    end
    
    methods

        %-------------------------------------------------
        % constructor
        %-------------------------------------------------
        
        function self = regularwavefunction2d(order,kwave,origin)

            % set default for origin
            if nargin < 3
                origin = 0;
            end

            % call parent constructor
            self = self@wavefunction2d(order,kwave,origin);
                       
        end
        
        %-----------------------------------------------
        % function that determines radial behaviour
        %-----------------------------------------------

        function val = radial_function(self,r)

            val = besselj(abs(self.order),r);
            
        end
        
        %-----------------------------------------------
        % derivative of function that determines radial behaviour
        %-----------------------------------------------

        function val = derivative_radial_function(self,r)
            
            val = besseljd(abs(self.order),r);
            
        end
        
        %-----------------------------------------------
        % return vector of coefficients for scatterer at
        % centre
        %-----------------------------------------------

        function cof = get_coefficients(self,centre,nmax)
    
            % check centre is self.origin... other centres are not currently
            % supported because they need translation addition theorem etc
            if centre~=self.origin                
                error('Only centre = self.origin currently supported.')                
            end

            % setup return vector
            cof = zeros(2*nmax+1,1);

            if abs(self.order) <= nmax
                % set entry corresponding to current order to 1
                cof(self.order+nmax+1) = 1;
            else
                % the order is larger than nmax and so the coefficient
                % vector does not have a place to put the 1
                warning('order > nmax so coefficient vector is incomplete.')
            end                           
            
        end
        
    end % end methods
    
end