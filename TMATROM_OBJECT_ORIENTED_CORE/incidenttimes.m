% Incident field that is a scalar multiple of another incident field.
%
%   Note: this is intended to be used only by incident.mtimes method.
%
% See also: point_source, plane_wave, incident.
%
% Stuart C. Hawkins - 30 November 2018

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


classdef incidenttimes < incident
    
    properties
        left
        scalar
        kwave
    end
    
    methods
        
        %------------------------------------------------
        % constructor
        %------------------------------------------------
        
        function self = incidenttimes(left,scalar)
        
            % check that input object is of class incident            
            if ~isa(left,'incident')
                error('left must be of class incident')
            end
            
            % check that input object is of class double
            if ~isa(scalar,'double')
                error('scalar must be of class double')
            end
            
            % store given objects... we will use these when we need to
            % evaluate etc
            self.kwave = left.kwave;
            self.left = left;
            self.scalar = scalar;
            
        end
   
        %------------------------------------------------
        % get coefficients
        %------------------------------------------------

        function cof = get_coefficients(self,centre,nmax)
           
            % combine output of left and right object methods
            cof = self.scalar * self.left.get_coefficients(centre,nmax);
            
        end
    
        %------------------------------------------------
        % evaluate
        %------------------------------------------------

        function val = evaluate(self,varargin)
            
            % combine output of left and right object methods
            val = self.scalar * self.left.evaluate(varargin{:});
            
        end
        
        %------------------------------------------------
        % evaluate gradient
        %------------------------------------------------

        function [dx,dy] = evaluateGradient(self,varargin)
           
            % use left and right object methods
            [dx,dy] = self.left.evaluateGradient(varargin{:});
                      
            % do scalar multiplication
            dx = self.scalar * dx;
            dy = self.scalar * dy;
            
        end
                        
    end
    
end
    