% Radiating incident field
%
%  Warning: this is an ABSTRACT base class... it is not possible
%  for objects of this class to be instantiated. 
%
%  Note: in future we will implement +, -, * for radiating_incident but
%  only uminus is currently implemented.
%
% See also: incident, point_source,
% radiatingwavefunction2d.
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


classdef radiating_incident < incident
    
    properties
    end
    
    %=================================================================
    % methods with standard access
    %=================================================================

    methods
        
        %-------------------------------------------------
        % constructor
        %-------------------------------------------------

        function self = radiating_incident(varargin)
            
            % simply call parent constructor
            self = self@incident(varargin{:});
            
        end
        
        %-------------------------------------------------
        % overloaded -
        %-------------------------------------------------
                
        function obj = uminus(self)
            
            % check that input object is of class incident
            if ~isa(self,'incident')
                error('self must be of class incident')
            end
        
            % create instance of incidentminus class
            obj = radiatingincidentuminus(self);
            
        end
        
    end
    
    %=================================================================
    % abstract methods
    %=================================================================

    % These must be provided by child classes

    methods(Abstract)
        
        % compute the far field
        val = evaluateFarField(self,points);
        
    end
    
end