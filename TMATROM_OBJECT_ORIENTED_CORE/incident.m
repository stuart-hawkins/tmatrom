% Incident field
%
%  Warning: this is an ABSTRACT base class... it is not possible
%  for objects of this class to be instantiated. 
%
% See also: plane_wave, point_source, regularwavefunction2d,
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


classdef incident < handle
    
    properties
    end
    
    %=================================================================
    % methods with standard access
    %=================================================================

    methods
    
        %-------------------------------------------------
        % overloaded +
        %-------------------------------------------------
        
        function obj = plus(self,other)
            
            % check that input object is of class incident
            if ~isa(self,'incident')
                error('self must be of class incident')
            end
        
            % check that input object is of class incident
            if ~isa(other,'incident')
                error('other must be of class incident')
            end

            % create instance of incidentplus class
            obj = incidentplus(self,other);
            
        end
    
        %-------------------------------------------------
        % overloaded -
        %-------------------------------------------------
        
        function obj = minus(self,other)
            
            % check that input object is of class incident
            if ~isa(self,'incident')
                error('self must be of class incident')
            end
        
            % check that input object is of class incident
            if ~isa(other,'incident')
                error('other must be of class incident')
            end

            % create instance of incidentminus class
            obj = incidentminus(self,other);
            
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
            obj = incidentuminus(self);
            
        end
        
        %-------------------------------------------------
        % overloaded *
        %-------------------------------------------------

        % This provides facility to multiply by a scalar
        
        function obj = mtimes(obj1,obj2)
            
            % check that one of the input objects is of class incident
            if ~isa(obj1,'incident') && ~isa(obj2,'incident')
                % This should never happen... this method will only be
                % called if one of the arguments is of class incident
                error('one of arguments must be of class incident')
            end
        
            % check that one of the input object is a double
            if ~isa(obj1,'double') && ~isa(obj2,'double')
                error('one of the arguments must be of class double')
            end

            % create instance of incidenttimes class
            if isa(obj1,'incident')
                obj = incidenttimes(obj1,obj2);
            else
                obj = incidenttimes(obj2,obj1);
            end
            
        end
        
    end
    
    %=================================================================
    % abstract methods
    %=================================================================

    % These must be provided by child classes

    methods(Abstract=true)
        
        cof = get_coefficients(self,centre,nmax);
        val = evaluate(self,points,mask);
        [dx,dy] = evaluateGradient(self,points,mask);
                
    end
    
end