% Create T-matrix using solver
%
%   T = tmatrix(n,k,s) returns the T-matrix of order n, wavenumber k, and
%   origin 0. The scatterer details are contained in the object S which
%   must be an instance of the 'solver' class.
%
%   T = tmatrix(n,k,s,x) returns the T-matrix of order n, wavenumber k, and
%   origin x. 
%
% See also: tmatrix, regularwavefunctionexpansion, 
% radiatingwavefunctionexpansion, plane_wave, point_source.
%
% Stuart C. Hawkins - 6 March 2018

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


function tmat = ghtmatrix(order,kwave,solver,origin)

% set default for origin
if nargin < 4
    origin = 0;
end

% check that solver is of the correct type
if ~isa(solver,'solver')
    error('solver must be an instance of the solver class')
end

% - - - - - - - - - - - - - - - - -
% setup the quadrature points
% - - - - - - - - - - - - - - - - -

% get quadrature points and weights
points = pi*(0:2*order+1)/(order+1);
weight = pi/(order+1);

% ensure these are column vector
points = points(:);

% - - - - - - - - - - - - - - - - -
% solve the scattering problems
% - - - - - - - - - - - - - - - - -

% setup a cell array of incident fields
for n = -order:order
    inc{n + order + 1} = regularwavefunction2d(n, ...
        kwave,origin);
end

% set the incident field in the solver
solver.setIncidentField(inc);

% solve the scattering problems
solver.solve();

% get the farfield
farfield = solver.getFarField(points,1:2*order+1);

% assume the farfield was computed with origin 0... we need
% to adjust for the modified origin
if origin ~= 0
    
    % Use (21) in Dufva et al, Progress in Electromagnetics
    % Research B, Vol 4, 79-99, 2008
    
    % compute scaling factor
    sigma = exp(1i*kwave*real(conj(origin)*exp(1i*points)));
    
    % transform the far field
    farfield = spdiags(sigma(:),0,length(sigma),length(sigma)) * farfield;
    
end

% - - - - - - - - - - - - - - - - -
% compute the T-matrix
% - - - - - - - - - - - - - - - - -

% get the order of the wavefunctions in a row vector
n = (-order:order);

% compute the T-matrix using (16) in Ganesh and Hawkins
% ANZIAM J. 51 C215--C230 (2010)
matrix = weight * scaling(n,kwave) * exp(1i*points*n)' * farfield;

% put some basic information in the comment string
comment = sprintf('Computed using ghtmatrix v%d with solver %s.',tmatrom_version(),class(solver));

% creat the T-matrix object
tmat = tmatrix(order,kwave,matrix,origin,comment);

end

%-----------------------------------------
% scaling function
%-----------------------------------------

% gives the coefficient in (16) Ganesh and Hawkins
% ANZIAM J. 51 C215--C230 (2010)

function val = scaling(n,kwave)

vec = 0.25*(1+1i)*sqrt(kwave/pi)*1i.^abs(n);

val = spdiags(vec(:),0,length(vec),length(vec));

end
