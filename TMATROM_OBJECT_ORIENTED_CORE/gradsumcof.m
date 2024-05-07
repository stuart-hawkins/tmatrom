% Evaluate gradient of wavefunction expansion series.
%
%   [dx,dy] = gradsumcof(x,x0,k,c,'H') returns the gradient (dx,dy) of the
%   radiating wavefunction expansion with coefficients c, centre x0 and
%   wavenumber k at points x.
%
%   [dx,dy] = gradsumcof(x,x0,k,c,'J') returns the gradient (dx,dy) of the
%   regular wavefunction expansion with coefficients c, centre x0 and
%   wavenumber k at points x.
%
%   [dx,dy] = gradsumcof(x,x0,k,c,'F') returns the gradient (dx,dy) of the
%   far field of the radiating wavefunction expansion with coefficients c,
%   centre x0 and wavenumber k at points abs(x) on the unit circle.
%
% Note: in the above vectors in the plane are represented by
% complex numbers.
%
% Stuart C. Hawkins - 7 May 2024

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


function [dx,dy] = gradsumcof(points,centre,kwave,cof,type)

% use dsumcof to evaluate the partial derivatives of the expansion with
% respect to polar coordinates r and theta
[dr,dtheta,er,etheta,rad] = dsumcof(points,centre,kwave,cof,type);

% calculate the partial derivatives with respect to x and y
dx = real(er) .* dr + real(etheta)./rad .* dtheta;
dy = imag(er) .* dr + imag(etheta)./rad .* dtheta;

