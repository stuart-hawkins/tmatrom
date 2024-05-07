% Evaluate partial derivatives of wavefunction expansion series.
%
%   [dr,dtheta,er,etheta,rad] = dsumcof(x,x0,k,c,'H') returns the [partial
%   derivatives dr and dtheta of the radiating wavefunction expansion with
%   coefficients c, centre x0 and wavenumber k at points x. The unit
%   vectors er and etheta associated with dr and dtheta are also computed.
%
%   [...] = dsumcof(x,x0,k,c,'J') returns the values z of the regular 
%   wavefunction expansion with coefficients c, centre x0 and wavenumber k 
%   at points x.
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


function [dr,dtheta,er,etheta,rad] = dsumcof(points,centre,kwave,cof,type)

if strcmp(type,'F')
    error('Type F not supported')
end

%-------------------------------------------------
% setup
%-------------------------------------------------

% make sure the coefficient vector is a column vector
cof=cof(:);

% determine the maximum order from the length of the coefficient vector
nmax=0.5*(length(cof)-1);

% create a vector of indexes.... helps to vectorize the computation
n=-nmax:nmax;

%-------------------------------------------------
% turn points into a vector
%-------------------------------------------------

% get the shape of points so we can restore it later
[np,mp]=size(points);

% reshape into a vector
p=reshape(points-centre,np*mp,1);

%-------------------------------------------------
% compute the field
%-------------------------------------------------

% convert to polar coordinates
theta=angle(p);
rad=abs(p);

% make a matrix from n and rad
[nd,rd]=meshgrid(n,kwave*rad);

% get Bessel/Hankel/far field values as appropriate
if strcmp(type,'J')
    
    bess=besselj(abs(nd),rd);
    bessd=kwave*besseljd(abs(nd),rd);
    
elseif strcmp(type,'H')
    
    bess=besselh(abs(nd),rd);
    bessd=kwave*besselhd(abs(nd),rd);
    
end

% compute the angular part
Y=exp(1i*theta*n);
Yd=Y*diag(1i*n);

% put it together
M = bess.*Y;
Mdr = bessd.*Y;
Mdtheta = bess.*Yd;

%-------------------------------------------------
% make the return value the same shape as the original
% array of points
%-------------------------------------------------

% compute the sum of the wavefunctions and reshape
dr=reshape(Mdr*cof,np,mp);
dtheta=reshape(Mdtheta*cof,np,mp);

if nargout>2
    er = reshape(p./abs(p),np,mp);
    etheta = reshape(exp(1i*pi/2)*er,np,mp);
    rad = reshape(rad,np,mp);
end