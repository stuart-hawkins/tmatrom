% Interlace vectors
%
%   a = interlace(a1,a2,...,an) returns the vector a of length n*m created 
%   by interlacing the entries of the vectors a1,...,an having length m.
%
% Example:
%
%   a = interlace([1 2 3],[4 5 6],[7 8 9]) returns a = [1 4 7 2 5 8 3 6 9];
%
% Stuart C. Hawkins - 23 November 2017

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


function val = interlace(varargin)

% get the number if input vectors
n = nargin;

% get the length of the first input vector
m = length(varargin{1});

% check that the lengths of the other input vectors are the same as the
% lengths of the first
for k=2:nargin
    if length(varargin{k}) ~= m
        error('input arrays must have same length')
    end
end

% initialise return array
val = zeros(1,n*m);

% insert the entries from the inputs into the return array
for k=1:n
    val(k:n:end) = varargin{k};
end

