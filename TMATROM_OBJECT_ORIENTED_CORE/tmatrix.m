% T-matrix class
%
%   T = tmatrix(fname) loads a T-matrix from file fname. The file is
%   assumed to be of Matlab *.mat format.
%
%   T = tmatrix(fname,type) loads a T-matrix from file fname. The file type
%   is specified by type; supported values are '-matlab' and '-tmatrom'.
%
%   T = tmatrix(n,k,M) creates a T-matrix object of order n, wavenumber k, 
%   and matrix M. The origin is set to be 0.
%
%   T = tmatrix(n,k,M,x) creates a T-matrix object of order n, wavenumber k, 
%   matrix M and origin x.
%
%   T = tmatrix(n,k,M,x,str) creates a T-matrix object of order n, wavenumber k, 
%   matrix M, origin x and comment str. 
%
% Also:
%
%   val = T.error() gives a measure of the error in the T-matrix based on a
%   symmetry relation. See Ganesh and Hawkins ANZIAM J. Vol 51 
%   Pages C215--C230 (2010) for details.
% 
%   T.setOrigin(x) sets the origin of the T-matrix to x. This is a virtual
%   origin that is only used in interactions with wave functions. This does
%   not change the T-matrix.
%
%   T.setComments(str) stores str as comments associated with the T-matrix.
%
%   T.getComments() prints any comments associated with the T-matrix.
%
%   val = T.getComments() returns any comments associated with the T-matrix
%   in val.
%
%   T.save(fname) saves the T-matrix in .mat format.
%
%   T.save(fname,'-matlab') saves the T-matrix in .mat format.
%
%   T.save(fname,'-tmatrom') saves the T-matrix in ASCII .tmat format.
%
% Example:
%
%   p = plane_wave(0,k);
%   u = regularwavefunctionexpansion(n,0,p);
%   T = tmatrix('sample.mat');
%   v = T * u;
%
%   Now v is a radiating wavefunction expansion for the scattered field
%   induced by the plane wave p. The T-matrix is loaded from file
%   sample.mat.
%
% See also: regularwavefunctionexpansion, radiatingwavefunctionexpansion,
% plane_wave, point_source, ghtmatrix.
%
% Stuart C. Hawkins - 6 March 2018

% Copyright 2014, 2015, 2016, 2017, 2018, 2022, 2023 Stuart C. Hawkins and M. Ganesh.
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


classdef tmatrix < handle
    
    properties
        
        % given properties
        order
        kwave
        origin
        matrix
        comments
    end
    
    methods

        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = tmatrix(varargin)
            
            if nargin==1
                
                % then we are loading from a file
                fname = varargin{1};
                
                % load the T-matrix
                self.load(fname);
            
            elseif nargin==2
                
                % then we are loading from a file
                fname = varargin{1};
                opts = varargin{2};
                
                % load the T-matrix
                self.load(fname,opts);
            
            else
                
                % then we are computing the T-matrix from a matrix
                self.order = varargin{1};
                self.kwave = varargin{2};
                self.matrix = varargin{3};
                
                % set default for origin
                if nargin < 4 || isempty(varargin{4})
                    self.origin = 0;
                else
                    self.origin = varargin{4};
                end
                
                % read comments if given
                if nargin < 5
                    self.comments = '';
                else
                    self.comments = varargin{5};
                end
            
            end
            
        end
        
        %-----------------------------------------
        % error check
        %-----------------------------------------

        % Check the symmetry relation (17) in Ganesh and Hawkins
        % ANZIAM J. 51 C215--C230 (2010)

        function val = error(self,opts)
           
            val = max(max(abs(self.matrix + self.matrix' ...
                + 2 * self.matrix' * self.matrix)));
            
            if nargin>1
                val = val / max(max(abs(self.matrix)));
            end
            
        end
                    
        %-----------------------------------------
        % multiply tmatrix x regular wave expansion
        %-----------------------------------------

        function val = mtimes(self,expansion)
            
            % - - - - - - - - - - - - - - - - - 
            % check the T-matrix and the wave 
            % expansion are compatible
            % - - - - - - - - - - - - - - - - - 

            if ~isa(expansion,'regularwavefunctionexpansion')
                
                error('expansion must be a regularwavefunctionexpansion')
                
            end
            
            if self.kwave ~= expansion.kwave
                
                error('T-matrix and expansion wavenumbers do not match.')
                
            end
            
            if self.origin ~= expansion.origin
                
                error('T-matrix and expansion centers do not match.')
                
            end
           
            if self.order ~= expansion.order
                
                error('T-matrix and expansion orders do not match.')
                
            end
            
            % - - - - - - - - - - - - - - - - - 
            % do product 
            % - - - - - - - - - - - - - - - - - 

            % create a radiating wave function expansion with coefficients
            % obtained by matrix multiplication with the T-matrix
            val = radiatingwavefunctionexpansion(self.order,self.origin,...
                self.kwave,self.matrix * expansion.coefficients(:));
            
        end
        
        %-----------------------------------------
        % set the origin
        %-----------------------------------------

        function setOrigin(self,origin)
            
            self.origin = origin;
            
        end
        
        %-----------------------------------------
        % save
        %-----------------------------------------

        function save(self,fname,opts)
            
            if nargin<3
                opts = '-matlab';
            end
            
            if strcmp(opts,'-matlab')
            
                % put the class properties in a struct
                data = struct('order',self.order,'kwave',self.kwave,...
                    'origin',self.origin,'matrix',self.matrix,...
                    'version',tmatrom_version(),'comments',self.comments);
                
                % save the struct
                save(fname,'-struct','data');
            
            elseif strcmp(opts,'-raw')
                
                save(fname,'-ascii','-double',self.matrix)
                
            elseif strcmp(opts,'-tmatrom')
                
                % open file
                fid = fopen(fname,'w');
                
                % write out data
                fprintf(fid,'%d\n',self.order);
                fprintf(fid,'%0.15e\n',self.kwave);
                fprintf(fid,'%0.15e\n',real(self.origin));
                fprintf(fid,'%0.15e\n',imag(self.origin));
                fprintf(fid,'%0.15e\n',real(self.matrix(:)));
                fprintf(fid,'%0.15e\n',imag(self.matrix(:)));
                fprintf(fid,'%0.15f\n',tmatrom_version);                
                fprintf(fid,'%s\n',self.comments);
                
                % close file
                fclose(fid);
                
            else
                
                error('Filetype %s not recognised',opts)
                
            end
                
        end        
        
        %-----------------------------------------
        % load
        %-----------------------------------------

        function load(self,fname,opts)
            
            if nargin<3
                opts = '-matlab';
            end
            
            if strcmp(opts,'-matlab')
            
                % put the class properties into a struct
                data = load(fname);
                
                % ** in future released we might need to check the version here
                % eg if class properties change **
                
                % set the class properties from the struct
                self.order = data.order;
                self.kwave = data.kwave;
                self.origin = data.origin;
                self.matrix = data.matrix;
                self.comments = data.comments;
            
            elseif strcmp(opts,'-tmatrom')
                
                % open file
                fid = fopen(fname,'r');
                
                % read data
                self.order = fscanf(fid,'%d',1);
                self.kwave = fscanf(fid,'%f',1);
                self.origin = fscanf(fid,'%f',1);
                temp_real = fscanf(fid,'%f',(2*self.order+1)^2);
                temp_imag = fscanf(fid,'%f',(2*self.order+1)^2);
                self.matrix = reshape(temp_real+1i*temp_imag,...
                    2*self.order+1,2*self.order+1);
                tmatrom_version = fscanf(fid,'%f',1);                
                
                % get comments... trim first and last character because
                % they will be white space
                temp = fscanf(fid,'%c');
                self.comments = temp(2:end);
                
                % close file
                fclose(fid);

            else
                
                error('Filetype %s not recognised',opts)
                
            end
            
        end
        
        %-----------------------------------------
        % set comments
        %-----------------------------------------

        function setComments(self,str)
            
            self.comments = str;
            
        end

        %-----------------------------------------
        % get comments
        %-----------------------------------------

        function varargout = getComments(self)
            
            if nargout == 0
                
                fprintf('%s\n',self.comments);
                
            else
                
                varargout{1} = self.comments;
                
            end
            
        end
        
    end % end methods
    
end

