%    SuperResolution: SuperResolution evaluation software
%    Copyright (C) 2011  S. Villena, M. Vega, D. Babacan, J. Mateos, 
%                        R. Molina and  A. K. Katsaggelos
%
%    If you use this software to evaluate any of the methods, please cite 
%    the corresponding papers (see manual).
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function C = shift_matrix(dx,dy)

% Get image dimensions
[M N dims] = size(dx);
nopixels = M*N;

dindices = [1:nopixels]' + dx(:)*M + dy(:) ;
% Derin: These cause repetitions
dindices(find(dindices>nopixels)) = nopixels; 
dindices(find(dindices<1)) = 1;

C = sparse([1:nopixels], [dindices'],[ones(size(dindices))], nopixels, nopixels);


