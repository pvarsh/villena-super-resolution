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

function A = dwnsmpl_matrix(M,N,opt)
%Creates downsampling matrix of size M/d*N/d by M/d*N/d

% Populate the locations to sample

nopixels = M*N;

dindices = reshape(1:nopixels,M,N);
dindices = dindices(1:opt.res:M,1:opt.res:N);
dindices = dindices(:);

if length(dindices) ~= opt.m*opt.n,
    error('Function dwnsmpl_matrix: matrix-vector sizes don''t match');
end

A = sparse([1:opt.m*opt.n], [dindices'],[ones(size(dindices))], opt.m*opt.n, nopixels);


