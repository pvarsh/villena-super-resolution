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

function [y,W] = create_data(opt)

M = opt.M;
N = opt.N;

H = circconvmatx2(opt.h,M,N);
A = dwnsmpl_matrix(M,N,opt);
W = sparse([]);

for i=1:opt.L,
    dx = opt.sx_true(i);
    dy = opt.sy_true(i);
    theta = opt.theta_true(i);
    %C = warp_matrix(dx,dy,theta,M,N);
    C = warp_matrix_bilinear(dx,dy,theta,M,N);
    
    W = [W; A*H*C];
end

y = W*opt.xtrue(:) + opt.sigma/255*randn(opt.L* opt.m*opt.n,1); 



   
