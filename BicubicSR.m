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

function  [x,out] = BicubicSR(y,opt)


N = opt.N;
M = opt.M;
nopix = opt.N*opt.M;

[ys] = unwrapLR(y,opt);
x = imresize(ys{1}, opt.res, 'bicubic');
x=x(:);

if ~opt.Real,
    MSE = sum(sum( (x(:) - opt.xtrue(:)).^2 ) )/nopix;
    PSNR = 10*log10(1/MSE);

    out.MSE = MSE;
    out.PSNR = PSNR;
else
    out = [];
end
