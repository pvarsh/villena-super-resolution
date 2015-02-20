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

function [C,varargout] = warp_matrix_bilinear(sx,sy,theta,M,N)

if M ~= N,
    error('image must be square');
end

% clear
% clc
% N = 6;
% M = 6;
% 
% theta = pi/2;
% sx = 0;
% sy = 0;

% X: horizontal
% Y: vertical

[X Y] = meshgrid([-N/2:N/2-1],[-N/2:N/2-1]);

indices = [X(:)'; Y(:)'; ones(1,N*M)];

S = [cos(theta) -sin(theta) sx;
     sin(theta) cos(theta)  sy;];
     
new_indices = S*indices;

dx = new_indices(1,:)-indices(1,:);
dy = new_indices(2,:)-indices(2,:);

dx = dx(:);
dy = dy(:);

a = dx - floor(dx);
b = dy - floor(dy);

Da = spdiags(a,0,N*M,N*M);
Db = spdiags(b,0,N*M,N*M);

dx = reshape(dx,[M,N]);
dy = reshape(dy,[M,N]);

% CHECK: To prevent errors coming from integer shifts
dx = dx + 1e-6;
dy = dy + 1e-6;

Lbl = shift_matrix(floor(dx),ceil(dy));
Lbr = shift_matrix(ceil(dx),ceil(dy));
Ltl = shift_matrix(floor(dx),floor(dy));
Ltr = shift_matrix(ceil(dx),floor(dy));

if sum(abs(a)) == 0 && sum(abs(b)) == 0,
    C = shift_matrix(floor(dx),floor(dy));
else
    %C = Db*(1-Da)*Lbl + Da*Db*Lbr + (1-Da)*(1-Db)*Ltl + (1-Db)*Da*Ltr;
    C = spdiags(b.*(1-a), 0, N*M,N*M)*Lbl;
    C = C+ spdiags(b.*a, 0, N*M,N*M)*Lbr;
    C = C+ spdiags((1-b).*(1-a), 0, N*M,N*M)*Ltl;
    C = C+ spdiags((1-b).*a, 0, N*M,N*M)*Ltr;

end

varargout{1} = Lbl;
varargout{2} = Lbr;
varargout{3} = Ltl;
varargout{4} = Ltr;
varargout{5} = a;
varargout{6} = b;

