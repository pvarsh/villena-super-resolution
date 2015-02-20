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

function H = circconvmatx2(h, M,N)

[mh, nh] = size(h);
if mh~=nh,
    error('blur kernel must be square');
end

blockheight = N;
blockwidth = M;

for i=((nh+1)/2):nh,
    
    %Take the kernel
    h0 = h(:,i);
    h0 = h0(:)';
    h0 = fliplr(h0);
    
    h0 = circshift( padarray( h0, [0 N-nh], 'post' ), [0 -(nh+1)/2+1] );
    
    if sum(find(abs(h0)>0)) == 0,
        H0 = sparse(M,M);
    else
        H0 = sparse(circulant(h0'));
    end
    if i==(nh+1)/2,
        H = H0;
    else
        H = [H H0];
    end

end

for i=(nh+1)/2-1:-1:1,
    
    %Take the kernel
    h0 = h(:,i);
    h0 = h0(:)';
    h0 = fliplr(h0);
    
    h0 = circshift( padarray( h0, [0 N-nh], 'post' ), [0 -(nh+1)/2+1] );
    
    if sum(find(abs(h0)>0)) == 0,
        H0 = sparse(M,M);
    else
        H0 = sparse(circulant(h0'));
    end
    
    H = [H0 H];
   
end


H = [H, sparse(zeros(N,N*(N-nh)))];

H = circshift( H, [0 -N*((nh+1)/2-1)] );
H0 = H;
for i=1:N-1,
    H = [H; circshift( H0, [0 N*i] )];
end

   
    