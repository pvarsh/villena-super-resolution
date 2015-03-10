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

%----------
% Parameters: (defaults used in this package in [])
%               x:  first image (all registration done to first image)
%               k:  frame number in loaded sequence
%               yk: kth image
%               Sigma: don't know [0]
%               A: don't know [1]
%               H: don't know [1]
%               Lambdapk: 3x3 matrix [zeros(3)]
%               betak: don't know [1]
%               opt: options object
%                   - opt.N: super-resolved dimension 1
%                   - opt.M: super-resolved dimenison 2
%                   - theta
%                   - opt.sx
%                   - opt.sy
%                   - opt.theta_init
%                   - opt.sx_init
%                   - opt.sy_init
%                   - opt.LK_maxit: max iterations
%                   - opt.debug: 
%                   - opt.LK_thr: threshold?
%                   - opt.method: [variational]
% Calls:
%               warp_matrix_bilinear()
% Notes: on first run opt.theta_init, opt.theta, opt.sx_init, opt_sx, opt.sy_init, opt.sy
%        are all [0,0].
% Returns: newsk = [theta, sx, sy]]
function [newsk, Lambdak,varargout] = LKvar(x, k, yk, Sigma, A,H, Lambdapk, betak, opt)
                 

nopix = opt.N*opt.M;
N = opt.N;
M = opt.M;
[X Y] = meshgrid([-N/2:N/2-1],[-N/2:N/2-1]);
X = X(:);
Y = Y(:);

newsk = [opt.theta(k), opt.sx(k), opt.sy(k)]';
sk = newsk;
debug = 1;

if sum(sum(abs(Lambdapk))) == 0,
    Lambdapk_inv = zeros(3);
else
    Lambdapk_inv = inv(Lambdapk);
end
ccc_Lambdak_min = Lambdapk;

s_init = [opt.theta_init(k); opt.sx_init(k); opt.sy_init(k)];

ccc_e_min = 9.9e10;
sk_min = sk; 

for i=1:opt.LK_maxit,
    
    thetak = sk(1);
    dx = sk(2);
    dy = sk(3);
    
    [C,Lbl,Lbr,Ltl,Ltr,a,b]  = warp_matrix_bilinear(dx,dy,thetak,M,N);

%     disp(A);
%     disp(H);
%     disp(C);
%     disp(size(A));
%     disp(size(H));
%     disp(size(C));
%     disp(size(x));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ccc_varargout = A*H*C*x;
     if i==1
         ccc_varargout_min = ccc_varargout;
         Lambdak = ccc_Lambdak_min;
     end
     ccc_e = norm(yk-ccc_varargout,2);
    if ccc_e <= ccc_e_min
        ccc_e_min = ccc_e;
        sk_min = newsk;
        ccc_varargout_min = ccc_varargout;
        ccc_Lambdak_min = Lambdak ;
        i_max = i;
        if opt.debug,
            fprintf('%d:',i_max);
        end
        
   else
        newsk = sk_min;
        varargout{1} = ccc_varargout_min ;%A*H*C*x;
        Lambdak = ccc_Lambdak_min ; 
       break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    
    
    P1 = spdiags( -X*sin(thetak) - Y*cos(thetak), 0 , nopix,nopix);
    P2 = spdiags( +X*cos(thetak) - Y*sin(thetak), 0 , nopix,nopix);
    
    O{2} = spdiags((1-b), 0, nopix,nopix)*(Ltr-Ltl) + spdiags(b, 0, nopix,nopix)*(Lbr-Lbl); %M1
    O{3} = spdiags((1-a), 0, nopix,nopix)*(Lbl-Ltl) + spdiags(a, 0, nopix,nopix)*(Lbr-Ltr); %M2
        
    O{1} = P1*O{2} + P2*O{3};
    
    O{1} = A*H*O{1};
    O{2} = A*H*O{2};
    O{3} = A*H*O{3};
    
    O11 = O{1}'*O{1};
    O22 = O{2}'*O{2};
    O33 = O{3}'*O{3};
    O12 = O{1}'*O{2};
    O13 = O{1}'*O{3};
    O23 = O{2}'*O{3};

    O1x = O{1}*x;
    O2x = O{2}*x;
    O3x = O{3}*x;
         
    if opt.debug,
        err_img = (yk-A*H*C*x);
%         figure(1), imshow(reshape(O1x,[opt.n,opt.m])),title('O1x');
%         figure(2), imshow(reshape(O2x,[opt.n,opt.m])),title('O2x');
%         figure(3), imshow(reshape(O3x,[opt.n,opt.m])),title('O3x');
        figure(14), imshow(reshape(yk,[opt.n,opt.m])),title('Original');
        figure(15), imshow(reshape(A*H*C*x,[opt.n,opt.m])),title('Estimated');
        figure(16), imagesc(reshape(err_img,[opt.n,opt.m])), title('Error');
    end
    
    Phik = [ x'*O11*x, x'*O12*x, x'*O13*x;
             x'*O12*x, x'*O22*x, x'*O23*x;
             x'*O13*x, x'*O23*x, x'*O33*x];            

    
    Lambdak_inv = Lambdapk_inv + betak*Phik;
    
    rhs = Lambdapk_inv*s_init + betak*Phik*sk + betak*[(yk-A*H*C*x)'*O1x;(yk-A*H*C*x)'*O2x;(yk-A*H*C*x)'*O3x];
         
    if strcmp(opt.method,'variational'),
    
        Psik = [trace(O11*Sigma), trace(O12*Sigma), trace(O13*Sigma);
                trace(O12*Sigma), trace(O22*Sigma), trace(O23*Sigma);
                trace(O13*Sigma), trace(O23*Sigma), trace(O33*Sigma)];
        
        Lambdak_inv = Lambdak_inv + Psik;
        
        rhs = rhs + betak*Psik*sk - betak*[trace( (A*H*C)'*O{1}*Sigma ); trace( (A*H*C)'*O{2}*Sigma ); trace( (A*H*C)'*O{3}*Sigma )];
        
    end
    
    if rank(Lambdak_inv) < 3,
        if i==1,
            Lambdak = [];
        end
        break;
    end
        
    Lambdak = inv(Lambdak_inv);
    newsk = Lambdak * rhs;
    
   

    
    if norm(newsk-sk)/norm(sk) < opt.LK_thr,
        sk = newsk;
        break;
    end
    
%     fprintf(2,'Estimated registration parameters at iteration i %d, conv = %g\n',i, norm(newsk-sk)/norm(sk));
%     fprintf(2,'sx = %g, sy = %g, theta = %g\n',newsk(2),newsk(3),newsk(1));
    
    sk = newsk;      
    clear O;
%     opt.cancel = get(opt.ref_cancel,'Value');
%      if opt.cancel,
%         break;
%     end
end

newsk = sk_min;
if opt.debug,
    fprintf(2,'LKvar converged in %d iterations with err = %g at iteration = %d\n',i, ccc_e_min, i_max);
end
varargout{1} = ccc_varargout_min ;%A*H*C*x;
Lambdak = ccc_Lambdak_min ; 


