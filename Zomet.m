%    SuperResolution: SuperResolution evaluation software
%    Copyright (C) 2011  S. Villena, M. Vega, D. Babacan, J. Mateos, 
%                        R. Molina and  A. K. Katsaggelos
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

% Implementation of RobustSR from Assaf Zomet, Alex Rav-Acha and Shmuel
% Peleg 
function [x,out,varargout] = Zomet(y,opt,HRaxis,handles)

fprintf(opt.LogFile, 'Running Zomet\n');
fprintf(opt.LogFile, '--------------------------------------------------------------------\n');

N = opt.N;
M = opt.M;

nopix = N*M;
[ys,yvec] = unwrapLR(y,opt);

H = circconvmatx2(opt.h,M,N);
A = dwnsmpl_matrix(M,N,opt);
W = sparse([]);

%% % Initialize HR image
if opt.InitImgExists == 1,
    x = opt.x0;
elseif opt.InitWithAvgImg,
    fprintf(opt.LogFile,'\nCalculate initial average image\n');
    for k=1:opt.L,
        dx = opt.sx(k);
        dy = opt.sy(k);
        theta = opt.theta(k);
        [C,Lbl,Lbr,Ltl,Ltr,a,b]  = warp_matrix_bilinear(dx,dy,theta,M,N);
        W = [W; A*H*C];
    end

    [avg_img] = get_avg_img(y,W,opt);
    x = avg_img;
    if opt.ShowImages, 
        axes(HRaxis), imshow(reshape(x,[opt.M,opt.N]));
    end
    
else
    [ys] = unwrapLR(y,opt);
    x = imresize(ys{1}, opt.res, 'bicubic');
    x=x(:);
end

maxPSNR = 0;
for i=1:opt.maxit,
    % set(handles.text_Number_iteration,'String',i); 
    oldx = x;
    
    % Compute the backprojection term
    T = zeros([nopix, opt.L]);
    for k=1:opt.L,
        [C,Lbl,Lbr,Ltl,Ltr,a,b]  = warp_matrix_bilinear(opt.sx(k),opt.sy(k),opt.theta(k),M,N);
        ykhat = A*H*C*x;
        
        G = ykhat-yvec{k};

        T(:,k) = (A*H*C)'*G;
    end
    
    % Take the median of T
    Med = median(T, 2);
    
    x = x - opt.L * opt.lambda * Med;
    
    xconv = norm(x - oldx)/norm(oldx);
    
    if ~opt.Real,
        MSEs(i) = sum(sum( (x(:) - opt.xtrue(:)).^2 ) )/nopix;
        PSNRs(i) = 10*log10(1/MSEs(i));
         set(handles.text_mse2,'String',num2str(MSEs(i)));
         set(handles.text_psnr2,'String',num2str(PSNRs(i),'%f dB'));

        if PSNRs(i) > maxPSNR,
            maxPSNR = PSNRs(i);
            maxPSNR_x = x;
            maxPSNR_it = i;
        end

        fprintf(opt.LogFile,'\n xconv = %g, MSE = %g, PSNR = %g',xconv, MSEs(i),PSNRs(i));
    end
    
    if opt.ShowImages, 
        set(handles.text_Msr,'String',M);
        set(handles.text_Nsr,'String',N);
        axes(HRaxis), imshow(reshape(x,[opt.M,opt.N]))
    end
    
    if opt.WriteImages & mod(i,10) == 0,
        %save(sprintf('x_var%d_it%d.mat',strcmp(opt.method,'variational'), i), 'x','MSEs', 'opt');
        imwrite(reshape(x,[opt.M,opt.N]),sprintf('x_RegErr_Zomet_sigma%g_it%d.png',opt.sigma, i));
    end
    
    if opt.KeepHistory,
        history_x{i} = x;
    end
    
    if xconv < opt.thr,
        break;
    end
    if ~opt.Real & PSNRs(i) < 0,
        break;
    end
    % opt.cancel = get(opt.ref_cancel,'Value');
    % if opt.cancel,
    %     break;
    % end
end

if ~opt.Real,
    out.MSEs = MSEs;
    out.MSE = MSEs(end);
    out.PSNRs = PSNRs;
    out.PSNR = PSNRs(end);
    out.maxPSNR  = maxPSNR ;
    out.maxPSNR_x = reshape(maxPSNR_x,M,N);
    out.maxPSNR_it = maxPSNR_it;
end
out.it = i;    
  
if opt.KeepHistory,
    varargout{1} = history_x;
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    