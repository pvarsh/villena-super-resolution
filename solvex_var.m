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


%% TODO: show images as iterations update
function [x,out,varargout] = solvex_var(y,opt,HRaxis,handles)

if strcmp(opt.method,'variational'),
    fprintf(opt.LogFile, 'Running SR with optimal distributions\n');
else
    fprintf(opt.LogFile, 'Running SR with degenerate distributions\n');
end
fprintf(opt.LogFile, '--------------------------------------------------------------------\n');

nopix = opt.N*opt.M;
nu = 1e-4;
N = opt.N;
M = opt.M;

dx = [0 0 0;
     -1 1 0;
      0 0 0];
dy = [0 -1 0;
      0 1 0;
      0 0 0];

Dx = circconvmatx2(dx, opt.M, opt.N);
Dy = circconvmatx2(dy, opt.M, opt.N);

H = circconvmatx2(opt.h,opt.M,opt.N);
A = dwnsmpl_matrix(opt.M,opt.N,opt);
W = sparse([]);

[X Y] = meshgrid([-N/2:N/2-1],[-N/2:N/2-1]);
X = X(:);
Y = Y(:);
  
if opt.InitImgExists == 1,
	opt.sx = opt.sx0;
	opt.sy = opt.sy0;
	opt.theta = opt.theta0;
else
	opt.sx = opt.sx_init;
	opt.sy = opt.sy_init;
	opt.theta = opt.theta_init;
end



%% If the registration is not accurate, estimate the registration first 
if ~opt.Real,
    if opt.EstimateRegistration,
        [ys,yvecs] = unwrapLR(y,opt);
        x = imresize(ys{1}, opt.res, 'bicubic'); % Start with a bicubic interpolated image to get a reference frame
        x=x(:);

        fprintf(opt.LogFile,'Estimating the registration parameters\n');

        %for k=1:opt.L,
        for k=2:opt.L,    
            fprintf(opt.LogFile,'\n********************************************\n');
            fprintf(opt.LogFile,'Image %d\n',k);
            [newsk, Lambdak,yhatk] = LKvar(x, k, yvecs{k}, 0 , A,H, zeros(3), 1, opt);

            if ~opt.Real,
                sk_true = [opt.theta_true(k); opt.sx_true(k); opt.sy_true(k)];
                sk_true = sk_true(:);
            end

            fprintf(opt.LogFile,'Initial registration   sx = %g, sy = %g, theta = %g\n',opt.sx_init(k),opt.sy_init(k),opt.theta_init(k));
            fprintf(opt.LogFile,'Previous registration  sx = %g, sy = %g, theta = %g\n',opt.sx(k),opt.sy(k),opt.theta(k));
            fprintf(opt.LogFile,'Estimated registration sx = %g, sy = %g, theta = %g\n',newsk(2),newsk(3),newsk(1));
            if ~opt.Real,
                fprintf(opt.LogFile,'True registration      sx = %g, sy = %g, theta = %g\n',opt.sx_true(k),opt.sy_true(k),opt.theta_true(k));
                fprintf(opt.LogFile,'Estimation error NMSE = %g\n', norm(newsk-sk_true)^2);
            end

            opt.theta(k) = newsk(1);
            opt.sx(k) = newsk(2)-.5;
            opt.sy(k) = newsk(3)-.5;
            Lambdas{k} = Lambdak;
            %             pause
        end
    end
end
for k=1:opt.L,
    dx = opt.sx(k);
    dy = opt.sy(k);
    theta = opt.theta(k);
    [C,Lbl,Lbr,Ltl,Ltr,a,b]  = warp_matrix_bilinear(dx,dy,theta,M,N);
    W = [W; A*H*C];
end

%% % Initialize
if opt.InitImgExists == 1,
    x = opt.x0;
elseif opt.InitWithAvgImg,
    fprintf(opt.LogFile,'\nCalculating initial average image\n');
    [avg_img] = get_avg_img(y,W,opt);
    x = avg_img;
else
    [ys] = unwrapLR(y,opt);
    x = imresize(ys{1}, opt.res, 'bicubic');
    x=x(:);
end

% if opt.interv_ShowImages , 
%     axes(HRaxis), imshow(reshape(x,[opt.M,opt.N]));
    
% end


%% Calculate initial hyperparameters
[ys,yvecs] = unwrapLR(y,opt);

u = (Dx*x(:)).^2 + (Dy*x(:)).^2;
u = u+eps;
T = spdiags(1./sqrt(u), 0, nopix, nopix);
alpha_im = (nopix/2 )/sum(sum(sqrt(u)));

e = y-W*x; 

if strcmp(opt.noise_estimate, 'SEPARATE')
    [temp, es] = unwrapLR(e,opt);
    for k=1:opt.L,
        betak(k) = length(es{k}) / sum(es{k}.^2);
    end
     if opt.FixedBeta
                      betak(1:end)=1/opt.sigma;
     end
else
    e = sum(e.^2);
    betak = length(y) / e * ones(opt.L,1);
     if opt.FixedBeta
                      betak=1/opt.sigma;
     end
end


if  opt.FixedParameters,
    if opt.CalcParamFromTrue,
        
        u_true = (Dx*opt.xtrue(:)).^2 + (Dy*opt.xtrue(:)).^2;
        %T = spdiags(1./sqrt(u_true), 0, nopix, nopix);
        alpha_im = (nopix/2 )/sum(sum(sqrt(u_true)));
        
        e = y-W*opt.xtrue(:);
        
        if strcmp(opt.noise_estimate, 'SEPARATE')
            [temp, es] = unwrapLR(e,opt);
            for k=1:opt.L,
                betak(k) = length(es{k}) / sum(es{k}.^2);
            end
        else
            e = sum(e.^2);
            betak = length(y) / e * ones(opt.L,1);
        end
    else
        betak = opt.fixed_betak;
        alpha_im = opt.fixed_alpha_im;
    end
end

if ~opt.Real,
    MSE = sum(sum( (x(:) - opt.xtrue(:)).^2 ) )/nopix;
    PSNR = 10*log10(1/MSE);
    fprintf(opt.LogFile,'betak = \n');
    fprintf(opt.LogFile,'%g ',betak);
    fprintf(opt.LogFile,'\nalpha_im = %g, MSE = %g, PSNR = %g\n',alpha_im,MSE,PSNR);
else
    fprintf(opt.LogFile,'betak = \n');
    fprintf(opt.LogFile,'%g ',betak);
    fprintf(opt.LogFile,'\nalpha_im = %g\n',alpha_im);
end
%% Initial registration covariance matrices
for k=1:opt.L,
             
    Lambdas_p{k} = [1e-6  0.0   0.0;
                    0.0   1e-6  0.0;
                    0.0   0.0   1e-6]; 

	Lambdas_p{k} = [1e-6 1e-8 1e-8;
                  1e-8 1e-6  1e-8;
                  1e-8 1e-8 1e-6];

    Lambdas_p{k} = zeros(3);              
    
    
    Lambdas{k} = Lambdas_p{k};
    
end
              
xconv = 1;
maxPSNR = 0;
maxPSNR_x = x;
maxPSNR_it = 1;
maxPSNR_sx = opt.sx;
maxPSNR_sy = opt.sy;
maxPSNR_theta = opt.theta;


for i=1:opt.maxit,
    disp(['>>>> SR iteration ', num2str(i), ' (maxit = ', num2str(opt.maxit), ')...']);
    % opt.cancel = get(opt.ref_cancel,'Value');
    % if opt.cancel,
    %    break;
    % end
    fprintf(opt.LogFile,'\n-------------------------------------------------------------------\n');
    fprintf(opt.LogFile,'Iteration %d\n',i);
   % set(handles.edit_Niteration,'String',num2str(i));
    % set(handles.text_Number_iteration,'String',i); 

    oldx = x(:);

    %% Estimate the image
    % PCG
    
    Sigma_inv = sparse(nopix,nopix);
    W = sparse([]);
    rhs = 0;

    for k=1:opt.L,
        
        [C,Lbl,Lbr,Ltl,Ltr,a,b]  = warp_matrix_bilinear(opt.sx(k),opt.sy(k),opt.theta(k),M,N);
        W = [W; A*H*C];
        
        rhs = rhs + betak(k)*(A*H*C)'*yvecs{k};
        Sigma_inv = Sigma_inv + betak(k)*(A*H*C)'*A*H*C;
        
        % Create the terms related to the uncertainty arising from registration
        if strcmp(opt.method,'variational'),
            P1 = spdiags( -X*sin(opt.theta(k)) - Y*cos(opt.theta(k)), 0 , nopix,nopix);
            P2 = spdiags( +X*cos(opt.theta(k)) - Y*sin(opt.theta(k)), 0 , nopix,nopix);

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
            clear O;

            Lambda = Lambdas{k};
            Sigma_inv = Sigma_inv + betak(k) * Lambda(1,1)*O11 + betak(k) *Lambda(2,2)*O22 + betak(k) *Lambda(3,3)*O33 ...
                    + betak(k) *2*Lambda(1,2)*O12 + betak(k) *2*Lambda(1,3)*O13 + betak(k) *2*Lambda(2,3)*O23;
            
        end 
    end
    if strcmp(opt.method,'variational'),
        fprintf(opt.LogFile,'Finished calculating Sigma registration terms\n');
    end
    
    Sigma_inv = Sigma_inv + alpha_im * Dx'*T*Dx + alpha_im * Dy'*T*Dy;

    fprintf(opt.LogFile,'Running PCG to estimate the image\n');
    [x,flag,relres,iter] = pcgmod(Sigma_inv,rhs,opt.PCG_thr,opt.PCG_maxit,[],[],x, opt.PCG_minit);

    fprintf(opt.LogFile,'PCG returned the solution in %d iterations, flag = %d\n',iter, flag);
    
    if opt.ApproxSigma,
        Sigma = 1./spdiags(Sigma_inv,0);
        Sigma = spdiags(Sigma, 0, nopix,nopix);
    else
        Sigma = inv(Sigma_inv);
    end
    
    %% Estimate the registration parameters
    
    if opt.EstimateRegistration, 
        fprintf(opt.LogFile,'Estimating the registration parameters\n');
        
        %for k=1:opt.L,
        for k=2:opt.L,
            fprintf(opt.LogFile,'\n********************************************\n');
            fprintf(opt.LogFile,'Image %d\n',k);
            [newsk, Lambdak,yhatk] = LKvar(x, k, yvecs{k}, Sigma, A,H, Lambdas_p{k}, betak(k), opt);
             
            if ~opt.Real, 
                sk_true = [opt.theta_true(k); opt.sx_true(k); opt.sy_true(k)];
                sk_true = sk_true(:);                        
            end
    
            fprintf(opt.LogFile,'Initial registration   sx = %g, sy = %g, theta = %g\n',opt.sx_init(k),opt.sy_init(k),opt.theta_init(k)/pi*180);
            fprintf(opt.LogFile,'Previous registration  sx = %g, sy = %g, theta = %g\n',opt.sx(k),opt.sy(k),opt.theta(k)/pi*180);
            fprintf(opt.LogFile,'Estimated registration sx = %g, sy = %g, theta = %g\n',newsk(2),newsk(3),newsk(1)/pi*180);
            if ~opt.Real, 
                fprintf(opt.LogFile,'True registration      sx = %g, sy = %g, theta = %g\n',opt.sx_true(k)/pi*180,opt.sy_true(k)/pi*180,opt.theta_true(k)/pi*180);
                fprintf(opt.LogFile,'Estimation error NMSE = %g\n', norm(newsk-sk_true)^2);
            end

 	        opt.theta(k) = newsk(1);
            opt.sx(k) = newsk(2);
            opt.sy(k) = newsk(3);
            if ~isempty(Lambdak),
                Lambdas{k} = Lambdak;
            end
%             pause

            
        end
    end
    
    %% Estimate the hyperparameters
    if  opt.FixedParameters == 0,
        fprintf(opt.LogFile,'Estimating the hyperparameters\n');
        
        if strcmp(opt.method,'variational'),
            Sigma_mat = reshape(spdiags(Sigma, 0),opt.N,opt.M);
            DxSigma = Sigma_mat + circshift(Sigma_mat,[0 -1]);
            DySigma = Sigma_mat + circshift(Sigma_mat,[-1 0]);
            if opt.DIVIDE_U,
                u = (Dx*x).^2 + (Dy*x).^2 + 1/nopix*DxSigma(:) + 1/nopix*DySigma(:);
                fprintf(opt.LogFile,'avg DxSigma = %g, avg DySigma = %g\n',1/nopix*mean(DxSigma(:)),1/nopix*mean(DySigma(:)));
            else
                u = (Dx*x).^2 + (Dy*x).^2 + DxSigma(:) + DySigma(:);
                fprintf(opt.LogFile,'avg DxSigma = %g, avg DySigma = %g\n',mean(DxSigma(:)),mean(DySigma(:)));
            end
            
            
        else
            u = (Dx*x).^2 + (Dy*x).^2;
        end
        
        T = spdiags(1./sqrt(u), 0, nopix, nopix);
        alpha_im = (nopix/2)/sum(sum(sqrt(u)));
        
        e = y-W*x;
        
        % Different noise variances??
        if strcmp(opt.noise_estimate, 'SEPARATE')
            [temp, es] = unwrapLR(e,opt);
            
            for k=1:opt.L,
                esk(k) = sum(es{k}.^2);
                
                if strcmp(opt.method,'variational'),
                    
                    [C,Lbl,Lbr,Ltl,Ltr,a,b]  = warp_matrix_bilinear(opt.sx(k),opt.sy(k),opt.theta(k),M,N);
                    
                    traceBkSigma = sum( spdiags( (A*H*C)'*(A*H*C)* Sigma, 0 ) );
                    
                    P1 = spdiags( -X*sin(opt.theta(k)) - Y*cos(opt.theta(k)), 0 , nopix,nopix);
                    P2 = spdiags( +X*cos(opt.theta(k)) - Y*sin(opt.theta(k)), 0 , nopix,nopix);
                    
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
                    clear O;
                    
                    Lambda = Lambdas{k};
                    xOx = Lambda(1,1)*x'*O11*x + Lambda(2,2)*x'*O22*x + Lambda(3,3)*x'*O33*x ...
                        + 2*Lambda(1,2)*x'*O12*x + 2*Lambda(1,3)*x'*O13*x + 2*Lambda(2,3)*x'*O23*x;
                    traceOSigma = sum ( spdiags( Lambda(1,1)*O11*Sigma + Lambda(2,2)*O22*Sigma + Lambda(3,3)*O33*Sigma ...
                        + 2*Lambda(1,2)*O12*Sigma + 2*Lambda(1,3)*O13*Sigma + 2*Lambda(2,3)*O23*Sigma, 0) );
                    if opt.verbose,
                        fprintf(opt.LogFile,'k = %d, e = %g, xOx = %g, traceBkSigma = %g, traceOSigma = %g\n', k , esk(k), xOx, traceBkSigma, traceOSigma);
                    end
                    esk(k) = esk(k) + traceBkSigma + traceOSigma + xOx;
                    betak(k) = (length(es{k})) / esk(k);
                else
                    if opt.verbose,
                        fprintf(opt.LogFile,'k = %d, e = %g\n', k , esk(k));
                    end
                    betak(k) = (length(es{k})) / esk(k);
                end
                 if opt.FixedBeta
                      betak(1:end)=1/opt.sigma;
                 end
            end
        else
            e = sum(e.^2);
            betak = (length(y)) / e * ones(opt.L,1);
             if opt.FixedBeta
                      betak=1/opt.sigma;
             end
        end
    end
       
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
            maxPSNR_sx = opt.sx;
            maxPSNR_sy = opt.sy;
            maxPSNR_theta = opt.theta;
        end

        svec = [opt.theta(:); opt.sx(:); opt.sy(:)];
        svec_true = [opt.theta_true(:); opt.sx_true(:); opt.sy_true(:)];
        %sNMSEs(i) = (norm(svec-svec_true)^2)/(norm(svec_true)^2);
        sNMSEs(i) = (norm(svec-svec_true)^2);
    end
    
    
    if opt.WriteImages & mod(i,3) == 0,
        disp(['Writing image' num2str(i)]);
        xout = uint8(reshape(x, [opt.M, opt.N]));
        % imwrite(xout,sprintf('x_RegErr_var%d_sigma%g_init%d_it%d.png',strcmp(opt.method,'variational'), opt.sigma, opt.DIVIDE_U, i));
        imwrite(xout, [opt.outpath '/' opt.out_file_prefix '_iter_' num2str(i) '.png']);
    end
            
    if ~opt.Real, 
        fprintf(opt.LogFile,'betak = \n');
        fprintf(opt.LogFile,'%g ',betak);
        fprintf(opt.LogFile,'\nalpha_im = %g,xconv = %g, PSNR = %g, MSE = %g, sNMSE = %g\n',alpha_im,xconv, PSNRs(i), MSEs(i), sNMSEs(i));
    else
        fprintf(opt.LogFile,'betak = \n');
        fprintf(opt.LogFile,'%g ',betak);
        fprintf(opt.LogFile,'\nalpha_im = %g,xconv = %g\n',alpha_im,xconv);
    end
%     pause
    
    if opt.ShowImages,
        HRaxis = axes();
        % set(handles.text_Msr,'String',M);
        % set(handles.text_Nsr,'String',N);
        axes(HRaxis), imshow(reshape(x,[opt.M,opt.N]));
    end
    
    if opt.KeepHistory,
        history_x{i} = x;
        history_alpha{i} = alpha_im;
        history_betak{i} = betak;
        history_theta{i} = opt.theta;
        history_sx{i} = opt.sx;
        history_sy{i} = opt.sy;
        history_Lambda{i} = Lambdas;
    end
    
    if xconv < opt.thr,
        break;
    end
    % opt.cancel = get(opt.ref_cancel,'Value');
    % if opt.cancel,
    %     break;
    % end
end

out.betak = betak;
out.alpha_im = alpha_im;
out.Lambdas = Lambdas;
out.xconv = xconv;
out.theta = opt.theta;
out.sx = opt.sx;
out.sy = opt.sy;
out.it = i;

if ~opt.Real, 
    out.MSEs = MSEs;
    out.MSE = MSEs(end);
    out.PSNRs = PSNRs;
    out.PSNR = PSNRs(end);
    out.sNMSEs = sNMSEs;
    out.sNMSE = sNMSEs(end);
    out.maxPSNR  = maxPSNR ;
    out.maxPSNR_x = reshape(maxPSNR_x,M,N);
    out.maxPSNR_it = maxPSNR_it;
    out.maxPSNR_sx = maxPSNR_sx;
    out.maxPSNR_sy = maxPSNR_sy;
    out.maxPSNR_theta = maxPSNR_theta;
end

if opt.KeepHistory,
    varargout{1} = history_x;
    varargout{2} = history_alpha;
    varargout{3} = history_betak;
    varargout{4} = history_theta;
    varargout{5} = history_sx;
    varargout{6} = history_sy;
    varargout{7} = history_Lambda;
end

function e = local_func(x,W,y,T,alpha_im,betak,Dx,Dy)
e = y-W*x';
e = sum(e.^2);
Dxx = Dx*x';
Dyy = Dy*x';
e = betak*e + alpha_im*Dxx'*T*Dxx + alpha_im*Dyy'*T*Dyy;
return

function g = local_grad(x,W,y,T,alpha_im,betak,Dx,Dy)
r = y-W*x';
g = W*r;
Dxx = Dx*x';
Dyy = Dy*x';
g = -2.*g + 2*alpha_im*Dx'*T*Dxx + 2*alpha_im*Dy'*T*Dyy;
g = g';
return
