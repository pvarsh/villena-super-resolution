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

MFSR = get(handles.edit_MFSR,'String');
MFSR = round(abs(str2num(MFSR)));
if isempty(MFSR)
    MFSR = 2;
end
set(handles.edit_MFSR,'String',num2str(MFSR));
handles.opt.P = MFSR;

handles.output = hObject;
guidata(hObject, handles);
handles.opt.sx = handles.opt.sx_init;
handles.opt.sy = handles.opt.sy_init;
handles.opt.theta = handles.opt.theta_init;
handles.opt.P = double(handles.opt.res);
handles.opt.alpha = 0.6;
handles.opt.beta = 1/100;
handles.opt.lambda = 0.005;
handles.opt.maxit = 40;
handles.opt.InitImgExists = 0;
handles.opt.ShowImages = 1;
handles.opt.WriteImages = 0;
handles.opt.KeepHistory = 0;
handles.opt.thr = 1e-5;
handles.opt.PCG_maxit = 100;
handles.opt.LK_maxit = 100;
handles.opt.debug = 0;
handles.opt.LK_thr =1e-4;
handles.opt.noise_estimate = 'SEPARATE';
handles.opt.PCG_thr =1e-10;
handles.opt.PCG_minit = 10;
handles.opt.verbose = 1;
handles.opt.FixedBeta = 0 ; 
handles.opt.lambda_prior = 1;
handles.opt.method = 'variational';

handles.opt.ApproxSigma = 1;

handles.opt.DIVIDE_U = 0;
handles.opt.FixedParameters = 0;