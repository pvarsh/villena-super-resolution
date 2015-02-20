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

function varargout = MatrixReader(varargin)
% MATRIXREADER M-file for MatrixReader.fig
%      MATRIXREADER, by itself, creates a new MATRIXREADER or raises the existing
%      singleton*.
%
%      H = MATRIXREADER returns the handle to a new MATRIXREADER or the handle to
%      the existing singleton*.
%
%      MATRIXREADER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATRIXREADER.M with the given input arguments.
%
%      MATRIXREADER('Property','Value',...) creates a new MATRIXREADER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MatrixReader_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MatrixReader_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MatrixReader

% Last Modified by GUIDE v2.5 18-Nov-2011 11:33:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MatrixReader_OpeningFcn, ...
                   'gui_OutputFcn',  @MatrixReader_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MatrixReader is made visible.
function MatrixReader_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MatrixReader (see VARARGIN)
set(handles.pushbutton_ok,'Enable','off');

% Choose default command line output for MatrixReader
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MatrixReader wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MatrixReader_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_size as text
%        str2double(get(hObject,'String')) returns contents of edit_size as a double

set(handles.pushbutton_save,'Enable','off');
set(handles.pushbutton_ok,'Enable','off');
sizeh = get(gcbo,'String');
sizeh = str2num(sizeh);

if isempty(sizeh)
    warndlg({'Invalid Size value.','Using 3 '},'Invalid Value');
    set(gcbo,'String','3');
    sizeh = 3;
end

sizeh = round(abs((sizeh)));
if sizeh - floor(sizeh/2)*2 == 0,
    sizeh = sizeh +1;
    set(gcbo,'String',num2str(sizeh));

end

set(handles.uitable,'Data',zeros(sizeh));
set(handles.uitable,'Visible','on');
set(handles.pushbutton_save,'Enable','on');
set(handles.pushbutton_ok,'Enable','on');
%set(handles.uitable,'ColumnWidth',num2cell(25*ones(1,sizeh)));
handles.sizeh=sizeh;
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f_name,n_route,index]=uiputfile('*.mat','Save Matrix');
f_name_out ='h_aux.mat';
if index == 1
    h = get(handles.uitable,'Data');
    if sum(h(:)) ~= 0
         h = h/sum(h(:));
    end
    save([n_route f_name],'h');
    save('h_aux.mat','h');
   % save([n_route f_name_out],'h');
end


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f_name,n_route,index]=uigetfile({'*.mat'},'Load Matrix');

if index == 1
    load([n_route f_name]);
    [m,n] =size(h);
    if m ~= n
        warndlg({'Matrix must be square.'},'Invalid Value');
    else
        set(handles.uitable,'visible','on');
        set(handles.uitable,'Data',h);
        set(handles.edit_size,'String',m);
        handles.sizeh = m;
        set(handles.pushbutton_save,'Enable','on');

    end
end
set(handles.pushbutton_ok,'Enable','on');
handles.output = hObject;
guidata(hObject, handles);



% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = get(handles.uitable,'Data');
%size(h) 
%handles.sizeh
if h(:) == 0
 warndlg({'All the values are equal to zero.','Setting them to one.'},'Invalid Value');  
 h=ones(size(h));
  set(handles.uitable,'Data',h);
 
end

 if sum(h(:)) ~= 0
         h = h/sum(h(:));
         
 end
    
%h = h/sum(h(:));
save('tempSR/h_aux.mat','h');
%save('h_aux.mat','h');
handles.output = hObject;
guidata(hObject, handles);
close(handles.figure1);


% --- Executes during object creation, after setting all properties.
function uitable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%load('tempSR/h_aux.mat');
%set(gcbo,'Data',h);


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in uitable.
function uitable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
