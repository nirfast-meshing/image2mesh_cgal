function varargout = image2mesh_gui(varargin)
% IMAGE2MESH_GUI MATLAB code for image2mesh_gui.fig
%      IMAGE2MESH_GUI, by itself, creates a new IMAGE2MESH_GUI or raises the existing
%      singleton*.
%
%      H = IMAGE2MESH_GUI returns the handle to a new IMAGE2MESH_GUI or the handle to
%      the existing singleton*.
%
%      IMAGE2MESH_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE2MESH_GUI.M with the given input arguments.
%
%      IMAGE2MESH_GUI('Property','Value',...) creates a new IMAGE2MESH_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before image2mesh_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to image2mesh_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help image2mesh_gui

% Last Modified by GUIDE v2.5 13-Jun-2011 13:22:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @image2mesh_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @image2mesh_gui_OutputFcn, ...
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


% --- Executes just before image2mesh_gui is made visible.
function image2mesh_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to image2mesh_gui (see VARARGIN)

% Choose default command line output for image2mesh_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes image2mesh_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = image2mesh_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function cell_size_Callback(hObject, eventdata, handles)
% hObject    handle to cell_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_size as text
%        str2double(get(hObject,'String')) returns contents of cell_size as a double


% --- Executes during object creation, after setting all properties.
function cell_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cell_radius_edge_Callback(hObject, eventdata, handles)
% hObject    handle to cell_radius_edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_radius_edge as text
%        str2double(get(hObject,'String')) returns contents of cell_radius_edge as a double


% --- Executes during object creation, after setting all properties.
function cell_radius_edge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_radius_edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function facet_size_Callback(hObject, eventdata, handles)
% hObject    handle to facet_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of facet_size as text
%        str2double(get(hObject,'String')) returns contents of facet_size as a double


% --- Executes during object creation, after setting all properties.
function facet_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to facet_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function facet_angle_Callback(hObject, eventdata, handles)
% hObject    handle to facet_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of facet_angle as text
%        str2double(get(hObject,'String')) returns contents of facet_angle as a double


% --- Executes during object creation, after setting all properties.
function facet_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to facet_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function facet_distance_Callback(hObject, eventdata, handles)
% hObject    handle to facet_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of facet_distance as text
%        str2double(get(hObject,'String')) returns contents of facet_distance as a double


% --- Executes during object creation, after setting all properties.
function facet_distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to facet_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function specialregion_size_Callback(hObject, eventdata, handles)
% hObject    handle to specialregion_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of specialregion_size as text
%        str2double(get(hObject,'String')) returns contents of specialregion_size as a double


% --- Executes during object creation, after setting all properties.
function specialregion_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to specialregion_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function specialregion_label_Callback(hObject, eventdata, handles)
% hObject    handle to specialregion_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of specialregion_label as text
%        str2double(get(hObject,'String')) returns contents of specialregion_label as a double


% --- Executes during object creation, after setting all properties.
function specialregion_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to specialregion_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function infilename_Callback(hObject, eventdata, handles)
% hObject    handle to infilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of infilename as text
%        str2double(get(hObject,'String')) returns contents of infilename as a double


% --- Executes during object creation, after setting all properties.
function infilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to infilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browsebutton.
function browsebutton_Callback(hObject, eventdata, handles)
% hObject    handle to browsebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
