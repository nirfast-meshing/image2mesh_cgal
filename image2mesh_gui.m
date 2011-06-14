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

% Last Modified by GUIDE v2.5 14-Jun-2011 09:33:33

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
[fn, pathname] = uigetfile( ...
    {'*.bmp;*.jpg;*.tif;*.gif;*.mha','Image Files (*.bmp,*.jpg,*.tif,*.gif,*.mha)';'*.*','All Files (*.*)'}, ...
   'Pick a file');
if isequal(fn,0)
    error('\nYou need to select an image file!');
end
handles.inputfn = fullfile(pathname,fn);
guidata(hObject,handles);
UpdateInputFileInfo(hObject,eventdata,handles);

function UpdateInputFileInfo(hObject,eventdata,handles)
set(handles.infilename,'String',handles.inputfn);
UpdateImageInformation(hObject,eventdata,handles);
UpdateOutputFn(hObject,eventdata,handles);

function UpdateOutputFn(hObject,eventdata,handles)
s=get(handles.infilename,'String');
s=remove_extension(s);
s=add_extension(s,'.ele');
set(handles.outputfn,'String',s)

function UpdateImageInformation(hObject,eventdata,handles)
[mask info] = GetImageStack(get(handles.infilename,'String'),[]);
regions = unique(mask(:));
if ~isempty(info)
    s{1} = sprintf('Rows: %d, Cols: %d, Slices: %d',info.Dimensions(1), ...
        info.Dimensions(2), info.Dimensions(3));
    if isfield(info,'PixelDimensions')
        set(handles.xpixel,'String',num2str(info.PixelDimensions(1)));
        set(handles.ypixel,'String',num2str(info.PixelDimensions(2)));
        set(handles.zpixel,'String',num2str(info.PixelDimensions(3)));
        s{2} = sprintf('Offset: [%.1f %.1f %.1f]',info.Offset(1),info.Offset(2),info.Offset(3));
    else
    end
    set(handles.imageinfotxt,'String',s);
end
s={};
s{1} = sprintf('Region IDs (%d total):',length(regions));
fmt = ['%d,' repmat(' %d,',1,length(regions)-1)];
s{2} = sprintf(fmt,regions);
set(handles.specialregiontxt,'String',s);


% --- Executes on key press with focus on infilename and none of its controls.
function infilename_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to infilename (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key,'return')
    guidata(hObject,handles);
    UpdateInputFileInfo(hObject,eventdata,handles);
end







function xpixel_Callback(hObject, eventdata, handles)
% hObject    handle to xpixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xpixel as text
%        str2double(get(hObject,'String')) returns contents of xpixel as a double


% --- Executes during object creation, after setting all properties.
function xpixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xpixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ypixel_Callback(hObject, eventdata, handles)
% hObject    handle to ypixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ypixel as text
%        str2double(get(hObject,'String')) returns contents of ypixel as a double


% --- Executes during object creation, after setting all properties.
function ypixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ypixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zpixel_Callback(hObject, eventdata, handles)
% hObject    handle to zpixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zpixel as text
%        str2double(get(hObject,'String')) returns contents of zpixel as a double


% --- Executes during object creation, after setting all properties.
function zpixel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in callimage2mesh_cgal.
function callimage2mesh_cgal_Callback(hObject, eventdata, handles)
% hObject    handle to callimage2mesh_cgal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

param.facet_angle = str2double(get(handles.facet_angle,'String'));
param.facet_distance = str2double(get(handles.facet_distance,'String'));
param.facet_size = str2double(get(handles.facet_size,'String'));
param.medfilter = 0;
param.pad = 0;
param.cell_size = str2double(get(handles.cell_size,'String'));
param.cell_radius_edge = str2double(get(handles.cell_radius_edge,'String'));
param.special_subdomain_label = str2num(get(handles.specialregion_label,'String'));
param.special_subdomain_size  = str2num(get(handles.specialregion_size, 'String'));
s = get(handles.infilename,'String');
[foo ext]=remove_extension(s);
if ~strcmpi(ext,'.mha')
    sx=str2double(get(handles.xpixel,'String'));
    sy=str2double(get(handles.ypixel,'String'));
    sz=str2double(get(handles.zpixel,'String'));
    if isempty(sx) || isempty(sy) || isempty(sz) || ...
            isnan(sx) || isnan(sy) || isnan(sz)
        errordlg('Pixel size information is invalid!');
        error('Pixel size information is invalid!');
    end
    param.xpixelsize = sx;
    param.ypixelsize = sy;
    param.zpixelsize = sz;
end
h=helpdlg('Please wait...','Mesh Generator Running!');
image2mesh_cgal(s,param,get(handles.outputfn,'String'));
if ishandle(h)
    close(h);
end

function outputfn_Callback(hObject, eventdata, handles)
% hObject    handle to outputfn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputfn as text
%        str2double(get(hObject,'String')) returns contents of outputfn as a double


% --- Executes during object creation, after setting all properties.
function outputfn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputfn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
