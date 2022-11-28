function varargout = input_set(varargin)
% INPUTSET MATLAB code for InputSet.fig
%      INPUTSET, by itself, creates a new INPUTSET or raises the existing
%      singleton*.
%
%      H = INPUTSET returns the handle to a new INPUTSET or the handle to
%      the existing singleton*.
%
%      INPUTSET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUTSET.M with the given input arguments.
%
%      INPUTSET('Property','Value',...) creates a new INPUTSET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before input_set_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to input_set_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InputSet

% Last Modified by GUIDE v2.5 20-Dec-2020 19:54:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @input_set_OpeningFcn, ...
    'gui_OutputFcn',  @input_set_OutputFcn, ...
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


% --- Executes just before InputSet is made visible.
function input_set_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InputSet (see VARARGIN)

% Choose default command line output for InputSet
handles.output = hObject;

handles.previous_handles = varargin{1};
handles.current_handles = handles.previous_handles;
type = varargin{2};

if ~isempty(type)
    % Update handles structure
    if type == 'zonotope'
        handles.ZonotopeCenter = handles.previous_handles.ZonotopeCenter;
        set(handles.txtZonotopeCenter, 'String', handles.ZonotopeCenter)
        handles.ZonotopeGM = handles.previous_handles.ZonotopeGM;
        set(handles.txtZonotopeGM, 'String', handles.ZonotopeGM)
        set(handles.rbZonotope,'Value', 1)
        set(handles.rbInterval,'Value', 0)
        set(handles.txtIntervalCenter, 'Enable', 'off')
        set(handles.txtIntervalWidth, 'Enable', 'off')
        set(handles.popIntervalCenter, 'Enable', 'off')
        set(handles.popIntervalWidth, 'Enable', 'off')
        set(handles.txtZonotopeCenter, 'Enable', 'on')
        set(handles.txtZonotopeGM, 'Enable', 'on')
        set(handles.popZonotopeCenter, 'Enable', 'on')
        set(handles.popZonotopeGM, 'Enable', 'on')
    elseif type == 'interval'
        handles.IntervalCenter = handles.previous_handles.IntervalCenter;
        set(handles.txtIntervalCenter, 'String', handles.IntervalCenter)
        handles.IntervalWidth = handles.previous_handles.IntervalWidth;
        set(handles.txtIntervalWidth, 'String', handles.IntervalWidth)
        set(handles.rbInterval,'Value', 1)
        set(handles.rbZonotope,'Value', 0)
        set(handles.txtIntervalCenter, 'Enable', 'on')
        set(handles.txtIntervalWidth, 'Enable', 'on')
        set(handles.popIntervalCenter, 'Enable', 'on')
        set(handles.popIntervalWidth, 'Enable', 'on')
        set(handles.txtZonotopeCenter, 'Enable', 'off')
        set(handles.txtZonotopeGM, 'Enable', 'off')
        set(handles.popZonotopeCenter, 'Enable', 'off')
        set(handles.popZonotopeGM, 'Enable', 'off')
    end
end

WS_vars = evalin('base', 'who');
WS_vars = [{''}; WS_vars];

if isempty(WS_vars)
    set(handles.popIntervalCenter, 'String', 'No Worskpace Variables')
    set(handles.popIntervalWidth, 'String', 'No Worskpace Variables')
    set(handles.popZonotopeCenter, 'String', 'No Worskpace Variables')
    set(handles.popZonotopeGM, 'String', 'No Worskpace Variables')
else
    set(handles.popIntervalCenter, 'String', WS_vars)
    set(handles.popIntervalCenter, 'Value', 1)
    set(handles.popIntervalWidth, 'String', WS_vars)
    set(handles.popIntervalWidth, 'Value', 1)
    set(handles.popZonotopeCenter, 'String', WS_vars)
    set(handles.popZonotopeCenter, 'Value', 1)
    set(handles.popZonotopeGM, 'String', WS_vars)
    set(handles.popZonotopeGM, 'Value', 1)
end
guidata(hObject, handles);
uiwait

% UIWAIT makes InputSet wait for user response (see UIRESUME)
% uiwait(handles.InputSet);

% --- Outputs from this function are returned to the command line.
function varargout = input_set_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.current_handles;

delete(hObject);


% --- Executes on button press in rbInterval.
function rbInterval_Callback(hObject, eventdata, handles)
% hObject    handle to rbInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbInterval
set(hObject, 'Enable', 'off')
set(handles.rbZonotope, 'Enable', 'on')
set(handles.rbZonotope, 'Value', 0)

set(handles.txtIntervalCenter, 'Enable', 'on')
set(handles.txtIntervalWidth, 'Enable', 'on')
set(handles.popIntervalCenter, 'Enable', 'on')
set(handles.popIntervalWidth, 'Enable', 'on')
set(handles.txtZonotopeCenter, 'Enable', 'off')
set(handles.txtZonotopeGM, 'Enable', 'off')
set(handles.popZonotopeCenter, 'Enable', 'off')
set(handles.popZonotopeGM, 'Enable', 'off')


function txtIntervalCenter_Callback(hObject, eventdata, handles)
% hObject    handle to txtIntervalCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIntervalCenter as text
%        str2double(get(hObject,'String')) returns contents of txtIntervalCenter as a double

center = get(hObject, 'String');
handles.IntervalCenter = center;
handles.current_handles.IntervalCenter = handles.IntervalCenter;
handles.IntervalCenterMethod = 'txt';
guidata(hObject, handles);
set(handles.popIntervalCenter, 'Value', 1)


% --- Executes during object creation, after setting all properties.
function txtIntervalCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIntervalCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtIntervalWidth_Callback(hObject, eventdata, handles)
% hObject    handle to txtIntervalWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIntervalWidth as text
%        str2double(get(hObject,'String')) returns contents of txtIntervalWidth as a double
width = get(hObject, 'String');
handles.IntervalWidth = width;
handles.current_handles.IntervalWidth = handles.IntervalWidth;
handles.IntervalWidthMethod = 'txt';
guidata(hObject, handles);
set(handles.popIntervalWidth, 'Value', 1)


% --- Executes during object creation, after setting all properties.
function txtIntervalWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIntervalWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbZonotope.
function rbZonotope_Callback(hObject, eventdata, handles)
% hObject    handle to rbZonotope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbZonotope

set(hObject, 'Enable', 'off')
set(handles.rbInterval, 'Enable', 'on')
set(handles.rbInterval, 'Value', 0)

set(handles.txtIntervalCenter, 'Enable', 'off')
set(handles.txtIntervalWidth, 'Enable', 'off')
set(handles.popIntervalCenter, 'Enable', 'off')
set(handles.popIntervalWidth, 'Enable', 'off')
set(handles.txtZonotopeCenter, 'Enable', 'on')
set(handles.txtZonotopeGM, 'Enable', 'on')
set(handles.popZonotopeCenter, 'Enable', 'on')
set(handles.popZonotopeGM, 'Enable', 'on')


function txtZonotopeCenter_Callback(hObject, eventdata, handles)
% hObject    handle to txtZonotopeCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtZonotopeCenter as text
%        str2double(get(hObject,'String')) returns contents of txtZonotopeCenter as a double

center = get(hObject, 'String');
handles.ZonotopeCenter = center;
handles.current_handles.ZonotopeCenter = handles.ZonotopeCenter;
handles.ZonotopeCenterMethod = 'txt';
guidata(hObject, handles);
set(handles.popZonotopeCenter, 'Value', 1)


% --- Executes during object creation, after setting all properties.
function txtZonotopeCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZonotopeCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtZonotopeGM_Callback(hObject, eventdata, handles)
% hObject    handle to txtZonotopeGM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtZonotopeGM as text
%        str2double(get(hObject,'String')) returns contents of txtZonotopeGM as a double

GM = get(hObject, 'String');
handles.ZonotopeGM = GM;
handles.current_handles.ZonotopeGM = handles.ZonotopeGM;
handles.ZonotopeGMMethod = 'txt';
guidata(hObject, handles);
set(handles.popZonotopeGM, 'Value', 1)


% --- Executes during object creation, after setting all properties.
function txtZonotopeGM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZonotopeGM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbCancel.
function pbCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pbCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = {};
guidata(hObject, handles)
uiresume


% --- Executes on button press in pbOk.
function pbOk_Callback(hObject, eventdata, handles)
% hObject    handle to pbOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.rbInterval, 'Value')
    if isfield(handles, 'IntervalCenter') && isfield(handles, 'IntervalWidth')
        center_str = handles.IntervalCenter;
        center = evalin('base', center_str);
        
        width_str = handles.IntervalWidth;
        width = evalin('base', width_str);
        
        %check for dimensions
        if length(center) ~= length(width)
            uiwait(msgbox('Dimensions must agree', 'Error', 'error', 'modal'))
            return
        end
        U.center = center_str;
        U.width = width_str;
        U.type = 'interval';
    else
        uiwait(msgbox('Please enter all options for Interval', 'Error', 'error', 'modal'))
        return
    end
elseif get(handles.rbZonotope, 'Value')
    if isfield(handles, 'ZonotopeCenter') && isfield(handles, 'ZonotopeGM')
        center_str = handles.ZonotopeCenter;
        center = evalin('base', center_str);
        
        GM_str = handles.ZonotopeGM;
        GM = evalin('base', GM_str);
        
        %check for dimensions
        if length(center) ~= length(GM)
            uiwait(msgbox('Dimensions must agree', 'Error', 'error', 'modal'))
            return
        end
        U.center = center_str; 
        U.width = GM_str;
        U.type = 'zonotope';
    else
        uiwait(msgbox('Please enter all options for Zonotope', 'Error', 'error', 'modal'))
        return
    end
else
    uiwait(msgbox('please choose one of the two options: Interval or Zonotope', 'Error', 'error', 'modal'))
    return
end

try
    handles.output = U;
    guidata(hObject, handles)
    uiresume
catch
end


% --- Executes on selection change in popIntervalCenter.
function popIntervalCenter_Callback(hObject, eventdata, handles)
% hObject    handle to popIntervalCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popIntervalCenter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popIntervalCenter
contents = cellstr(get(hObject,'String'));
center = contents{get(hObject,'Value')};
center_value = evalin('base', center);
set(handles.txtIntervalCenter,'String', mat2str(center_value));
handles.IntervalCenterMethod = 'pop';
handles.IntervalCenter = mat2str(center_value);
handles.current_handles.IntervalCenter = handles.IntervalCenter;

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popIntervalCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popIntervalCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popIntervalWidth.
function popIntervalWidth_Callback(hObject, eventdata, handles)
% hObject    handle to popIntervalWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popIntervalWidth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popIntervalWidth
contents = cellstr(get(hObject,'String'));
width = contents{get(hObject,'Value')};
width_value = evalin('base', width);
set(handles.txtIntervalWidth,'String', mat2str(width_value));
handles.IntervalWidthMethod = 'pop';
handles.IntervalWidth = mat2str(width_value);
handles.current_handles.IntervalWidth = handles.IntervalWidth;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popIntervalWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popIntervalWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popZonotopeCenter.
function popZonotopeCenter_Callback(hObject, eventdata, handles)
% hObject    handle to popZonotopeCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popZonotopeCenter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popZonotopeCenter
contents = cellstr(get(hObject,'String'));
center = contents{get(hObject,'Value')};
center_value = evalin('base', center);
set(handles.txtZonotopeCenter,'String', mat2str(center_value));
handles.ZonotopeCenterMethod = 'pop';
handles.ZonotopeCenter = mat2str(center_value);
handles.current_handles.ZonotopeCenter = handles.ZonotopeCenter;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popZonotopeCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popZonotopeCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popZonotopeGM.
function popZonotopeGM_Callback(hObject, eventdata, handles)
% hObject    handle to popZonotopeGM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popZonotopeGM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popZonotopeGM
contents = cellstr(get(hObject,'String'));
GM = contents{get(hObject,'Value')};
GM_value = evalin('base', GM);
set(handles.txtZonotopeGM,'String', mat2str(GM_value));
handles.ZonotopeGMMethod = 'pop';
handles.ZonotopeGM = mat2str(GM_value);
handles.current_handles.ZonotopeGM = handles.ZonotopeGM;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popZonotopeGM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popZonotopeGM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_zonotope.
function pb_zonotope_Callback(hObject, eventdata, handles)
% hObject    handle to pb_zonotope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im = 'Info_zonotope.png';
infoBox({[path_im, im]});
uiwait;

% --- Executes on button press in pb_interval.
function pb_interval_Callback(hObject, eventdata, handles)
% hObject    handle to pb_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im = 'Info_Interval.png';
infoBox({[path_im, im]});
uiwait;
