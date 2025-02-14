function varargout = infoBox(varargin)
% infoBox - GUI to display information to inform the user with a picture
%
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
%
% Syntax:
%    varargout = infoBox(varargin)
%
% Inputs:
%    varargin - Cell array with one cell containing the path to 
%               the picture that should be displayed on the Info Box
%
% Outputs:
%    varargout - ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       ???
% Last update:   13-February-2015
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aux_infoBox_OpeningFcn, ...
                   'gui_OutputFcn',  @aux_infoBox_OutputFcn, ...
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


% Auxiliary functions -----------------------------------------------------

function aux_infoBox_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for infoBox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if nargin ~= 4
   errordlg('Wrong number of input arguments','Error Dialog','modal'); 
end

% position ok panel
set(handles.ok_panel,'Units','pixels');
ok_panel_pos = get(handles.ok_panel, 'Position');

set(hObject,'Units','pixels');
pos = get(hObject,'Position');

% position image
path = cell2mat(varargin{1,1});
h = imread(path);

if size(h,1) > 700 
    pos(1,3) = size(h,2)/1.5;
    pos(1,4) = size(h,1)/1.5;
else
    pos(1,3) = size(h,2);
    pos(1,4) = size(h,1);
end

set(hObject,'Position',pos+50);
set(hObject,'Units','characters');

% set axes
posAxes = pos;
posAxes(1,1) = 25;
posAxes(1,2) = ok_panel_pos(1,4)+5;
ok_panel_pos(1,3) = pos(1,3)+55;
posAxes(1,3) = posAxes(1,3);
posAxes(1,4) = posAxes(1,4);
set(handles.axes,'Units','pixels');
set(handles.axes,'Position',posAxes);

% set axes ok panel
set(handles.axes,'Units','characters');
set(handles.ok_panel,'Position',ok_panel_pos);
set(handles.ok_panel,'Units','characters');

% set axes pb
set(handles.pb,'Units','pixels');
posPb = get(handles.pb,'Position');
posPb(1,1) = round((pos(1,3) - posPb(1,3))/2)+20;
set(handles.pb,'Position',posPb);
set(handles.pb,'Units','characters');
imshow(h);

% turn off ticks
set(handles.axes,'XTick',[]);
set(handles.axes,'YTick',[]);


function varargout = aux_infoBox_OutputFcn(hObject, eventdata, handles) 
% output function
varargout{1} = handles.output;


function aux_pb_Callback(hObject, eventdata, handles)
    delete(handles.infoBox)


% --- Executes on slider movement.
function aux_slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function aux_slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% ------------------------------ END OF CODE ------------------------------
