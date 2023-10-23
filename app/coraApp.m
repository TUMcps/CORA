

% Table of Content:
% - Initialization code
% - General functions
% - Linear System: 
%       - Custom Functions - Linear
%       - Callback Functions - Linear
%       - Finish Linear
% - Nonlinear System:
%       - Custom Functions - Nonlinear
%       - Callback Functions - Nonlinear 
% - Hybrid System:
%       - Custom Functions - Hybrid
%       - Callback Functions - Hybrid


% --- Initialization code -------------------------------------------------                         


% Initialization function
function varargout = coraApp(varargin)
% coraApp MATLAB code for coraApp.fig
%      coraApp, by itself, creates a new coraApp or raises the existing
%      singleton*.
%
%      H = coraApp returns the handle to a new coraApp or the handle to
%      the existing singleton*.
%
%      coraApp('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in coraApp.M with the given input arguments.
%
%      coraApp('Property','Value',...) creates a new coraApp or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before coraApp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coraApp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help coraApp

% Last Modified by GUIDE v2.5 15-Jul-2021 18:09:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @coraApp_OpeningFcn, ...
    'gui_OutputFcn',  @coraApp_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargin > 1 && ismember( ...
        varargin{1},{ ...
        'txtIC_RSim_Hybrid_CreateFcn' ...
        'txtIC_RSim_Linear_CreateFcn' ...
        'txtIC_RSim_Nonlinear_CreateFcn'...
        })
    % TL: not sure what these are. Throw an error when gui_mainfcn is
    % called but they don't seem to be required anyway (?), skipping ...
else
    % call main func
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});

    % clear output
    if nargout == 0
        clear varargout
    end
end

% Opening function --- Executes before coraApp is made visible
function coraApp_OpeningFcn(hObject, eventdata, handles, varargin)

handles.dynamicsEq_Nonlinear_error = 0;
handles.txt_equation = {};
handles.output = hObject;

handles.linear_plots = {};
handles.nonlinear_plots = {};
handles.hybrid_plots = {};

handles.dimensions_Linear = {};
handles.dimensions_Nonlinear = {};
handles.dimensions_Hybrid = {};

handles.nrOfInputs_Linear = {};
handles.nrOfInputs_Nonlinear = {};
handles.nrOfInputs_Hybrid = {};

handles.initial_set_linear = {};
handles.initial_set_linear_type = {};
handles.input_set_linear = {};
handles.input_set_linear_type = {};
handles.initial_set_nonlinear = {};
handles.initial_set_nonlinear_type = {};
handles.input_set_nonlinear = {};
handles.input_set_nonlinear_type = {};
handles.initial_set_hybrid = {};
handles.initial_set_hybrid_type = {};
handles.input_set_hybrid = {};
handles.input_set_hybrid_type = {};

% TL: some control handles are missing (?), init dummy handles ..
dummyUIControl = uicontrol(gcf,'String','10');
set(dummyUIControl,'Visible','off');    
handles.txtNCI_RSim_Linear = dummyUIControl;
handles.txtNCI_RSim_Nonlinear = dummyUIControl;
handles.txtNCI_RSim_Hybrid = dummyUIControl;
handles.pbUpNCI_RSim_Linear = dummyUIControl;
handles.pbDownNCI_RSim_Linear = dummyUIControl;
handles.pbUpNCI_RSim_Nonlinear = dummyUIControl;
handles.pbDownNCI_RSim_Nonlinear = dummyUIControl;
handles.pbUpNCI_RSim_Hybrid = dummyUIControl;
handles.pbDownNCI_RSim_Hybrid = dummyUIControl;

WS_vars = evalin('base', 'who');
WS_vars = [{''}; WS_vars];

%show the linear panel only in the begining
set(handles.panelLinear, 'visible', 'on')
set(handles.panelNonLinear, 'visible', 'off')
set(handles.panelHybrid, 'visible', 'off')

%show main systems
main_systems  = {'Linear System', 'Nonlinear System', 'Hybrid System'};

set(handles.popupMain, 'String', main_systems);
set(handles.popupMain, 'Value', 1);

%show workspace variables
if isempty(WS_vars)
    %linear
    set(handles.popUTrans_Linear, 'String', 'No Workspace Variables')
    set(handles.popX0_Sim_Linear, 'String', 'No Workspace Variables')
    set(handles.popU_Sim_Linear, 'String', 'No Workspace Variables')
    set(handles.popA_Linear, 'String', 'No Workspace Variables')
    set(handles.popB_Linear, 'String', 'No Workspace Variables')
    set(handles.popC_Linear, 'String', 'No Workspace Variables')
    set(handles.popD_Linear, 'String', 'No Workspace Variables')
    set(handles.popE_Linear, 'String', 'No Workspace Variables')
    set(handles.popF_Linear, 'String', 'No Workspace Variables')
    
    %nonlinear
    set(handles.popUTrans_Nonlinear, 'String', 'No Workspace Variables')
    set(handles.popX0_Sim_Nonlinear, 'String', 'No Workspace Variables')
    set(handles.popU_Sim_Nonlinear, 'String', 'No Workspace Variables')
    
    %Hybrid
    set(handles.popX0_Sim_Hybrid, 'String', 'No Workspace Variables')
    set(handles.popU_Sim_Hybrid, 'String', 'No Workspace Variables')
    
else
    %linear
    set(handles.popUTrans_Linear, 'String', WS_vars)
    set(handles.popUTrans_Linear, 'Value', 1)
    
    set(handles.popX0_Sim_Linear, 'String', WS_vars)
    set(handles.popX0_Sim_Linear, 'Value', 1)
    set(handles.popU_Sim_Linear, 'String', WS_vars)
    set(handles.popU_Sim_Linear, 'Value', 1)
    
    set(handles.popA_Linear, 'String', WS_vars)
    set(handles.popA_Linear, 'Value', 1)
    set(handles.popB_Linear, 'String', WS_vars)
    set(handles.popB_Linear, 'Value', 1)
    set(handles.popC_Linear, 'String', WS_vars)
    set(handles.popC_Linear, 'Value', 1)
    set(handles.popD_Linear, 'String', WS_vars)
    set(handles.popD_Linear, 'Value', 1)
    set(handles.popE_Linear, 'String', WS_vars)
    set(handles.popE_Linear, 'Value', 1)
    set(handles.popF_Linear, 'String', WS_vars)
    set(handles.popF_Linear, 'Value', 1)
    
    %nonlinear
    set(handles.popUTrans_Nonlinear, 'String', WS_vars)
    set(handles.popUTrans_Nonlinear, 'Value', 1)
    
    set(handles.popX0_Sim_Nonlinear, 'String', WS_vars)
    set(handles.popX0_Sim_Nonlinear, 'Value', 1)
    set(handles.popU_Sim_Nonlinear, 'String', WS_vars)
    set(handles.popU_Sim_Nonlinear, 'Value', 1)
    
    %hybrid
    set(handles.popX0_Sim_Hybrid, 'String', WS_vars)
    set(handles.popX0_Sim_Hybrid, 'Value', 1)
    set(handles.popU_Sim_Hybrid, 'String', WS_vars)
    set(handles.popU_Sim_Hybrid, 'Value', 1)
end

%show logo
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
imfile1 = 'CORAlogo.png';
axes(handles.panelHeader, 'position', [0.6,0.3,0.6,0.7]);
[cora_logo, ~, alpha] = imread([path_im, imfile1]);
image(cora_logo, 'AlphaData', alpha)
axis off
axis image

% refresh icon button
[a,~]=imread('Refresh.jpg');
[r,c,~]=size(a); 
x=ceil(r/35); 
y=ceil(c/35); 
g=a(1:x:end,1:y:end,:);
g(g==255)=5.5*255;
set(handles.refresh,'CData',g)

%show equations
axes(handles.panelSystem_Linear, 'position', [0.35,-0.2,0.209,0.7]);
imfile2 = 'equation.png';
[eq_img, map] = imread([path_im, imfile2]);
alpha = ones(size(eq_img,1), size(eq_img,2));
alpha(eq_img(:,:,1) == 255 & eq_img(:,:,2) == 255 & eq_img(:,:,3) == 255) = 0;

h = imshow(eq_img, map);
set(h, 'AlphaData', alpha)
axis off
axis image

%display algorithm in linear
alg = {'standard', 'wrapping-free', 'fromStart', 'adaptive'};
set(handles.popLA_Linear, 'String', alg);
set(handles.popLA_Linear, 'Value', 4);

set(handles.txtZOrder_Linear, 'visible', 'off')
set(handles.pbUpZOrder_Linear, 'visible', 'off')
set(handles.pbDownZOrder_Linear, 'visible', 'off')
set(handles.txtTaylor_Linear, 'visible', 'off')
set(handles.pbUpTaylor_Linear, 'visible', 'off')
set(handles.pbDownTaylor_Linear, 'visible', 'off')
set(handles.popRT_Linear, 'visible', 'off')
set(handles.zotxt, 'visible', 'off')
set(handles.tttxt, 'visible', 'off')
set(handles.rttxt, 'visible', 'off')
set(handles.pb_infoZO_Linear, 'visible', 'off')
set(handles.pb_infoTT_Linear, 'visible', 'off')
set(handles.pb_infoRT_Linear, 'visible', 'off')
set(handles.txtTimeStep_Linear, 'visible', 'off')
set(handles.tstxt, 'visible', 'off')

%display algorithm in Hybrid
alg = {'standard', 'wrapping-free', 'fromStart', 'adaptive'};
set(handles.popLA_Hybrid, 'String', alg);
set(handles.popLA_Hybrid, 'Value', 4);

set(handles.OptionsError_Hybrid, 'visible', 'off')
set(handles.txtOptionError_Hybrid, 'visible', 'off')
set(handles.txtOptionErrorUp_Hybrid, 'visible', 'off')
set(handles.txtOptionErrorDown_Hybrid, 'visible', 'off')
set(handles.txtZOrder_Hybrid, 'visible', 'off')
set(handles.pbUpZOrder_Hybrid, 'visible', 'off')
set(handles.pbDownZOrder_Hybrid, 'visible', 'off')
set(handles.txtTaylor_Hybrid, 'visible', 'off')
set(handles.pbUpTaylor_Hybrid, 'visible', 'off')
set(handles.pbDownTaylor_Hybrid, 'visible', 'off')
set(handles.popRT_Hybrid, 'visible', 'off')
set(handles.zotxt_Hybrid, 'visible', 'off')
set(handles.tttxt_Hybrid, 'visible', 'off')
set(handles.rttxt_Hybrid, 'visible', 'off')
set(handles.pb_infoZO_Hybrid, 'visible', 'off')
set(handles.pb_infoTT_Hybrid, 'visible', 'off')
set(handles.pb_infoRT_Hybrid, 'visible', 'off')
set(handles.txtTimeStep_Hybrid, 'visible', 'off')
set(handles.tstxt_Hybrid, 'visible', 'off')
set(handles.txtGO_Hybrid, 'Enable', 'off')
set(handles.GOtxt_Hybrid, 'Enable', 'off')

%linear
set(handles.popLA_Hybrid, 'visible', 'off')
set(handles.txtA_Hybrid, 'visible', 'off')
set(handles.pbInfoLA_Hybrid, 'visible', 'off')

%nonlinear
set(handles.txtTO_Hybrid, 'visible', 'off')
set(handles.txtTO2_Hybrid, 'visible', 'off')
set(handles.txtIO_Hybrid, 'visible', 'off')
set(handles.txtIO2_Hybrid, 'visible', 'off')
set(handles.txtEO_Hybrid, 'visible', 'off')
set(handles.txtEO2_Hybrid, 'visible', 'off')
set(handles.txtCP_Hybrid, 'visible', 'off')

set(handles.rbConLin_Hybrid, 'visible', 'off')
set(handles.txtConLinTO_Hybrid, 'visible', 'off')
set(handles.pbUpConLinTO_Hybrid, 'visible', 'off')
set(handles.pbDownConLinTO_Hybrid, 'visible', 'off')
set(handles.txtIO_Hybrid_lin, 'visible', 'off')
set(handles.txtEO_Hybrid_lin, 'visible', 'off')
set(handles.pbUpIO_Hybrid_lin, 'visible', 'off')
set(handles.pbDownIO_Hybrid_lin, 'visible', 'off')
set(handles.pbUpEO_Hybrid_lin, 'visible', 'off')
set(handles.pbDownEO_Hybrid_lin, 'visible', 'off')
set(handles.pb_infoCPIO_Hybrid_lin, 'visible', 'off')
set(handles.pb_infoCPEE_Hybrid_lin, 'visible', 'off')
set(handles.txtIO2_Hybrid_lin, 'visible', 'off')
set(handles.txtEO2_Hybrid_lin, 'visible', 'off')

set(handles.pb_infoCL_Hybrid, 'visible', 'off')
set(handles.pb_infoCLTO_Hybrid, 'visible', 'off')
set(handles.pb_infoCP_Hybrid, 'visible', 'off')
set(handles.pb_infoCPTO_Hybrid, 'visible', 'off')
set(handles.pb_infoCPIO_Hybrid, 'visible', 'off')
set(handles.pb_infoCPEE_Hybrid, 'visible', 'off')
set(handles.txtConPolTO_Hybrid, 'visible', 'off')
set(handles.pbUpConPolTO_Hybrid, 'visible', 'off')
set(handles.pbDownConPolTO_Hybrid, 'visible', 'off')
set(handles.pbUpIO_Hybrid, 'visible', 'off')
set(handles.pbDownIO_Hybrid, 'visible', 'off')
set(handles.pbUpEO_Hybrid, 'visible', 'off')
set(handles.pbDownEO_Hybrid, 'visible', 'off')

%display plot colors
colors = {'none', 'red', 'green', 'blue', 'black', 'yellow', 'magenta', 'cyan', 'white', 'gray'};

%linear
set(handles.popColorReach_Linear, 'String', colors);
set(handles.popColorReach_Linear, 'Value', 10);

set(handles.popEdgeColorReach_Linear, 'String', colors);
set(handles.popEdgeColorReach_Linear, 'Value', 1);

set(handles.popColorInitial_Linear, 'String', colors);
set(handles.popColorInitial_Linear, 'Value', 9);

set(handles.popEdgeColorInitial_Linear, 'String', colors);
set(handles.popEdgeColorInitial_Linear, 'Value', 5);

set(handles.popColorSimulation_Linear, 'String', colors);
set(handles.popColorSimulation_Linear, 'Value', 5);

%nonlinear
set(handles.popColorReach_Nonlinear, 'String', colors);
set(handles.popColorReach_Nonlinear, 'Value', 10);

set(handles.popEdgeColorReach_Nonlinear, 'String', colors);
set(handles.popEdgeColorReach_Nonlinear, 'Value', 1);

set(handles.popColorInitial_Nonlinear, 'String', colors);
set(handles.popColorInitial_Nonlinear, 'Value', 9);

set(handles.popEdgeColorInitial_Nonlinear, 'String', colors);
set(handles.popEdgeColorInitial_Nonlinear, 'Value', 5);

set(handles.popColorSimulation_Nonlinear, 'String', colors);
set(handles.popColorSimulation_Nonlinear, 'Value', 5);

%hybrid
set(handles.popColorReach_Hybrid, 'String', colors);
set(handles.popColorReach_Hybrid, 'Value', 10);

set(handles.popEdgeColorReach_Hybrid, 'String', colors);
set(handles.popEdgeColorReach_Hybrid, 'Value', 1);

set(handles.popColorInitial_Hybrid, 'String', colors);
set(handles.popColorInitial_Hybrid, 'Value', 9);

set(handles.popEdgeColorInitial_Hybrid, 'String', colors);
set(handles.popEdgeColorInitial_Hybrid, 'Value', 5);

set(handles.popColorSimulation_Hybrid, 'String', colors);
set(handles.popColorSimulation_Hybrid, 'Value', 5);

%display plot line styles
line_styles = ...
    {'-', '-+', '-o', '-*', '-.', '-x', '-s', '-d', '-^', '->', '-<', '-p', '-h', ...
    ':+', ':o', ':*', ':.', ':x', ':s', ':d', ':^', ':>', ':<', ':p', ':h', ...
    '--+', '--', '--*', '--.', '--x', '--s', '--d', '--^', '-->', '--<', '--p', '--h', ...
    '-.+', '-.*', '-..', '-.x', '-.s', '-.d', '-.^', '-.>', '-.<', '-.p', '-.h'};

%linear
set(handles.popLineSimulation_Linear, 'String', line_styles);
set(handles.popLineSimulation_Linear, 'Value', 1);

%nonlinear
set(handles.popLineSimulation_Nonlinear, 'String', line_styles);
set(handles.popLineSimulation_Nonlinear, 'Value', 1);

%hybrid
set(handles.popLineSimulation_Hybrid, 'String', line_styles);
set(handles.popLineSimulation_Hybrid, 'Value', 1);

% Update handles structure
guidata(hObject, handles);

% Output Function --- Outputs are returned to the command line
function varargout = coraApp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function popupMain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTE
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupMain.
function popupMain_Callback(hObject, eventdata, handles)
% hObject    handle to popupMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupMain contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMain

contents = cellstr(get(hObject,'String'));
currentOption = contents{get(hObject,'Value')};

switch currentOption
    case 'Linear System'
        set(handles.panelLinear, 'visible', 'on')
        set(handles.panelNonLinear, 'visible', 'off')
        set(handles.panelHybrid, 'visible', 'off')
        set(handles.panelSystem_Linear, 'visible', 'on')
        set(handles.panelSettings_Linear, 'visible', 'off')
        set(handles.panelSimulation_Linear, 'visible', 'off')
        set(handles.panelPlots_Linear, 'visible', 'off')
        set(handles.pbSystem, 'Enable', 'off')
        set(handles.pbSimulation, 'Enable', 'on')
        set(handles.pbSettings, 'Enable', 'on')
        set(handles.pbPlots, 'Enable', 'on')
        
    case 'Nonlinear System'
        set(handles.panelLinear, 'visible', 'off')
        set(handles.panelNonLinear, 'visible', 'on')
        set(handles.panelHybrid, 'visible', 'off')
        set(handles.panelSystem_Nonlinear, 'visible', 'on')
        set(handles.panelSettings_Nonlinear, 'visible', 'off')
        set(handles.panelSimulation_Nonlinear, 'visible', 'off')
        set(handles.panelPlots_Nonlinear, 'visible', 'off')
        set(handles.pbSystem, 'Enable', 'off')
        set(handles.pbSimulation, 'Enable', 'on')
        set(handles.pbSettings, 'Enable', 'on')
        set(handles.pbPlots, 'Enable', 'on')
        
    case 'Hybrid System'
        set(handles.panelLinear, 'visible', 'off')
        set(handles.panelNonLinear, 'visible', 'off')
        set(handles.panelHybrid, 'visible', 'on')
        set(handles.panelSystem_Hybrid, 'visible', 'on')
        set(handles.panelSettings_Hybrid, 'visible', 'off')
        set(handles.panelSimulation_Hybrid, 'visible', 'off')
        set(handles.panelPlots_Hybrid, 'visible', 'off')
        set(handles.pbSystem, 'Enable', 'off')
        set(handles.pbSimulation, 'Enable', 'on')
        set(handles.pbSettings, 'Enable', 'on')
        set(handles.pbPlots, 'Enable', 'on')
end


% --- Executes on button press in pbRun.
function pbRun_Callback(hObject, eventdata, handles)

contents = cellstr(get(handles.popupMain,'String'));
currentOption = contents{get(handles.popupMain,'Value')};

switch currentOption
    case 'Linear System'
        script = sprintf('%s.m', 'Linear');
    case 'Nonlinear System'
        script = sprintf('%s.m', 'Nonlinear');
    case 'Hybrid System'
        script = sprintf('%s.m', 'Hybrid');
end

file = script;
path = [CORAROOT, filesep, 'models', filesep];
dinfo = dir(path);
dinfo(ismember( {dinfo.name}, {'.', '..'})) = [];
allFileNames = {dinfo(:).name};

if ~any(ismember(allFileNames, 'auxiliary'))
    mkdir(fullfile(path,'auxiliary'))
end

addpath([path, 'auxiliary'])
path = [CORAROOT, filesep, 'models', filesep, 'auxiliary'];
file_system = fullfile(path, file);
handles.fileSystem = file_system;
guidata(hObject, handles) %update handles

%generate the linear file to execute
[id, err] = generate_file(handles.fileSystem, handles, currentOption,{});

if id == -1
    handles = rmfield(handles, 'fileSystem');
    guidata(hObject, handles)
    return
else
    if err == 1
        return
    else
        run(handles.fileSystem)
    end
end


% --- Execues on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
%if get(handles.rbLoadSModel_Linear, 'Value')

contents = cellstr(get(handles.popupMain,'String'));
currentOption = contents{get(handles.popupMain,'Value')};

switch currentOption
    case 'Linear System'
        script = sprintf('%s.m', 'Linear');
    case 'Nonlinear System'
        script = sprintf('%s.m', 'Nonlinear');
    case 'Hybrid System'
        script = sprintf('%s.m', 'Hybrid');
end

[file, path] = uiputfile(script);
file_system = fullfile(path, file);
handles.fileSystem = file_system;
guidata(hObject, handles) %update handles

%generate the linear file to execute
[id, err] = generate_file(handles.fileSystem, handles, currentOption,file);
if id == -1
    handles = rmfield(handles, 'fileSystem');
    guidata(hObject, handles)
    return
end

if err == 1
    return
end



% --- Executes on button press in pbSystem.
function pbSystem_Callback(hObject, eventdata, handles)
% hObject    handle to pbSystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.popupMain,'String'));
currentOption = contents{get(handles.popupMain,'Value')};

switch currentOption
    case 'Linear System'
        set(handles.panelSystem_Linear, 'visible', 'on')
        set(handles.panelSettings_Linear, 'visible', 'off')
        set(handles.panelSimulation_Linear, 'visible', 'off')
        set(handles.panelPlots_Linear, 'visible', 'off')
        
    case 'Nonlinear System'
        set(handles.panelSystem_Nonlinear, 'visible', 'on')
        set(handles.panelSettings_Nonlinear, 'visible', 'off')
        set(handles.panelSimulation_Nonlinear, 'visible', 'off')
        set(handles.panelPlots_Nonlinear, 'visible', 'off')
        
    case 'Hybrid System'
        set(handles.panelSystem_Hybrid, 'visible', 'on')
        set(handles.panelSettings_Hybrid, 'visible', 'off')
        set(handles.panelSimulation_Hybrid, 'visible', 'off')
        set(handles.panelPlots_Hybrid, 'visible', 'off')
end

set(hObject, 'Enable', 'off')
set(handles.pbSimulation, 'Enable', 'on')
set(handles.pbSettings, 'Enable', 'on')
set(handles.pbPlots, 'Enable', 'on')


function pbSettings_Callback(hObject, eventdata, handles)
% hObject    handle to pbSimulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.popupMain,'String'));
currentOption = contents{get(handles.popupMain,'Value')};

switch currentOption
    case 'Linear System'
        set(handles.panelSystem_Linear, 'visible', 'off')
        set(handles.panelSettings_Linear, 'visible', 'on')
        set(handles.panelSimulation_Linear, 'visible', 'off')
        set(handles.panelPlots_Linear, 'visible', 'off')
        
    case 'Nonlinear System'
        set(handles.panelSystem_Nonlinear, 'visible', 'off')
        set(handles.panelSettings_Nonlinear, 'visible', 'on')
        set(handles.panelSimulation_Nonlinear, 'visible', 'off')
        set(handles.panelPlots_Nonlinear, 'visible', 'off')
        
    case 'Hybrid System'
        set(handles.panelSystem_Hybrid, 'visible', 'off')
        set(handles.panelSettings_Hybrid, 'visible', 'on')
        set(handles.panelSimulation_Hybrid, 'visible', 'off')
        set(handles.panelPlots_Hybrid, 'visible', 'off')
end

set(hObject, 'Enable', 'off')
set(handles.pbSystem, 'Enable', 'on')
set(handles.pbSimulation, 'Enable', 'on')
set(handles.pbPlots, 'Enable', 'on')


function pbSimulation_Callback(hObject, eventdata, handles)
% hObject    handle to pbSimulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.popupMain,'String'));
currentOption = contents{get(handles.popupMain,'Value')};

switch currentOption
    case 'Linear System'
        set(handles.panelSystem_Linear, 'visible', 'off')
        set(handles.panelSettings_Linear, 'visible', 'off')
        set(handles.panelSimulation_Linear, 'visible', 'on')
        set(handles.panelPlots_Linear, 'visible', 'off')
        
    case 'Nonlinear System'
        set(handles.panelSystem_Nonlinear, 'visible', 'off')
        set(handles.panelSettings_Nonlinear, 'visible', 'off')
        set(handles.panelSimulation_Nonlinear, 'visible', 'on')
        set(handles.panelPlots_Nonlinear, 'visible', 'off')
        
    case 'Hybrid System'
        set(handles.panelSystem_Hybrid, 'visible', 'off')
        set(handles.panelSettings_Hybrid, 'visible', 'off')
        set(handles.panelSimulation_Hybrid, 'visible', 'on')
        set(handles.panelPlots_Hybrid, 'visible', 'off')
end

set(hObject, 'Enable', 'off')
set(handles.pbSystem, 'Enable', 'on')
set(handles.pbSettings, 'Enable', 'on')
set(handles.pbPlots, 'Enable', 'on')


function pbPlots_Callback(hObject, eventdata, handles)
% hObject    handle to pbPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.popupMain,'String'));
currentOption = contents{get(handles.popupMain,'Value')};

switch currentOption
    case 'Linear System'
        set(handles.panelSystem_Linear, 'visible', 'off')
        set(handles.panelSettings_Linear, 'visible', 'off')
        set(handles.panelSimulation_Linear, 'visible', 'off')
        set(handles.panelPlots_Linear, 'visible', 'on')
        
    case 'Nonlinear System'
        set(handles.panelSystem_Nonlinear, 'visible', 'off')
        set(handles.panelSettings_Nonlinear, 'visible', 'off')
        set(handles.panelSimulation_Nonlinear, 'visible', 'off')
        set(handles.panelPlots_Nonlinear, 'visible', 'on')
        
    case 'Hybrid System'
        set(handles.panelSystem_Hybrid, 'visible', 'off')
        set(handles.panelSettings_Hybrid, 'visible', 'off')
        set(handles.panelSimulation_Hybrid, 'visible', 'off')
        set(handles.panelPlots_Hybrid, 'visible', 'on')
end

set(hObject, 'Enable', 'off')
set(handles.pbSystem, 'Enable', 'on')
set(handles.pbSimulation, 'Enable', 'on')
set(handles.pbSettings, 'Enable', 'on')


% --- Linear System -------------------------------------------------------


% --- Custom Functions - Linear -------------------------------------------


function matrices = create_matrices_linear(handles)

if get(handles.checkA_Linear, 'Value')
    mat_str = handles.A;
    A = evalin('base', mat_str);
    matrices.A = A;
end

if get(handles.checkB_Linear, 'Value')
    mat_str = handles.B;
    B = evalin('base', mat_str);
    matrices.B = B;
end

if get(handles.checkC_Linear, 'Value')
    mat_str = handles.C;
    C = evalin('base', mat_str);
    matrices.c = C;
end

if get(handles.checkD_Linear, 'Value')
    mat_str = handles.D;
    D = evalin('base', mat_str);
    matrices.C = D;
end

if get(handles.checkE_Linear, 'Value')
    mat_str = handles.E;
    E = evalin('base', mat_str);
    matrices.D = E;
end

if get(handles.checkF_Linear, 'Value')
    mat_str = handles.F;
    F = evalin('base', mat_str);
    matrices.K = F;
end



% --- Callback Functions - Linear -----------------------------------------


% --- Executes on button press in rbLoadSModel_Linear.
function rbLoadSModel_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to rbLoadSModel_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbLoadSModel_Linear
handles.dimensions_Linear = [];
guidata(hObject, handles)

set(handles.popDim1Plot_Linear, 'String', ' ');
set(handles.popDim1Plot_Linear, 'Value', 1);
set(handles.popDim2Plot_Linear, 'String', ' ');
set(handles.popDim2Plot_Linear, 'Value', 1);

if get(hObject, 'Value')
    set(handles.pbLoadSModel_Linear, 'Enable', 'on')
    set(hObject, 'Enable', 'off')
    set(handles.rbCreateMat_Linear, 'Enable', 'on')
    set(handles.rbCreateMat_Linear, 'Value', 0)
    
    set(handles.checkA_Linear, 'Enable', 'off')
    set(handles.checkB_Linear, 'Enable', 'off')
    set(handles.checkC_Linear, 'Enable', 'off')
    set(handles.checkD_Linear, 'Enable', 'off')
    set(handles.checkE_Linear, 'Enable', 'off')
    set(handles.checkF_Linear, 'Enable', 'off')
    
    set(handles.checkA_Linear, 'Value', 0)
    set(handles.checkB_Linear, 'Value', 0)
    set(handles.checkC_Linear, 'Value', 0)
    set(handles.checkD_Linear, 'Value', 0)
    set(handles.checkE_Linear, 'Value', 0)
    set(handles.checkF_Linear, 'Value', 0)
    
    set(handles.textA_Linear, 'Enable', 'off')
    set(handles.textB_Linear, 'Enable', 'off')
    set(handles.textC_Linear, 'Enable', 'off')
    set(handles.textD_Linear, 'Enable', 'off')
    set(handles.textE_Linear, 'Enable', 'off')
    set(handles.textF_Linear, 'Enable', 'off')
    
    set(handles.textA_Linear, 'String', ' ')
    set(handles.textB_Linear, 'String', ' ')
    set(handles.textC_Linear, 'String', ' ')
    set(handles.textD_Linear, 'String', ' ')
    set(handles.textE_Linear, 'String', ' ')
    set(handles.textF_Linear, 'String', ' ')
    
    set(handles.popA_Linear, 'Enable', 'off')
    set(handles.popB_Linear, 'Enable', 'off')
    set(handles.popC_Linear, 'Enable', 'off')
    set(handles.popD_Linear, 'Enable', 'off')
    set(handles.popE_Linear, 'Enable', 'off')
    set(handles.popF_Linear, 'Enable', 'off')
end


% --- Executes on button press in rbCreateMat_Linear.
function rbCreateMat_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to rbCreateMat_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbCreateMat_Linear
set(handles.pbLoadSModel_Linear, 'string', 'Load File')
handles.dimensions_Linear = [];
handles.nrOfInputs_Linear = [];
guidata(hObject, handles)

set(handles.popDim1Plot_Linear, 'String', ' ');
set(handles.popDim1Plot_Linear, 'Value', 1);
set(handles.popDim2Plot_Linear, 'String', ' ');
set(handles.popDim2Plot_Linear, 'Value', 1);
    
if get(hObject, 'Value')
    set(hObject, 'Enable', 'off')
    set(handles.rbLoadSModel_Linear, 'Enable', 'on')
    set(handles.rbLoadSModel_Linear, 'Value', 0)
    set(handles.pbLoadSModel_Linear, 'Enable', 'off')
    set(handles.checkA_Linear, 'Enable', 'on')
    set(handles.checkB_Linear, 'Enable', 'on')
    set(handles.checkC_Linear, 'Enable', 'on')
    set(handles.checkD_Linear, 'Enable', 'on')
    set(handles.checkE_Linear, 'Enable', 'on')
    set(handles.checkF_Linear, 'Enable', 'on')
end


% --- Executes on button press in pbLoadSModel_Linear.
function pbLoadSModel_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadSModel_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[spaceEx_file, spaceEx_path] = uigetfile('*.xml');
try
    handles.spaceExFile = fullfile(spaceEx_path, spaceEx_file);
    spaceex2cora(handles.spaceExFile);
catch
    return
end

file = strsplit(handles.spaceExFile, filesep);
file = file{end};
file = strsplit(file, '.');
file = file{1};
linSys =  eval(strcat(file, '()'));

if strcmp(class(linSys),'linearSys')
    dimensions = string(1:linSys.dim);
    dimensions1 = ['', dimensions, 'time'];
    dimensions2 = ['', dimensions];

    handles.dimensions_Linear = size(linSys.A,1);
    handles.nrOfInputs_Linear = linSys.nrOfInputs;
    guidata(hObject, handles)

    set(handles.popDim1Plot_Linear, 'String', dimensions1);
    set(handles.popDim1Plot_Linear, 'Value', 1);
    set(handles.popDim2Plot_Linear, 'String', dimensions2);
    set(handles.popDim2Plot_Linear, 'Value', 1);
    set(hObject, 'string', spaceEx_file)
    guidata(hObject, handles) %update handles
else
    uiwait(msgbox('The loaded SpaceEx file is not linear', 'Error', 'error', 'modal'))
    return
end

% --- Executes during object creation, after setting all properties.
function textA_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function textA_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to textA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = get(hObject, 'String');
Av = evalin('base', A);
dimensions = string(1:size(Av,1));
dimensions1 = ['', dimensions, 'time'];
dimensions2 = ['', dimensions];

handles.dimensions_Linear = size(Av,1);
guidata(hObject, handles)

set(handles.popDim1Plot_Linear, 'String', dimensions1);
set(handles.popDim1Plot_Linear, 'Value', 1);
set(handles.popDim2Plot_Linear, 'String', dimensions2);
set(handles.popDim2Plot_Linear, 'Value', 1);

handles.A = A;
set(handles.popA_Linear, 'Value', 1)
guidata(hObject, handles);



% --- Executes on button press in checkA_Linear.
function checkA_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkA_Linear

if get(hObject, 'Value')
    set(handles.textA_Linear, 'Enable', 'on')
    set(handles.popA_Linear, 'Enable', 'on')
else
    set(handles.textA_Linear, 'Enable', 'off')
    set(handles.popA_Linear, 'Enable', 'off')
end

% --- Executes during object creation, after setting all properties.
function textB_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textB_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function textB_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to textB_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

B = get(hObject, 'String');
handles.B = B;
guidata(hObject, handles);
set(handles.popB_Linear, 'Value', 1)

% --- Executes on button press in checkB_Linear.
function checkB_Linear_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    set(handles.textB_Linear, 'Enable', 'on')
    set(handles.popB_Linear, 'Enable', 'on')
else
    set(handles.textB_Linear, 'Enable', 'off')
    set(handles.popB_Linear, 'Enable', 'off')
end


% --- Executes during object creation, after setting all properties.
function textC_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textC_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function textC_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to textC_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

C = get(hObject, 'String');
handles.C = C;
guidata(hObject, handles);
set(handles.popC_Linear, 'Value', 1)


% --- Executes on button press in checkC_Linear.
function checkC_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkC_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkC_Linear
if get(hObject, 'Value')
    set(handles.textC_Linear, 'Enable', 'on')
    set(handles.popC_Linear, 'Enable', 'on')
else
    set(handles.textC_Linear, 'Enable', 'off')
    set(handles.popC_Linear, 'Enable', 'off')
end


% --- Executes during object creation, after setting all properties.
function textD_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textD_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function textD_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to textD_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

D = get(hObject, 'String');
handles.D = D;
guidata(hObject, handles);
set(handles.popD_Linear, 'Value', 1)

% --- Executes on button press in checkD_Linear.
function checkD_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkD_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkD_Linear
if get(hObject, 'Value')
    set(handles.textD_Linear, 'Enable', 'on')
    set(handles.popD_Linear, 'Enable', 'on')
else
    set(handles.textD_Linear, 'Enable', 'off')
    set(handles.popD_Linear, 'Enable', 'off')
end


% --- Executes during object creation, after setting all properties.
function textE_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textF_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function textE_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to textB_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

E = get(hObject, 'String');
handles.E = E;
guidata(hObject, handles);
set(handles.popE_Linear, 'Value', 1)

% --- Executes on button press in checkE_Linear.
function checkE_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkE_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkE_Linear
if get(hObject, 'Value')
    set(handles.textE_Linear, 'Enable', 'on')
    set(handles.popE_Linear, 'Enable', 'on')
else
    set(handles.textE_Linear, 'Enable', 'off')
    set(handles.popE_Linear, 'Enable', 'off')
end


% --- Executes during object creation, after setting all properties.
function textF_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textF_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function textF_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to textF_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

F = get(hObject, 'String');
handles.F = F;
guidata(hObject, handles);
set(handles.popF_Linear, 'Value', 1)

% --- Executes on button press in checkF_Linear.
function checkF_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkF_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkF_Linear
if get(hObject, 'Value')
    set(handles.textF_Linear, 'Enable', 'on')
    set(handles.popF_Linear, 'Enable', 'on')
else
    set(handles.textF_Linear, 'Enable', 'off')
    set(handles.popF_Linear, 'Enable', 'off')
end


% --- Executes during object creation, after setting all properties.
function listWS_Linear_CreateFcn(hObject, ~, handles)
% hObject    handle to listWS_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listWS_Linear.
function listWS_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to listWS_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textB_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textB_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Finish Linear -------------------------------------------------------


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function txtTstart_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTstart_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTstart_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtTstart_Linear as a double


% --- Executes during object creation, after setting all properties.
function txtTstart_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTstart_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTfinal_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTfinal_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTfinal_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtTfinal_Linear as a double


% --- Executes during object creation, after setting all properties.
function txtTfinal_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTfinal_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtTimeStep_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTimeStep_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimeStep_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtTimeStep_Linear as a double


% --- Executes during object creation, after setting all properties.
function txtTimeStep_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimeStep_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtZOrder_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtZOrder_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtZOrder_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtZOrder_Linear as a double
value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtZOrder_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZOrder_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownZOrder_Linear.
function pbDownZOrder_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownZOrder_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtZOrder_Linear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtZOrder_Linear, 'String', value);
end


% --- Executes on button press in pbUpZOrder_Linear.
function pbUpZOrder_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpZOrder_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtZOrder_Linear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtZOrder_Linear, 'String', value);


function txtTaylor_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTaylor_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTaylor_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtTaylor_Linear as a double
value = str2double(get(hObject, 'String'));
if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtTaylor_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTaylor_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbUpTaylor_Linear.
function pbUpTaylor_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpTaylor_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtTaylor_Linear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtTaylor_Linear, 'String', value);


% --- Executes on button press in pbDownTaylor_Linear.
function pbDownTaylor_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownTaylor_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtTaylor_Linear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtTaylor_Linear, 'String', value);
end


% --- Executes on selection change in popRT_Linear.
function popRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popRT_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popRT_Linear


% --- Executes during object creation, after setting all properties.
function popRT_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLA_Linear.
function popLA_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popLA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popLA_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLA_Linear

contents = cellstr(get(hObject, 'String'));
alg = contents{get(hObject,'Value')};

if strcmp(alg, 'adaptive')
    set(handles.OptionsError, 'visible', 'on')
    set(handles.txtOptionError, 'visible', 'on')
    set(handles.txtOptionErrorUp, 'visible', 'on')
    set(handles.txtOptionErrorDown, 'visible', 'on')
    set(handles.txtZOrder_Linear, 'visible', 'off')
    set(handles.pbUpZOrder_Linear, 'visible', 'off')
    set(handles.pbDownZOrder_Linear, 'visible', 'off')
    set(handles.txtTaylor_Linear, 'visible', 'off')
    set(handles.pbUpTaylor_Linear, 'visible', 'off')
    set(handles.pbDownTaylor_Linear, 'visible', 'off')
    set(handles.popRT_Linear, 'visible', 'off')
    set(handles.zotxt, 'visible', 'off')
    set(handles.tttxt, 'visible', 'off')
    set(handles.rttxt, 'visible', 'off')
    set(handles.pb_infoZO_Linear, 'visible', 'off')
    set(handles.pb_infoTT_Linear, 'visible', 'off')
    set(handles.pb_infoRT_Linear, 'visible', 'off')
    set(handles.txtTimeStep_Linear, 'visible', 'off')
    set(handles.tstxt, 'visible', 'off')
else
    set(handles.OptionsError, 'visible', 'off')
    set(handles.txtOptionError, 'visible', 'off')
    set(handles.txtOptionErrorUp, 'visible', 'off')
    set(handles.txtOptionErrorDown, 'visible', 'off')
    
    set(handles.txtZOrder_Linear, 'visible', 'on')
    set(handles.pbUpZOrder_Linear, 'visible', 'on')
    set(handles.pbDownZOrder_Linear, 'visible', 'on')
    set(handles.txtTaylor_Linear, 'visible', 'on')
    set(handles.pbUpTaylor_Linear, 'visible', 'on')
    set(handles.pbDownTaylor_Linear, 'visible', 'on')
    set(handles.popRT_Linear, 'visible', 'on')
    set(handles.zotxt, 'visible', 'on')
    set(handles.tttxt, 'visible', 'on')
    set(handles.rttxt, 'visible', 'on')
    set(handles.pb_infoZO_Linear, 'visible', 'on')
    set(handles.pb_infoTT_Linear, 'visible', 'on')
    set(handles.pb_infoRT_Linear, 'visible', 'on')
    set(handles.txtTimeStep_Linear, 'visible', 'on')
    set(handles.tstxt, 'visible', 'on')
    set(handles.popUTrans_Linear, 'visible', 'on')
    set(handles.txtu_Linear, 'visible', 'on')
    set(handles.txtUTrans_Linear, 'visible', 'on')
    
    
end


% --- Executes during object creation, after setting all properties.
function popLA_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbR0_Linear.
function pbR0_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbR0_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [R0, handles_initial_set] = input_set(handles.initial_set_linear,handles.initial_set_linear_type);
    handles.R0 = R0;
    handles.initial_set_linear = handles_initial_set;
    handles.initial_set_linear_type = R0.type;
    guidata(hObject, handles)
catch
end


% --- Executes on button press in pbU_Linear.
function pbU_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbU_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%run('input_set.m')
%input_set_h = findobj('Tag','input_set');
%input_set_handles = guidata(input_set_h);

try
    [U, handles_input_set] = input_set(handles.input_set_linear,handles.input_set_linear_type);
    handles.U = U;
    handles.input_set_linear = handles_input_set;
    handles.input_set_linear_type = U.type;
    guidata(hObject, handles)
catch
end


% --- Executes on selection change in popUTrans_Linear.
function popUTrans_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popUTrans_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUTrans_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUTrans_Linear
contents = cellstr(get(hObject,'String'));
u = contents{get(hObject,'Value')};
u_value = evalin('base', u);
set(handles.txtUTrans_Linear,'String', mat2str(u_value));
handles.uTrans = mat2str(u_value);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popUTrans_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUTrans_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in rbSim_Linear.
function rbSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to rbSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbSim_Linear
if get(hObject,'Value')
    set(handles.rbRSim_Linear, 'Value', 0)
    set(handles.rbSimRRT_Linear, 'Value', 0)
    
    set(handles.checkSimulationPlot_Linear, 'Enable', 'on')
    set(handles.popColorSimulation_Linear, 'Enable', 'on')
    set(handles.popLineSimulation_Linear, 'Enable', 'on')
    set(handles.text43, 'Enable', 'on')
    set(handles.text44, 'Enable', 'on')

    set(handles.txtX0_Sim_Linear, 'Enable', 'on')
    set(handles.txtU_Sim_Linear, 'Enable', 'on')
    set(handles.popX0_Sim_Linear, 'Enable', 'on')
    set(handles.popU_Sim_Linear, 'Enable', 'on')
    
    set(handles.txtNoP_RSim_Linear, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Linear, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Linear, 'Enable', 'off')
    set(handles.txtFV_RSim_Linear, 'Enable', 'off')
    set(handles.txtFIV_RSim_Linear, 'Enable', 'off')
    set(handles.txtNCI_RSim_Linear, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Linear, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Linear, 'Enable', 'off')
    set(handles.txtNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.rbYesEPS_SimRRT_Linear, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Linear, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Linear, 'Enable', 'off')
    
else
    
    set(handles.checkSimulationPlot_Linear, 'Enable', 'off')
    set(handles.popColorSimulation_Linear, 'Enable', 'off')
    set(handles.popLineSimulation_Linear, 'Enable', 'off')
    set(handles.text43, 'Enable', 'off')
    set(handles.text44, 'Enable', 'off')
    
    set(handles.txtX0_Sim_Linear, 'Enable', 'off')
    set(handles.txtU_Sim_Linear, 'Enable', 'off')
    set(handles.popX0_Sim_Linear, 'Enable', 'off')
    set(handles.popU_Sim_Linear, 'Enable', 'off')
    
end


function txtX0_Sim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtX0_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtX0_Sim_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtX0_Sim_Linear as a double
x0 = get(hObject, 'String');
handles.x0 = x0;
guidata(hObject, handles);
set(handles.popX0_Sim_Linear, 'Value', 1)


% --- Executes during object creation, after setting all properties.
function txtX0_Sim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtX0_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtU_Sim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtU_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtU_Sim_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtU_Sim_Linear as a double
u = get(hObject, 'String');
%set(handles.txtU_Sim_Linear,'String', u);
handles.uSim = u;
guidata(hObject, handles);
set(handles.popU_Sim_Linear, 'Value', 1)


% --- Executes during object creation, after setting all properties.
function txtU_Sim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtU_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbRSim_Linear.
function rbRSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to rbRSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbRSim_Linear
if get(hObject,'Value')
    set(handles.rbSim_Linear, 'Value', 0)
    set(handles.rbSimRRT_Linear, 'Value', 0)

    set(handles.checkSimulationPlot_Linear, 'Enable', 'on')
    set(handles.popColorSimulation_Linear, 'Enable', 'on')
    set(handles.popLineSimulation_Linear, 'Enable', 'on')
    set(handles.text43, 'Enable', 'on')
    set(handles.text44, 'Enable', 'on')

    set(handles.txtX0_Sim_Linear, 'Enable', 'off')
    set(handles.txtU_Sim_Linear, 'Enable', 'off')
    set(handles.popX0_Sim_Linear, 'Enable', 'off')
    set(handles.popU_Sim_Linear, 'Enable', 'off')
    
    set(handles.txtNoP_RSim_Linear, 'Enable', 'on')
    set(handles.pbUpNoP_RSim_Linear, 'Enable', 'on')
    set(handles.pbDownNoP_RSim_Linear, 'Enable', 'on')
    set(handles.txtFV_RSim_Linear, 'Enable', 'on')
    set(handles.txtFIV_RSim_Linear, 'Enable', 'on')
    set(handles.txtNCI_RSim_Linear, 'Enable', 'on')
    set(handles.pbUpNCI_RSim_Linear, 'Enable', 'on')
    set(handles.pbDownNCI_RSim_Linear, 'Enable', 'on')
    
    set(handles.txtNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.rbYesEPS_SimRRT_Linear, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Linear, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Linear, 'Enable', 'off')
    
else
    
    set(handles.checkSimulationPlot_Linear, 'Enable', 'off')
    set(handles.popColorSimulation_Linear, 'Enable', 'off')
    set(handles.popLineSimulation_Linear, 'Enable', 'off')
    set(handles.text43, 'Enable', 'off')
    set(handles.text44, 'Enable', 'off')
    
    set(handles.txtNoP_RSim_Linear, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Linear, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Linear, 'Enable', 'off')
    set(handles.txtFV_RSim_Linear, 'Enable', 'off')
    set(handles.txtFIV_RSim_Linear, 'Enable', 'off')
    set(handles.txtNCI_RSim_Linear, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Linear, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Linear, 'Enable', 'off')
    
end


function txtNoP_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtNoP_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtNoP_RSim_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtNoP_RSim_Linear as a double
value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtNoP_RSim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNoP_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownNoP_RSim_Linear.
function pbDownNoP_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNoP_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_RSim_Linear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNoP_RSim_Linear, 'String', value);
end


% --- Executes on button press in pbUpNoP_RSim_Linear.
function pbUpNoP_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNoP_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_RSim_Linear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNoP_RSim_Linear, 'String', value);


function txtFV_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtFV_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFV_RSim_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtFV_RSim_Linear as a double
value = str2double(get(hObject, 'String'));
if value > 1
    value = 1;
elseif value < 0
    value = 0;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtFV_RSim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFV_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

function txtFIV_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtFIV_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFIV_RSim_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtFIV_RSim_Linear as a double
value = str2double(get(hObject, 'String'));
if value > 1
    value = 1;
elseif value < 0
    value = 0;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtFIV_RSim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFIV_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txtNCI_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtNCI_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtNCI_RSim_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtNCI_RSim_Linear as a double
value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtNCI_RSim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNCI_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbSimRRT_Linear.
function rbSimRRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to rbSimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbSimRRT_Linear

if get(hObject,'Value')
    set(handles.rbRSim_Linear, 'Value', 0)
    set(handles.rbSim_Linear, 'Value', 0)

    set(handles.checkSimulationPlot_Linear, 'Enable', 'on')
    set(handles.popColorSimulation_Linear, 'Enable', 'on')
    set(handles.popLineSimulation_Linear, 'Enable', 'on')
    set(handles.text43, 'Enable', 'on')
    set(handles.text44, 'Enable', 'on')

    set(handles.txtX0_Sim_Linear, 'Enable', 'off')
    set(handles.txtU_Sim_Linear, 'Enable', 'off')
    set(handles.popX0_Sim_Linear, 'Enable', 'off')
    set(handles.popU_Sim_Linear, 'Enable', 'off')
    set(handles.txtNoP_RSim_Linear, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Linear, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Linear, 'Enable', 'off')
    set(handles.txtFV_RSim_Linear, 'Enable', 'off')
    set(handles.txtFIV_RSim_Linear, 'Enable', 'off')
    set(handles.txtNCI_RSim_Linear, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Linear, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Linear, 'Enable', 'off')
    
    set(handles.txtNoP_SimRRT_Linear, 'Enable', 'on')
    set(handles.pbUpNoP_SimRRT_Linear, 'Enable', 'on')
    set(handles.pbDownNoP_SimRRT_Linear, 'Enable', 'on')
    set(handles.rbNoEPS_SimRRT_Linear, 'Enable', 'on')
    set(handles.txtSF_SimRRT_Linear, 'Enable', 'on')
    
else
    
    set(handles.checkSimulationPlot_Linear, 'Enable', 'off')
    set(handles.popColorSimulation_Linear, 'Enable', 'off')
    set(handles.popLineSimulation_Linear, 'Enable', 'off')
    set(handles.text43, 'Enable', 'off')
    set(handles.text44, 'Enable', 'off')
   
    set(handles.txtNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Linear, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Linear, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Linear, 'Enable', 'off')
    
end


function txtNoP_SimRRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtNoP_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtNoP_SimRRT_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtNoP_SimRRT_Linear as a double
value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtNoP_SimRRT_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNoP_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownNoP_SimRRT_Linear.
function pbDownNoP_SimRRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNoP_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_SimRRT_Linear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNoP_SimRRT_Linear, 'String', value);
end

% --- Executes on button press in pbUpNoP_SimRRT_Linear.
function pbUpNoP_SimRRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNoP_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_SimRRT_Linear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNoP_SimRRT_Linear, 'String', value);


% --- Executes on button press in pbDownNCI_RSim_Linear.
function pbDownNCI_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNCI_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNCI_RSim_Linear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNCI_RSim_Linear, 'String', value);
end

% --- Executes on button press in pbUpNCI_RSim_Linear.
function pbUpNCI_RSim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNCI_RSim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNCI_RSim_Linear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNCI_RSim_Linear, 'String', value);


% --- Executes on button press in rbYesEPS_SimRRT_Linear.
function rbYesEPS_SimRRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to rbYesEPS_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbYesEPS_SimRRT_Linear
set(hObject, 'Enable', 'off')
set(handles.rbNoEPS_SimRRT_Linear, 'Enable', 'on')
set(handles.rbNoEPS_SimRRT_Linear, 'Value', 0)


% --- Executes on button press in rbNoEPS_SimRRT_Linear.
function rbNoEPS_SimRRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to rbNoEPS_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbNoEPS_SimRRT_Linear
set(hObject, 'Enable', 'off')
set(handles.rbYesEPS_SimRRT_Linear, 'Enable', 'on')
set(handles.rbYesEPS_SimRRT_Linear, 'Value', 0)


function txtSF_SimRRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtSF_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSF_SimRRT_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtSF_SimRRT_Linear as a double

value = str2double(get(hObject, 'String'));
if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtSF_SimRRT_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSF_SimRRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popX0_Sim_Linear.
function popX0_Sim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popX0_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popX0_Sim_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popX0_Sim_Linear
contents = cellstr(get(hObject,'String'));
x0 = contents{get(hObject,'Value')};
x0_value = evalin('base', x0);
set(handles.txtX0_Sim_Linear,'String', mat2str(x0_value));
handles.x0 = x0;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popX0_Sim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popX0_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popU_Sim_Linear.
function popU_Sim_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popU_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popU_Sim_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popU_Sim_Linear
contents = cellstr(get(hObject,'String'));
u = contents{get(hObject,'Value')};
u_value = evalin('base', u);
set(handles.txtU_Sim_Linear,'String', mat2str(u_value));
handles.uSim = u;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popU_Sim_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popU_Sim_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtUTrans_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to txtUTrans_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtUTrans_Linear as text
%        str2double(get(hObject,'String')) returns contents of txtUTrans_Linear as a double
u = get(hObject, 'String');
handles.uTrans = u;
set(handles.popUTrans_Linear, 'Value', 1);
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function txtUTrans_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtUTrans_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popDim1Plot_Linear.
function popDim1Plot_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popDim1Plot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popDim1Plot_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popDim1Plot_Linear


% --- Executes during object creation, after setting all properties.
function popDim1Plot_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDim1Plot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popDim2Plot_Linear.
function popDim2Plot_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popDim2Plot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popDim2Plot_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popDim2Plot_Linear


% --- Executes during object creation, after setting all properties.
function popDim2Plot_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDim2Plot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushAddPlot_Linear.
function pushAddPlot_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pushAddPlot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listPlots_Linear, 'Value', 1);

contents_dim1 = cellstr(get(handles.popDim1Plot_Linear, 'String'));
dim1 = contents_dim1{get(handles.popDim1Plot_Linear, 'Value')};
dim1 = dim1(~isspace(dim1));

contents_dim2 = cellstr(get(handles.popDim2Plot_Linear, 'String'));
dim2 = contents_dim2{get(handles.popDim2Plot_Linear, 'Value')};
dim2 = dim2(~isspace(dim2));

if isempty(dim1) || isempty(dim2)
    return
end

if dim1 == 'time'
    dim = sprintf('[%s, %i]', dim1, str2double(dim2));
else
    dim = sprintf('[%i, %i]', str2double(dim1), str2double(dim2));
end
handles.linear_plots = [handles.linear_plots; dim];
guidata(hObject, handles);

set(handles.listPlots_Linear, 'String', handles.linear_plots);


% --- Executes on button press in pushDeletePlot_Linear.
function pushDeletePlot_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pushDeletePlot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.listPlots_Linear,'String'));
value = get(handles.listPlots_Linear, 'Value');

if isempty(contents)
    return
end

contents(value,:) = [];

handles.linear_plots = contents';
guidata(hObject, handles)
set(handles.listPlots_Linear, 'Value', 1);
set(handles.listPlots_Linear, 'String', contents);


% --- Executes on selection change in listPlots_Linear.
function listPlots_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to listPlots_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listPlots_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listPlots_Linear


% --- Executes during object creation, after setting all properties.
function listPlots_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listPlots_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkReachPlot_Linear.
function checkReachPlot_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkReachPlot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkReachPlot_Linear


% --- Executes on button press in checkInitialPlot_Linear.
function checkInitialPlot_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkInitialPlot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkInitialPlot_Linear


% --- Executes on button press in checkSimulationPlot_Linear.
function checkSimulationPlot_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to checkSimulationPlot_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkSimulationPlot_Linear


% --- Executes on selection change in popColorSimulation_Linear.
function popColorSimulation_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popColorSimulation_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorSimulation_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorSimulation_Linear


% --- Executes during object creation, after setting all properties.
function popColorSimulation_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorSimulation_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLineSimulation_Linear.
function popLineSimulation_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popLineSimulation_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popLineSimulation_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLineSimulation_Linear


% --- Executes during object creation, after setting all properties.
function popLineSimulation_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLineSimulation_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popColorInitial_Linear.
function popColorInitial_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popColorInitial_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorInitial_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorInitial_Linear


% --- Executes during object creation, after setting all properties.
function popColorInitial_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorInitial_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popEdgeColorInitial_Linear.
function popEdgeColorInitial_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popEdgeColorInitial_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popEdgeColorInitial_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popEdgeColorInitial_Linear


% --- Executes during object creation, after setting all properties.
function popEdgeColorInitial_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popEdgeColorInitial_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popColorReach_Linear.
function popColorReach_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popColorReach_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorReach_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorReach_Linear


% --- Executes during object creation, after setting all properties.
function popColorReach_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorReach_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popEdgeColorReach_Linear.
function popEdgeColorReach_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popEdgeColorReach_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popEdgeColorReach_Linear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popEdgeColorReach_Linear


% --- Executes during object creation, after setting all properties.
function popEdgeColorReach_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popEdgeColorReach_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Nonlinear System ----------------------------------------------------

% --- Callback Functions - Nonlinear --------------------------------------


function txtTstart_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTstart_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTstart_Nonlinear as text
%        str2double(get(hObject,'String')) returns contents of txtTstart_Nonlinear as a double


% --- Executes during object creation, after setting all properties.
function txtTstart_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTstart_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtTfinal_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTfinal_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTfinal_Nonlinear as text
%        str2double(get(hObject,'String')) returns contents of txtTfinal_Nonlinear as a double


% --- Executes during object creation, after setting all properties.
function txtTfinal_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTfinal_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTimeStep_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTimeStep_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimeStep_Nonlinear as text
%        str2double(get(hObject,'String')) returns contents of txtTimeStep_Nonlinear as a double


% --- Executes during object creation, after setting all properties.
function txtTimeStep_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimeStep_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton51.
function pushbutton51_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu31.
function popupmenu31_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu31 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu31


% --- Executes during object creation, after setting all properties.
function popupmenu31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu32.
function popupmenu32_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu32 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu32


% --- Executes during object creation, after setting all properties.
function popupmenu32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton52.
function pushbutton52_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu33.
function popupmenu33_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu33 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu33


% --- Executes during object creation, after setting all properties.
function popupmenu33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu34.
function popupmenu34_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu34 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu34


% --- Executes during object creation, after setting all properties.
function popupmenu34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton18.
function radiobutton18_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton18


% --- Executes on button press in radiobutton19.
function radiobutton19_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton19


% --- Executes on selection change in popDim1Plot_Nonlinear.
function popDim1Plot_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popDim1Plot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popDim1Plot_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popDim1Plot_Nonlinear


% --- Executes during object creation, after setting all properties.
function popDim1Plot_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDim1Plot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popDim2Plot_Nonlinear.
function popDim2Plot_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popDim2Plot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popDim2Plot_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popDim2Plot_Nonlinear


% --- Executes during object creation, after setting all properties.
function popDim2Plot_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDim2Plot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushAddPlot_Nonlinear.
function pushAddPlot_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pushAddPlot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.listPlots_Nonlinear, 'Value', 1);

contents_dim1 = cellstr(get(handles.popDim1Plot_Nonlinear, 'String'));
dim1 = contents_dim1{get(handles.popDim1Plot_Nonlinear, 'Value')};
dim1 = dim1(~isspace(dim1));

contents_dim2 = cellstr(get(handles.popDim2Plot_Nonlinear, 'String'));
dim2 = contents_dim2{get(handles.popDim2Plot_Nonlinear, 'Value')};
dim2 = dim2(~isspace(dim2));

if isempty(dim1) || isempty(dim2)
    return
end
if dim1 == 'time'
    dim = sprintf('[%s, %i]', dim1, str2double(dim2));
else
    dim = sprintf('[%i, %i]', str2double(dim1), str2double(dim2));
end
handles.nonlinear_plots = [handles.nonlinear_plots, dim];
guidata(hObject, handles);

set(handles.listPlots_Nonlinear, 'String', handles.nonlinear_plots);

% --- Executes on selection change in listPlots_Nonlinear.
function listPlots_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to listPlots_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listPlots_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listPlots_Nonlinear


% --- Executes during object creation, after setting all properties.
function listPlots_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listPlots_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkReachPlot_Nonlinear.
function checkReachPlot_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to checkReachPlot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkReachPlot_Nonlinear


% --- Executes on button press in checkInitialPlot_Nonlinear.
function checkInitialPlot_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to checkInitialPlot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkInitialPlot_Nonlinear


% --- Executes on button press in checkSimulationPlot_Nonlinear.
function checkSimulationPlot_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to checkSimulationPlot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkSimulationPlot_Nonlinear


% --- Executes on button press in pbLoadSModel_Nonlinear.
function pbLoadSModel_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadSModel_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[spaceEx_file, spaceEx_path] = uigetfile('*.xml');
try
    handles.spaceExFile_Nonlinear = fullfile(spaceEx_path, spaceEx_file);
    spaceex2cora(handles.spaceExFile_Nonlinear);
catch
    return
end

file = strsplit(handles.spaceExFile_Nonlinear, filesep);
file = file{end};
file = strsplit(file, '.');
file = file{1};
nonlinSys =  eval(strcat(file, '()'));

if strcmp(class(nonlinSys),'nonlinearSys')
    dimensions = string(1:nonlinSys.dim);
    dimensions1 = ['', dimensions, 'time'];
    dimensions2 = ['', dimensions];

    handles.dimensions_Nonlinear = nonlinSys.dim;
    handles.nrOfInputs_Nonlinear = nonlinSys.nrOfInputs;

    set(handles.popDim1Plot_Nonlinear, 'String', dimensions1);
    set(handles.popDim1Plot_Nonlinear, 'Value', 1);
    set(handles.popDim2Plot_Nonlinear, 'String', dimensions2);
    set(handles.popDim2Plot_Nonlinear, 'Value', 1);
    set(hObject, 'string', spaceEx_file)
    guidata(hObject, handles) %update handles
else
uiwait(msgbox('The loaded SpaceEx file is not nonlinear', 'Error', 'error', 'modal'))
return
end


% --- Executes on button press in rbLoadSModel_Nonlinear.
function rbLoadSModel_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbLoadSModel_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dimensions_Nonlinear = [];
handles.nrOfInputs_Nonlinear = [];
guidata(hObject, handles)

set(handles.popDim1Plot_Nonlinear, 'String', ' ');
set(handles.popDim1Plot_Nonlinear, 'Value', 1);
set(handles.popDim2Plot_Nonlinear, 'String', ' ');
set(handles.popDim2Plot_Nonlinear, 'Value', 1);

set(hObject, 'Enable', 'off')
set(handles.pbLoadSModel_Nonlinear, 'Enable', 'on')
set(handles.rbEnterDynamicEquation_Nonlinear, 'Enable', 'on')
set(handles.rbEnterDynamicEquation_Nonlinear, 'Value', 0)
set(handles.txtEnterDynamicEquation_Nonlinear, 'Enable', 'off')
%set(handles.txtEnterDynamicEquation_Nonlinear, 'String', ' ')
set(handles.rbLoadDynamicEquation_Nonlinear, 'Enable', 'on')
set(handles.rbLoadDynamicEquation_Nonlinear, 'Value', 0)
set(handles.pbLoadDynamicEquation_Nonlinear, 'Enable', 'off')
set(handles.pbLoadDynamicEquation_Nonlinear, 'String', 'Load File')
set(handles.textFunction_Nonlinear, 'Enable', 'off')
set(handles.textEnd_Nonlinear, 'Enable', 'off')



% --- Executes on button press in radiobutton23.
function radiobutton23_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton23



function edit55_Callback(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit55 as text
%        str2double(get(hObject,'String')) returns contents of edit55 as a double


% --- Executes during object creation, after setting all properties.
function edit55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox43.
function checkbox43_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox43



function edit56_Callback(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit56 as text
%        str2double(get(hObject,'String')) returns contents of edit56 as a double


% --- Executes during object creation, after setting all properties.
function edit56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox44.
function checkbox44_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox44



function edit57_Callback(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit57 as text
%        str2double(get(hObject,'String')) returns contents of edit57 as a double


% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox45.
function checkbox45_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox45



function edit58_Callback(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit58 as text
%        str2double(get(hObject,'String')) returns contents of edit58 as a double


% --- Executes during object creation, after setting all properties.
function edit58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox46.
function checkbox46_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox46



function edit59_Callback(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit59 as text
%        str2double(get(hObject,'String')) returns contents of edit59 as a double


% --- Executes during object creation, after setting all properties.
function edit59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox47.
function checkbox47_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox47



function edit60_Callback(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit60 as text
%        str2double(get(hObject,'String')) returns contents of edit60 as a double


% --- Executes during object creation, after setting all properties.
function edit60_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox48.
function checkbox48_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox48


% --- Executes on selection change in listbox14.
function listbox14_Callback(hObject, eventdata, handles)
% hObject    handle to listbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox14 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox14


% --- Executes during object creation, after setting all properties.
function listbox14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton24.
function radiobutton24_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton24



function edit61_Callback(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit61 as text
%        str2double(get(hObject,'String')) returns contents of edit61 as a double


% --- Executes during object creation, after setting all properties.
function edit61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit62_Callback(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit62 as text
%        str2double(get(hObject,'String')) returns contents of edit62 as a double


% --- Executes during object creation, after setting all properties.
function edit62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton25.
function radiobutton25_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton25



function edit63_Callback(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit63 as text
%        str2double(get(hObject,'String')) returns contents of edit63 as a double


% --- Executes during object creation, after setting all properties.
function edit63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton60.
function pushbutton60_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton61.
function pushbutton61_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit64_Callback(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit64 as text
%        str2double(get(hObject,'String')) returns contents of edit64 as a double


% --- Executes during object creation, after setting all properties.
function edit64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit65_Callback(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit65 as text
%        str2double(get(hObject,'String')) returns contents of edit65 as a double


% --- Executes during object creation, after setting all properties.
function edit65_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit66_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double


% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton26.
function radiobutton26_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton26



function edit67_Callback(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit67 as text
%        str2double(get(hObject,'String')) returns contents of edit67 as a double


% --- Executes during object creation, after setting all properties.
function edit67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton62.
function pushbutton62_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton63.
function pushbutton63_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton64.
function pushbutton64_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton65.
function pushbutton65_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton27.
function radiobutton27_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton27


% --- Executes on button press in radiobutton28.
function radiobutton28_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton28


function edit68_Callback(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit68 as text
%        str2double(get(hObject,'String')) returns contents of edit68 as a double


% --- Executes during object creation, after setting all properties.
function edit68_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu67.
function popupmenu67_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu67 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu67


% --- Executes during object creation, after setting all properties.
function popupmenu67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu68.~on popupmenu68_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu68 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu68


% --- Executes during object creation, after setting all properties.
function popupmenu68_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popColorSimulation_Nonlinear.
function popColorSimulation_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popColorSimulation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorSimulation_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorSimulation_Nonlinear


% --- Executes during object creation, after setting all properties.
function popColorSimulation_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorSimulation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLineSimulation_Nonlinear.
function popLineSimulation_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popLineSimulation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popLineSimulation_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLineSimulation_Nonlinear


% --- Executes during object creation, after setting all properties.
function popLineSimulation_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLineSimulation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popColorInitial_Nonlinear.
function popColorInitial_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popColorInitial_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorInitial_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorInitial_Nonlinear


% --- Executes during object creation, after setting all properties.
function popColorInitial_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorInitial_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popEdgeColorInitial_Nonlinear.
function popEdgeColorInitial_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popEdgeColorInitial_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popEdgeColorInitial_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popEdgeColorInitial_Nonlinear


% --- Executes during object creation, after setting all properties.
function popEdgeColorInitial_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popEdgeColorInitial_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popColorReach_Nonlinear.
function popColorReach_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popColorReach_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function popColorReach_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorReach_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popEdgeColorReach_Nonlinear.
function popEdgeColorReach_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popEdgeColorReach_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function popEdgeColorReach_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popEdgeColorReach_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtZOrder_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtZOrder_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtZOrder_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZOrder_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownZOrder_Nonlinear.
function pbDownZOrder_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownZOrder_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtZOrder_Nonlinear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtZOrder_Nonlinear, 'String', value);
end


% --- Executes on button press in pbUpZOrder_Nonlinear.
function pbUpZOrder_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpZOrder_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtZOrder_Nonlinear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtZOrder_Nonlinear, 'String', value);


function txtTaylor_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtTaylor_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));
if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtTaylor_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTaylor_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbUpTaylor_Nonlinear.
function pbUpTaylor_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpTaylor_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtTaylor_Nonlinear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtTaylor_Nonlinear, 'String', value);


% --- Executes on button press in pbDownTaylor_Nonlinear.
function pbDownTaylor_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownTaylor_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtTaylor_Nonlinear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtTaylor_Nonlinear, 'String', value);
end


% --- Executes on selection change in popRT_Nonlinear.
function popRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popRT_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popRT_Nonlinear


% --- Executes during object creation, after setting all properties.
function popRT_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLA_Nonlinear.
function popLA_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popLA_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popLA_Nonlinear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLA_Nonlinear


% --- Executes during object creation, after setting all properties.
function popLA_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLA_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbR0_Nonlinear.
function pbR0_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbR0_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [R0, handles_initial_set] = input_set(handles.initial_set_nonlinear,handles.initial_set_nonlinear_type);
    handles.R0_Nonlinear = R0;
    handles.initial_set_nonlinear = handles_initial_set;
    handles.initial_set_nonlinear_type = R0.type;
    guidata(hObject, handles)
catch
end


% --- Executes on button press in pbU_Nonlinear.
function pbU_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbU_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [U, handles_input_set] = input_set(handles.input_set_nonlinear,handles.input_set_nonlinear_type);
    handles.U_Nonlinear = U;
    handles.input_set_nonlinear = handles_input_set;
    handles.input_set_nonlinear_type = U.type;
    guidata(hObject, handles)
catch
end


% --- Executes on selection change in popUTrans_Nonlinear.
function popUTrans_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popUTrans_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
u = contents{get(hObject,'Value')};
u_value = evalin('base', u);
set(handles.txtUTrans_Nonlinear,'String', mat2str(u_value));
handles.uTrans_Nonlinear = mat2str(u_value);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popUTrans_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUTrans_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtUTrans_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtUTrans_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
u = get(hObject, 'String');
handles.uTrans_Nonlinear = u;
guidata(hObject, handles);
set(handles.popUTrans_Nonlinear, 'Value', 1)


% --- Executes during object creation, after setting all properties.
function txtUTrans_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtUTrans_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in rbSim_Nonlinear.
function rbSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.rbRSim_Nonlinear, 'Value', 0)
    set(handles.rbSimRRT_Nonlinear, 'Value', 0)

    set(handles.checkSimulationPlot_Nonlinear, 'Enable', 'on')
    set(handles.popColorSimulation_Nonlinear, 'Enable', 'on')
    set(handles.popLineSimulation_Nonlinear, 'Enable', 'on')
    set(handles.text103, 'Enable', 'on')
    set(handles.text104, 'Enable', 'on')

    set(handles.txtX0_Sim_Nonlinear, 'Enable', 'on')
    set(handles.txtU_Sim_Nonlinear, 'Enable', 'on')
    set(handles.popX0_Sim_Nonlinear, 'Enable', 'on')
    set(handles.popU_Sim_Nonlinear, 'Enable', 'on')
    
    set(handles.txtNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtFV_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtFIV_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtNCI_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.rbYesEPS_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Nonlinear, 'Enable', 'off')
    
else
    
    set(handles.checkSimulationPlot_Nonlinear, 'Enable', 'off')
    set(handles.popColorSimulation_Nonlinear, 'Enable', 'off')
    set(handles.popLineSimulation_Nonlinear, 'Enable', 'off')
    set(handles.text103, 'Enable', 'off')
    set(handles.text104, 'Enable', 'off')

    set(handles.txtX0_Sim_Nonlinear, 'Enable', 'off')
    set(handles.txtU_Sim_Nonlinear, 'Enable', 'off')
    set(handles.popX0_Sim_Nonlinear, 'Enable', 'off')
    set(handles.popU_Sim_Nonlinear, 'Enable', 'off')
    
end



function txtX0_Sim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtX0_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x0 = get(hObject, 'String');
handles.x0_nonlinear = x0;
guidata(hObject, handles);
set(handles.popX0_Sim_Nonlinear, 'Value', 1)

% --- Executes during object creation, after setting all properties.
function txtX0_Sim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtX0_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtU_Sim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtU_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
u = get(hObject, 'String');
handles.u_nonlinear = u;
guidata(hObject, handles);
set(handles.popU_Sim_Nonlinear, 'Value', 1)

% --- Executes during object creation, after setting all properties.
function txtU_Sim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtU_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbRSim_Nonlinear.
function rbRSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbRSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.rbSim_Nonlinear, 'Value', 0)
    set(handles.rbSimRRT_Nonlinear, 'Value', 0)
    
    set(handles.checkSimulationPlot_Nonlinear, 'Enable', 'on')
    set(handles.popColorSimulation_Nonlinear, 'Enable', 'on')
    set(handles.popLineSimulation_Nonlinear, 'Enable', 'on')
    set(handles.text103, 'Enable', 'on')
    set(handles.text104, 'Enable', 'on')

    set(handles.txtX0_Sim_Nonlinear, 'Enable', 'off')
    set(handles.txtU_Sim_Nonlinear, 'Enable', 'off')
    set(handles.popX0_Sim_Nonlinear, 'Enable', 'off')
    set(handles.popU_Sim_Nonlinear, 'Enable', 'off')
    
    set(handles.txtNoP_RSim_Nonlinear, 'Enable', 'on')
    set(handles.pbUpNoP_RSim_Nonlinear, 'Enable', 'on')
    set(handles.pbDownNoP_RSim_Nonlinear, 'Enable', 'on')
    set(handles.txtFV_RSim_Nonlinear, 'Enable', 'on')
    set(handles.txtFIV_RSim_Nonlinear, 'Enable', 'on')
    set(handles.txtNCI_RSim_Nonlinear, 'Enable', 'on')
    set(handles.pbUpNCI_RSim_Nonlinear, 'Enable', 'on')
    set(handles.pbDownNCI_RSim_Nonlinear, 'Enable', 'on')
    
    set(handles.txtNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.rbYesEPS_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Nonlinear, 'Enable', 'off')
    
else
    
    set(handles.checkSimulationPlot_Nonlinear, 'Enable', 'off')
    set(handles.popColorSimulation_Nonlinear, 'Enable', 'off')
    set(handles.popLineSimulation_Nonlinear, 'Enable', 'off')
    set(handles.text103, 'Enable', 'off')
    set(handles.text104, 'Enable', 'off')

    set(handles.txtNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtFV_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtFIV_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtNCI_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Nonlinear, 'Enable', 'off')
    
end



function txtNoP_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtNoP_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtNoP_RSim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNoP_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownNoP_RSim_Nonlinear.
function pbDownNoP_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNoP_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_RSim_Nonlinear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNoP_RSim_Nonlinear, 'String', value);
end

% --- Executes on button press in pbUpNoP_RSim_Nonlinear.
function pbUpNoP_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNoP_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_RSim_Nonlinear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNoP_RSim_Nonlinear, 'String', value);


function txtFV_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtFV_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));
if value > 1
    value = 1;
elseif value < 0
    value = 0;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtFV_RSim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFV_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFIV_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtFIV_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));
if value > 1
    value = 1;
elseif value < 0
    value = 0;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtFIV_RSim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFIV_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtNCI_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtNCI_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtNCI_RSim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNCI_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbSimRRT_Nonlinear.
function rbSimRRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbSimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.rbRSim_Nonlinear, 'Value', 0)
    set(handles.rbSim_Nonlinear, 'Value', 0)
    
    set(handles.checkSimulationPlot_Nonlinear, 'Enable', 'on')
    set(handles.popColorSimulation_Nonlinear, 'Enable', 'on')
    set(handles.popLineSimulation_Nonlinear, 'Enable', 'on')
    set(handles.text103, 'Enable', 'on')
    set(handles.text104, 'Enable', 'on')

    set(handles.txtX0_Sim_Nonlinear, 'Enable', 'off')
    set(handles.txtU_Sim_Nonlinear, 'Enable', 'off')
    set(handles.popX0_Sim_Nonlinear, 'Enable', 'off')
    set(handles.popU_Sim_Nonlinear, 'Enable', 'off')
    set(handles.txtNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtFV_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtFIV_RSim_Nonlinear, 'Enable', 'off')
    set(handles.txtNCI_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Nonlinear, 'Enable', 'off')
    
    set(handles.txtNoP_SimRRT_Nonlinear, 'Enable', 'on')
    set(handles.pbUpNoP_SimRRT_Nonlinear, 'Enable', 'on')
    set(handles.pbDownNoP_SimRRT_Nonlinear, 'Enable', 'on')
    set(handles.rbNoEPS_SimRRT_Nonlinear, 'Enable', 'on')
    set(handles.txtSF_SimRRT_Nonlinear, 'Enable', 'on')
    
else
   
    set(handles.checkSimulationPlot_Nonlinear, 'Enable', 'off')
    set(handles.popColorSimulation_Nonlinear, 'Enable', 'off')
    set(handles.popLineSimulation_Nonlinear, 'Enable', 'off')
    set(handles.text103, 'Enable', 'off')
    set(handles.text104, 'Enable', 'off')

    set(handles.txtNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Nonlinear, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Nonlinear, 'Enable', 'off')
    
end



function txtNoP_SimRRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtNoP_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);



% --- Executes during object creation, after setting all properties.
function txtNoP_SimRRT_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNoP_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownNoP_SimRRT_Nonlinear.
function pbDownNoP_SimRRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNoP_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_SimRRT_Nonlinear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNoP_SimRRT_Nonlinear, 'String', value);
end



function pbUpNoP_SimRRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNoP_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = get(handles.txtNoP_SimRRT_Nonlinear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNoP_SimRRT_Nonlinear, 'String', value);


% --- Executes on button press in pbDownNCI_RSim_Nonlinear.
function pbDownNCI_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNCI_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNCI_RSim_Nonlinear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNCI_RSim_Nonlinear, 'String', value);
end

% --- Executes on button press in pbUpNCI_RSim_Nonlinear.
function pbUpNCI_RSim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNCI_RSim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNCI_RSim_Nonlinear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNCI_RSim_Nonlinear, 'String', value);


% --- Executes on button press in rbYesEPS_SimRRT_Nonlinear.
function rbYesEPS_SimRRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbYesEPS_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'off')
set(handles.rbNoEPS_SimRRT_Nonlinear, 'Enable', 'on')
set(handles.rbNoEPS_SimRRT_Nonlinear, 'Value', 0)

% --- Executes on button press in rbNoEPS_SimRRT_Nonlinear.
function rbNoEPS_SimRRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbNoEPS_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'off')
set(handles.rbYesEPS_SimRRT_Nonlinear, 'Enable', 'on')
set(handles.rbYesEPS_SimRRT_Nonlinear, 'Value', 0)

function txtSF_SimRRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtSF_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));
if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtSF_SimRRT_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSF_SimRRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popX0_Sim_Nonlinear.
function popX0_Sim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popX0_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
x0 = contents{get(hObject,'Value')};
x0_value = evalin('base', x0);
set(handles.txtX0_Sim_Nonlinear,'String', mat2str(x0_value));
handles.x0_nonlinear = x0;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popX0_Sim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popX0_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popU_Sim_Nonlinear.
function popU_Sim_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to popU_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
u = contents{get(hObject,'Value')};
u_value = evalin('base', u);
set(handles.txtU_Sim_Nonlinear,'String', mat2str(u_value));
handles.u_nonlinear = u;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popU_Sim_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popU_Sim_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txtBox_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtBox_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBox_Nonlinear as text
%        str2double(get(hObject,'String')) returns contents of txtBox_Nonlinear as a double


% --- Executes during object creation, after setting all properties.
function txtBox_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBox_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbEnterDynamicEquation_Nonlinear.
function rbEnterDynamicEquation_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbEnterDynamicEquation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dimensions_Nonlinear = [];
handles.nrOfInputs_Nonlinear = [];
handles.dynamicsEq_Nonlinear_error =0;
guidata(hObject, handles)

set(handles.popDim1Plot_Nonlinear, 'String', ' ');
set(handles.popDim1Plot_Nonlinear, 'Value', 1);
set(handles.popDim2Plot_Nonlinear, 'String', ' ');
set(handles.popDim2Plot_Nonlinear, 'Value', 1);

set(hObject, 'Enable', 'off')
set(handles.txtEnterDynamicEquation_Nonlinear, 'Enable', 'on')
set(handles.rbLoadSModel_Nonlinear, 'Enable', 'on')
set(handles.rbLoadSModel_Nonlinear, 'Value', 0)
set(handles.pbLoadSModel_Nonlinear, 'Enable', 'off')
set(handles.pbLoadSModel_Nonlinear, 'String', 'Load File')
set(handles.rbLoadDynamicEquation_Nonlinear, 'Enable', 'on')
set(handles.rbLoadDynamicEquation_Nonlinear, 'Value', 0)
set(handles.pbLoadDynamicEquation_Nonlinear, 'Enable', 'off')
set(handles.pbLoadDynamicEquation_Nonlinear, 'String', 'Load File')
set(handles.textFunction_Nonlinear, 'Enable', 'on')
set(handles.textEnd_Nonlinear, 'Enable', 'on')
txtEnterDynamicEquation_Nonlinear_Callback(handles.txt_equation, eventdata, handles)


function txtEnterDynamicEquation_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtEnterDynamicEquation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEnterDynamicEquation_Nonlinear as text
%        str2double(get(hObject,'String')) returns contents of txtEnterDynamicEquation_Nonlinear as a double
handles.txt_equation = hObject;
dynamics_equations = get(hObject, 'String');

if isempty(dynamics_equations)
    return
end

id = fopen([CORAROOT, filesep, 'models', filesep, 'auxiliary',filesep,'dynamicsEq_nonlinear.m'], 'wt');
fprintf(id, '%s\n', 'function f = dynamicsEq_nonlinear(x, u)');
for line = 1:size(dynamics_equations,1)
    fprintf(id, '%s', '    ');
    fprintf(id, '%s\n', dynamics_equations(line,:));
end
fprintf(id, '%s', 'end');
fclose(id);

handles.dynamicsEq_Nonlinear = [CORAROOT, filesep, 'models', filesep, 'auxiliary',filesep,'dynamicsEq_nonlinear.m'];

dynamics_equation = handles.dynamicsEq_Nonlinear;
dynamics_equation = strsplit(dynamics_equation, filesep);
dynamics_equation = dynamics_equation{end};
dynamics_equation = strsplit(dynamics_equation,'.');
dynamics_equation = dynamics_equation{1};
dynamicsEq_handle_str = ['@', dynamics_equation];
dynamicsEq_handle = eval('base', dynamicsEq_handle_str);
addpath([CORAROOT, filesep, 'models', filesep, 'auxiliary'])
handles.dynamicsEq_Nonlinear_error =0;

try
    [temp,dimensions]  = inputArgsLength_app(dynamicsEq_handle,2);
    nrOfInputs = max(1,temp(2));
catch
    message = sprintf('Please make sure the equations are written in the following form: \n f(1,1) = ... \n f(2,1) = ... \n\n staring with the letter f and all the parameters/variables are defined!');
    uiwait(msgbox(message, 'Error', 'error', 'modal'));
    return
end

handles.dynamicsEq_Nonlinear_error = 1;

%extract num_states, num_inputs
[temp,dimensions]  = inputArgsLength_app(dynamicsEq_handle,2);
nrOfInputs = max(1,temp(2));
dimensions = string(1:dimensions);
dimensions1 = ['', dimensions, 'time'];
dimensions2 = ['', dimensions];

handles.dimensions_Nonlinear = dimensions;

set(handles.popDim1Plot_Nonlinear, 'String', dimensions1);
set(handles.popDim1Plot_Nonlinear, 'Value', 1);
set(handles.popDim2Plot_Nonlinear, 'String', dimensions2);
set(handles.popDim2Plot_Nonlinear, 'Value', 1);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function txtEnterDynamicEquation_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEnterDynamicEquation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbLoadDynamicEquation_Nonlinear.
function rbLoadDynamicEquation_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbLoadDynamicEquation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dimensions_Nonlinear = [];
guidata(hObject, handles)

set(handles.popDim1Plot_Nonlinear, 'String', ' ');
set(handles.popDim1Plot_Nonlinear, 'Value', 1);
set(handles.popDim2Plot_Nonlinear, 'String', ' ');
set(handles.popDim2Plot_Nonlinear, 'Value', 1);

set(hObject, 'Enable', 'off')
set(handles.pbLoadDynamicEquation_Nonlinear, 'Enable', 'on')
set(handles.rbLoadSModel_Nonlinear, 'Enable', 'on')
set(handles.rbLoadSModel_Nonlinear, 'Value', 0)
set(handles.pbLoadSModel_Nonlinear, 'Enable', 'off')
set(handles.pbLoadSModel_Nonlinear, 'String', 'Load File')
set(handles.rbEnterDynamicEquation_Nonlinear, 'Enable', 'on')
set(handles.rbEnterDynamicEquation_Nonlinear, 'Value', 0)
set(handles.txtEnterDynamicEquation_Nonlinear, 'Enable', 'off')
%set(handles.txtEnterDynamicEquation_Nonlinear, 'String', ' ')
set(handles.textFunction_Nonlinear, 'Enable', 'off')
set(handles.textEnd_Nonlinear, 'Enable', 'off')




% --- Executes on button press in pbLoadDynamicEquation_Nonlinear.
function pbLoadDynamicEquation_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadDynamicEquation_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dynamicsEq_file, dynamicsEq_path] = uigetfile('*.m');
set(hObject, 'string', dynamicsEq_file)
handles.dynamicsEq_Nonlinear = fullfile(dynamicsEq_path, dynamicsEq_file);

dynamics_equation = handles.dynamicsEq_Nonlinear;
dynamics_equation = strsplit(dynamics_equation, filesep);
dynamics_equation = dynamics_equation{end};
dynamics_equation = strsplit(dynamics_equation,'.');
dynamics_equation = dynamics_equation{1};
dynamicsEq_handle_str = ['@', dynamics_equation];
try
    dynamicsEq_handle = eval('base', dynamicsEq_handle_str);
catch
    return
end


%extract num_states, num_inputs
[temp,dimensions_eq]  = inputArgsLength_app(dynamicsEq_handle,2);
nrOfInputs = max(1,temp(2));
dimensions = string(1:dimensions_eq);
dimensions1 = ['', dimensions, 'time'];
dimensions2 = ['', dimensions];

handles.dimensions_Nonlinear = dimensions_eq;
handles.dynamicsEq_Nonlinear_error = 1;
guidata(hObject, handles)

set(handles.popDim1Plot_Nonlinear, 'String', dimensions1);
set(handles.popDim1Plot_Nonlinear, 'Value', 1);
set(handles.popDim2Plot_Nonlinear, 'String', dimensions2);
set(handles.popDim2Plot_Nonlinear, 'Value', 1);

guidata(hObject, handles) %update handles



function txtMaxLinError_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtMaxLinError_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaxLinError_Nonlinear as text
%        str2double(get(hObject,'String')) returns contents of txtMaxLinError_Nonlinear as a double
value = str2double(get(hObject, 'String'));
if value <= 0
    value = 0.1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtMaxLinError_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMaxLinError_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownMaxLinError_Nonlinear.
function pbDownMaxLinError_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownMaxLinError_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtMaxLinError_Nonlinear, 'String');
value = str2double(value);

if value <= 1
    value = value - 0.1;
else
    value = value - 1;
end

if value > 0
    value = num2str(value);
    set(handles.txtMaxLinError_Nonlinear, 'String', value);
end

% --- Executes on button press in pbUpMaxLinError_Nonlinear.
function pbUpMaxLinError_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpMaxLinError_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtMaxLinError_Nonlinear, 'String');
value = str2double(value);

if value < 1
    value = 0;
end

value = value + 1;

value = num2str(value);
set(handles.txtMaxLinError_Nonlinear, 'String', value);


function txtRI_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtRI_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtRI_Nonlinear as text
%        str2double(get(hObject,'String')) returns contents of txtRI_Nonlinear as a double
value = int16(str2double(get(hObject, 'String')));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtRI_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtRI_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownRI_Nonlinear.
function pbDownRI_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownRI_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtRI_Nonlinear, 'String');
value = int16(str2double(value)) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtRI_Nonlinear, 'String', value);
end


% --- Executes on button press in pbUpRI_Nonlinear.
function pbUpRI_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpRI_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtRI_Nonlinear, 'String');
value = int16(str2double(value)) + 1;
value = num2str(value);
set(handles.txtRI_Nonlinear, 'String', value);


% --- Executes on button press in rbConLin_Nonlinear.
function rbConLin_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbConLin_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Enable', 'off')
set(handles.txtConLinTO_Nonlinear, 'Enable', 'on')
set(handles.pbUpConLinTO_Nonlinear, 'Enable', 'on')
set(handles.pbDownConLinTO_Nonlinear, 'Enable', 'on')

set(handles.rbConPol_Nonlinear, 'Enable', 'on')
set(handles.rbConPol_Nonlinear, 'Value', 0)
set(handles.txtConPolTO_Nonlinear, 'Enable', 'off')
set(handles.pbUpConPolTO_Nonlinear, 'Enable', 'off')
set(handles.pbDownConPolTO_Nonlinear, 'Enable', 'off')
set(handles.txtIO_Nonlinear, 'Enable', 'off')
set(handles.pbUpIO_Nonlinear, 'Enable', 'off')
set(handles.pbDownIO_Nonlinear, 'Enable', 'off')
set(handles.txtEO_Nonlinear, 'Enable', 'off')
set(handles.pbUpEO_Nonlinear, 'Enable', 'off')
set(handles.pbDownEO_Nonlinear, 'Enable', 'off')

value = int16(str2double(get(handles.txtConLinTO_Nonlinear, 'String')));

if value < 3
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'off')
else
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'on')
end


function txtConLinTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtConLinTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = int16(str2double(get(hObject, 'String')));

if value < 3
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'off')
else
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'on')
end

if value <= 3
    value = num2str(value);
    set(hObject, 'String', value);
end


% --- Executes during object creation, after setting all properties.
function txtConLinTO_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtConLinTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownConLinTO_Nonlinear.
function pbDownConLinTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownConLinTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConLinTO_Nonlinear, 'String');
value = str2double(value) - 1;

if value < 3
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'off')
else
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'on')
end

if value >=2
    value = num2str(value);
    set(handles.txtConLinTO_Nonlinear, 'String', value);
end


% --- Executes on button press in pbUpConLinTO_Nonlinear.
function pbUpConLinTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpConLinTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConLinTO_Nonlinear, 'String');
value = str2double(value) + 1;

if value < 3
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'off')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'off')
else
    set(handles.txtIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.txtEO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'on')
    set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'on')
end

if value <= 3
    value = num2str(value);
    set(handles.txtConLinTO_Nonlinear, 'String', value);
end


% --- Executes on button press in rbConPol_Nonlinear.
function rbConPol_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to rbConPol_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'off')
set(handles.txtConLinTO_Nonlinear, 'Enable', 'off')
set(handles.pbUpConLinTO_Nonlinear, 'Enable', 'off')
set(handles.pbDownConLinTO_Nonlinear, 'Enable', 'off')
set(handles.txtIO_Nonlinear_lin, 'Enable', 'off')
set(handles.txtEO_Nonlinear_lin, 'Enable', 'off')
set(handles.pbUpIO_Nonlinear_lin, 'Enable', 'off')
set(handles.pbDownIO_Nonlinear_lin, 'Enable', 'off')
set(handles.pbUpEO_Nonlinear_lin, 'Enable', 'off')
set(handles.pbDownEO_Nonlinear_lin, 'Enable', 'off')

set(handles.rbConLin_Nonlinear, 'Enable', 'on')
set(handles.rbConLin_Nonlinear, 'Value', 0)
set(handles.txtConPolTO_Nonlinear, 'Enable', 'on')
set(handles.pbUpConPolTO_Nonlinear, 'Enable', 'on')
set(handles.pbDownConPolTO_Nonlinear, 'Enable', 'on')
set(handles.txtIO_Nonlinear, 'Enable', 'on')
set(handles.pbUpIO_Nonlinear, 'Enable', 'on')
set(handles.pbDownIO_Nonlinear, 'Enable', 'on')
set(handles.txtEO_Nonlinear, 'Enable', 'on')
set(handles.pbUpEO_Nonlinear, 'Enable', 'on')
set(handles.pbDownEO_Nonlinear, 'Enable', 'on')


function txtConPolTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtConPolTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = int16(str2double(get(hObject, 'String')));

if value < 3
    value = 3;
end

if value > 4
    value = 4;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtConPolTO_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtConPolTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownConPolTO_Nonlinear.
function pbDownConPolTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownConPolTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConPolTO_Nonlinear, 'String');
value = str2double(value) - 1;

if value >=3
    value = num2str(value);
    set(handles.txtConPolTO_Nonlinear, 'String', value);
end


% --- Executes on button press in pbUpConPolTO_Nonlinear.
function pbUpConPolTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpConPolTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConPolTO_Nonlinear, 'String');
value = str2double(value) + 1;

if value <= 4
    value = num2str(value);
    set(handles.txtConPolTO_Nonlinear, 'String', value);
end


function txtIO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtIO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = int16(str2double(get(hObject, 'String')));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtIO_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownIO_Nonlinear.
function pbDownIO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownIO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Nonlinear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtIO_Nonlinear, 'String', value);
end

% --- Executes on button press in pbUpIO_Nonlinear.
function pbUpIO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpIO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Nonlinear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtIO_Nonlinear, 'String', value);


function txtEO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to txtEO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = int16(str2double(get(hObject, 'String')));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtEO_Nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownEO_Nonlinear.
function pbDownEO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownEO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtEO_Nonlinear, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtEO_Nonlinear, 'String', value);
end

% --- Executes on button press in pbUpEO_Nonlinear.
function pbUpEO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpEO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtEO_Nonlinear, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtEO_Nonlinear, 'String', value);


% --- Executes on selection change in popA_Linear.
function popA_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
A = contents{get(hObject,'Value')};

A_value = evalin('base', A);
dimensions = string(1:size(A_value,1));
dimensions1 = ['', dimensions, 'time'];
dimensions2 = ['', dimensions];

handles.dimensions_Linear = size(A_value,1);
guidata(hObject, handles)

set(handles.popDim1Plot_Linear, 'String', dimensions1);
set(handles.popDim1Plot_Linear, 'Value', 1);
set(handles.popDim2Plot_Linear, 'String', dimensions2);
set(handles.popDim2Plot_Linear, 'Value', 1);

set(handles.textA_Linear,'String', mat2str(A_value));
handles.A = A;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popA_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popD_Linear.
function popD_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popD_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
D = contents{get(hObject,'Value')};
D_value = evalin('base', D);
set(handles.textD_Linear,'String', mat2str(D_value));
handles.D = D;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popD_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popD_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popB_Linear.
function popB_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popB_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
B = contents{get(hObject,'Value')};
B_value = evalin('base', B);
set(handles.textB_Linear,'String', mat2str(B_value));
handles.B = B;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popB_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popB_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popE_Linear.
function popE_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popE_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
E = contents{get(hObject,'Value')};
E_value = evalin('base', E);
set(handles.textE_Linear,'String', mat2str(E_value));
handles.E = E;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popE_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popE_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popC_Linear.
function popC_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popC_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
C = contents{get(hObject,'Value')};
C_value = evalin('base', C);
set(handles.textC_Linear,'String', mat2str(C_value));
handles.C = C;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popC_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popC_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popF_Linear.
function popF_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to popF_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
F = contents{get(hObject,'Value')};
F_value = evalin('base', F);
set(handles.textF_Linear,'String', mat2str(F_value));
handles.F = F;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popF_Linear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popF_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_infoZO_Linear.
function pb_infoZO_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoZO_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_zonotopeOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pb_infoTT_Linear.
function pb_infoTT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoTT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_taylorTerms.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pb_infoRT_Linear.
function pb_infoRT_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoRT_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_reductionTechnique.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pbInfoLA_Linear.
function pbInfoLA_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pbInfoLA_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_linAlg.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pb_infoZO_Nonlinear.
function pb_infoZO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoZO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_zonotopeOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoTT_Nonlinear.
function pb_infoTT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoTT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_taylorTerms.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoRT_Nonlinear.
function pb_infoRT_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoRT_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_reductionTechnique.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoMLE_Nonlinear.
function pb_infoMLE_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoMLE_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'InfoSaveData.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoRI_Nonlinear.
function pb_infoRI_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoRI_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'InfoSaveData.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCL_Nonlinear.
function pb_infoCL_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCL_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_consLinearization.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCP_Nonlinear.
function pb_infoCP_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCP_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_consPolynomialization.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCLTO_Nonlinear.
function pb_infoCLTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCLTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_tensorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCPTO_Nonlinear.
function pb_infoCPTO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPTO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_tensorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCPIO_Nonlinear.
function pb_infoCPIO_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPIO_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_intermediateOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCPEE_Nonlinear.
function pb_infoCPEE_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPEE_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_errorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pushDeletePlot_Nonlinear.
function pushDeletePlot_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pushDeletePlot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.listPlots_Nonlinear,'String'));
value = get(handles.listPlots_Nonlinear, 'Value');

if isempty(contents)
    return
end

contents(value,:) = [];

handles.nonlinear_plots = contents';
guidata(hObject, handles)
set(handles.listPlots_Nonlinear, 'Value', 1);
set(handles.listPlots_Nonlinear, 'String', contents);



function txtTstart_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtTstart_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTstart_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtTstart_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtTstart_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTstart_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTfinal_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtTfinal_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTfinal_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtTfinal_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtTfinal_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTfinal_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtZOrder_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtZOrder_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtZOrder_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtZOrder_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtZOrder_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtZOrder_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownZOrder_Hybrid.
function pbDownZOrder_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownZOrder_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtZOrder_Hybrid, 'String');
value = str2double(value) - 1;
value = num2str(value);
set(handles.txtZOrder_Hybrid, 'String', value);


% --- Executes on button press in pbUpZOrder_Hybrid.
function pbUpZOrder_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpZOrder_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtZOrder_Hybrid, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtZOrder_Hybrid, 'String', value);


function txtTaylor_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtTaylor_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTaylor_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtTaylor_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtTaylor_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTaylor_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbUpTaylor_Hybrid.
function pbUpTaylor_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpTaylor_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtTaylor_Hybrid, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtTaylor_Hybrid, 'String', value);


% --- Executes on button press in pbDownTaylor_Hybrid.
function pbDownTaylor_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownTaylor_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtTaylor_Hybrid, 'String');
value = str2double(value) - 1;
value = num2str(value);
set(handles.txtTaylor_Hybrid, 'String', value);

% --- Executes on selection change in popRT_Hybrid.
function popRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popRT_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popRT_Hybrid


% --- Executes during object creation, after setting all properties.
function popRT_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLA_Hybrid.
function popLA_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popLA_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject, 'String'));
alg = contents{get(hObject,'Value')};

if strcmp(alg, 'adaptive')
    set(handles.OptionsError_Hybrid, 'visible', 'on')
    set(handles.txtOptionError_Hybrid, 'visible', 'on')
    set(handles.txtOptionErrorUp_Hybrid, 'visible', 'on')
    set(handles.txtOptionErrorDown_Hybrid, 'visible', 'on')
    set(handles.txtZOrder_Hybrid, 'visible', 'off')
    set(handles.pbUpZOrder_Hybrid, 'visible', 'off')
    set(handles.pbDownZOrder_Hybrid, 'visible', 'off')
    set(handles.txtTaylor_Hybrid, 'visible', 'off')
    set(handles.pbUpTaylor_Hybrid, 'visible', 'off')
    set(handles.pbDownTaylor_Hybrid, 'visible', 'off')
    set(handles.popRT_Hybrid, 'visible', 'off')
    set(handles.zotxt_Hybrid, 'visible', 'off')
    set(handles.tttxt_Hybrid, 'visible', 'off')
    set(handles.rttxt_Hybrid, 'visible', 'off')
    set(handles.pb_infoZO_Hybrid, 'visible', 'off')
    set(handles.pb_infoTT_Hybrid, 'visible', 'off')
    set(handles.pb_infoRT_Hybrid, 'visible', 'off')
    set(handles.txtTimeStep_Hybrid, 'visible', 'off')
    set(handles.tstxt_Hybrid, 'visible', 'off')
else
    set(handles.OptionsError_Hybrid, 'visible', 'off')
    set(handles.txtOptionError_Hybrid, 'visible', 'off')
    set(handles.txtOptionErrorUp_Hybrid, 'visible', 'off')
    set(handles.txtOptionErrorDown_Hybrid, 'visible', 'off')
    
    set(handles.txtZOrder_Hybrid, 'visible', 'on')
    set(handles.pbUpZOrder_Hybrid, 'visible', 'on')
    set(handles.pbDownZOrder_Hybrid, 'visible', 'on')
    set(handles.txtTaylor_Hybrid, 'visible', 'on')
    set(handles.pbUpTaylor_Hybrid, 'visible', 'on')
    set(handles.pbDownTaylor_Hybrid, 'visible', 'on')
    set(handles.popRT_Hybrid, 'visible', 'on')
    set(handles.zotxt_Hybrid, 'visible', 'on')
    set(handles.tttxt_Hybrid, 'visible', 'on')
    set(handles.rttxt_Hybrid, 'visible', 'on')
    set(handles.pb_infoZO_Hybrid, 'visible', 'on')
    set(handles.pb_infoTT_Hybrid, 'visible', 'on')
    set(handles.pb_infoRT_Hybrid, 'visible', 'on')
    set(handles.txtTimeStep_Hybrid, 'visible', 'on')
    set(handles.tstxt_Hybrid, 'visible', 'on')

end

% --- Executes during object creation, after setting all properties.
function popLA_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLA_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbR0_Hybrid.
function pbR0_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbR0_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    [R0, handles_initial_set] = input_set(handles.initial_set_hybrid,handles.initial_set_hybrid_type);
    handles.R0_Hybrid = R0;
    handles.initial_set_hybrid = handles_initial_set;
    handles.initial_set_hybrid_type = R0.type;
    guidata(hObject, handles)
catch
end


% --- Executes on button press in pbU_Hybrid.
function pbU_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbU_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.popULoc_Hybrid,'String'));
loc = contents{get(handles.popULoc_Hybrid,'Value')};
if isempty(loc)
    return
elseif loc ~= 'all'
    loc = str2double(loc)+1;
else
    loc = 1;
end

try
    [U, handles_input_set{loc}] = input_set(handles.input_set_hybrid{loc},handles.input_set_hybrid_type{loc});
    handles.U_Hybrid{loc} = U;
    handles.input_set_hybrid{loc} = handles_input_set{loc};
    handles.input_set_hybrid_type{loc} = U.type;
    guidata(hObject, handles)
catch
end



% --- Executes on button press in pbLoadSModel_Hybrid.
function pbLoadSModel_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadSModel_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[spaceEx_file, spaceEx_path] = uigetfile('*.xml');
try
    handles.spaceExFile_Hybrid = fullfile(spaceEx_path, spaceEx_file);
    spaceex2cora(handles.spaceExFile_Hybrid,0);
catch
    return
end
    

file = strsplit(handles.spaceExFile_Hybrid, filesep);
file = file{end};
file = strsplit(file, '.');
file = file{1};
hybridSys =  eval(strcat(file, '()'));

if ~strcmp(class(hybridSys),'linearSys') && ~strcmp(class(hybridSys),'nonlinearSys')
    dimensions = string(1:hybridSys.location{1,1}.contDynamics.dim);
    dimensions1 = ['', dimensions, 'time'];
    dimensions2 = ['', dimensions];

    handles.dimensions_Hybrid = hybridSys.location{1,1}.contDynamics.dim;
    handles.nrOfInputs_Hybrid = hybridSys.location{1, 1}.contDynamics.nrOfInputs;

    set(handles.popDim1Plot_Hybrid, 'String', dimensions1);
    set(handles.popDim1Plot_Hybrid, 'Value', 1);
    set(handles.popDim2Plot_Hybrid, 'String', dimensions2);
    set(handles.popDim2Plot_Hybrid, 'Value', 1);
    set(hObject, 'string', spaceEx_file)
    
    temp1 = false; temp2 = false;
    list = {'polytope';'conZonotope';'zonoGirard';'hyperplaneMap'; ...
            'pancake';'nondetGuard';'levelSet'};
    for i = 1:length(hybridSys.location)
        for j = 1:length(hybridSys.location{i}.transition)
            if isa(hybridSys.location{i}.transition{j}.guard,'levelSet')
                temp1 = true;
            else
                temp2 = true;
            end
        end
    end
    if temp1 && temp2
        uiwait(msgbox('Hybrid automata with mixed guard sets consisting of level set and other guards are not supported', 'Error', 'error', 'modal'))
        return
    elseif temp2
        set(handles.popMethod_Hybrid,'String',list(1:end-1));
        set(handles.pb_info_enclose_Hybrid, 'Enable', 'on')
        set(handles.checkbox_Hybrid, 'Enable', 'on')
        set(handles.checkpca_Hybrid, 'Enable', 'on')
        set(handles.checkflow_Hybrid, 'Enable', 'on')
        set(handles.txtGO_Hybrid, 'Enable', 'off')
        set(handles.GOtxt_Hybrid, 'Enable', 'off')
    else
        set(handles.popMethod_Hybrid,'String',list(end));
        set(handles.pb_info_enclose_Hybrid, 'Enable', 'off')
        set(handles.checkbox_Hybrid, 'Enable', 'off')
        set(handles.checkpca_Hybrid, 'Enable', 'off')
        set(handles.checkflow_Hybrid, 'Enable', 'off')
        set(handles.txtGO_Hybrid, 'Enable', 'off')
        set(handles.GOtxt_Hybrid, 'Enable', 'off')
    end
    guidata(hObject, handles) %update handles   

    %display locations in hybrid --------------------------------------
    loc_size = size(hybridSys.location,2);
    handles.U_hybrid = cell(1,loc_size+1);
    handles.input_set_hybrid = cell(1,loc_size+1);
    handles.input_set_hybrid_type = cell(1,loc_size+1);
    guidata(hObject, handles)

    locs = string(1:loc_size);
    
    set(handles.popStartLoc_Hybrid,'String',locs)
    set(handles.popFinalLoc_Hybrid,'String',['0', locs])

    if size(hybridSys.location,2) ~= 1
        locs = ['all', locs];
    end
    set(handles.popULoc_Hybrid,'String',locs)
    handles.locs_hybrid = locs;

    contents = cellstr(get(handles.popULoc_Hybrid,'String'));
    loc = contents{get(handles.popULoc_Hybrid,'Value')};

    set(handles.pbU_Hybrid,'String',['Compose Input for Loc ',num2str(loc)])

    %display linear and/or nonlinear settings -------------------------
    
    position = get(handles.txtZOrder_Hybrid, 'Position');
    if position(2) < 37
        position = get(handles.txtZOrder_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtZOrder_Hybrid, 'Position', position );
        position = get(handles.pbUpZOrder_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpZOrder_Hybrid, 'Position', position );
        position = get(handles.pbDownZOrder_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownZOrder_Hybrid, 'Position', position );
        position = get(handles.txtTaylor_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtTaylor_Hybrid, 'Position', position );
        position = get(handles.pbUpTaylor_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpTaylor_Hybrid, 'Position', position );
        position = get(handles.pbDownTaylor_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownTaylor_Hybrid, 'Position', position );
        position = get(handles.popRT_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.popRT_Hybrid, 'Position', position );
        position = get(handles.zotxt_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.zotxt_Hybrid, 'Position', position );
        position = get(handles.tttxt_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.tttxt_Hybrid, 'Position', position );
        position = get(handles.rttxt_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.rttxt_Hybrid, 'Position', position );
        position = get(handles.pb_infoZO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoZO_Hybrid, 'Position', position );
        position = get(handles.pb_infoTT_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoTT_Hybrid, 'Position', position );
        position = get(handles.pb_infoRT_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoRT_Hybrid, 'Position', position );
        position = get(handles.txtTimeStep_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtTimeStep_Hybrid, 'Position', position );
        position = get(handles.tstxt_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.tstxt_Hybrid, 'Position', position );
        position = get(handles.txtTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtTO_Hybrid, 'Position', position );
        position = get(handles.txtTO2_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtTO2_Hybrid, 'Position', position );
        position = get(handles.txtIO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtIO_Hybrid, 'Position', position );
        position = get(handles.txtIO2_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtIO2_Hybrid, 'Position', position );
        position = get(handles.txtEO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtEO_Hybrid, 'Position', position );
        position = get(handles.txtEO2_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtEO2_Hybrid, 'Position', position );
        position = get(handles.txtCP_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtCP_Hybrid, 'Position', position );
        position = get(handles.rbConLin_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.rbConLin_Hybrid, 'Position', position );
        position = get(handles.txtConLinTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtConLinTO_Hybrid, 'Position', position );
        position = get(handles.pbUpConLinTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpConLinTO_Hybrid, 'Position', position );
        position = get(handles.pbDownConLinTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownConLinTO_Hybrid, 'Position', position );
        position = get(handles.txtIO_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.txtIO_Hybrid_lin, 'Position', position );
        position = get(handles.txtEO_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.txtEO_Hybrid_lin, 'Position', position );
        position = get(handles.pbUpIO_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpIO_Hybrid_lin, 'Position', position );
        position = get(handles.pbDownIO_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownIO_Hybrid_lin, 'Position', position );
        position = get(handles.pbUpEO_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpEO_Hybrid_lin, 'Position', position );
        position = get(handles.pbDownEO_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownEO_Hybrid_lin, 'Position', position );
        position = get(handles.pb_infoCPIO_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCPIO_Hybrid_lin, 'Position', position );
        position = get(handles.pb_infoCPEE_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCPEE_Hybrid_lin, 'Position', position );
        position = get(handles.txtIO2_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.txtIO2_Hybrid_lin, 'Position', position );
        position = get(handles.txtEO2_Hybrid_lin, 'Position')+ [ 0 1 0 0];
        set(handles.txtEO2_Hybrid_lin, 'Position', position );
        position = get(handles.pb_infoCL_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCL_Hybrid, 'Position', position );
        position = get(handles.pb_infoCLTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCLTO_Hybrid, 'Position', position );
        position = get(handles.pb_infoCP_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCP_Hybrid, 'Position', position );
        position = get(handles.pb_infoCPTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCPTO_Hybrid, 'Position', position );
        position = get(handles.pb_infoCPIO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCPIO_Hybrid, 'Position', position );
        position = get(handles.pb_infoCPEE_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pb_infoCPEE_Hybrid, 'Position', position );
        position = get(handles.txtConPolTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.txtConPolTO_Hybrid, 'Position', position );
        position = get(handles.pbUpConPolTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpConPolTO_Hybrid, 'Position', position );
        position = get(handles.pbDownConPolTO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownConPolTO_Hybrid, 'Position', position );
        position = get(handles.pbUpIO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpIO_Hybrid, 'Position', position );
        position = get(handles.pbDownIO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownIO_Hybrid, 'Position', position );
        position = get(handles.pbUpEO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbUpEO_Hybrid, 'Position', position );
        position = get(handles.pbDownEO_Hybrid, 'Position')+ [ 0 1 0 0];
        set(handles.pbDownEO_Hybrid, 'Position', position );
    end

     for i = 1 : length(hybridSys.location)
        if isa(hybridSys.location{i}.contDynamics,'linearSys')
            
            set(handles.popLA_Hybrid, 'visible', 'on')
            set(handles.txtA_Hybrid, 'visible', 'on')
            set(handles.pbInfoLA_Hybrid, 'visible', 'on')
            set(handles.OptionsError_Hybrid, 'visible', 'on')
            set(handles.txtOptionError_Hybrid, 'visible', 'on')
            set(handles.txtOptionErrorUp_Hybrid, 'visible', 'on')
            set(handles.txtOptionErrorDown_Hybrid, 'visible', 'on')
            
            set(handles.txtZOrder_Hybrid, 'visible', 'off')
            set(handles.pbUpZOrder_Hybrid, 'visible', 'off')
            set(handles.pbDownZOrder_Hybrid, 'visible', 'off')
            set(handles.txtTaylor_Hybrid, 'visible', 'off')
            set(handles.pbUpTaylor_Hybrid, 'visible', 'off')
            set(handles.pbDownTaylor_Hybrid, 'visible', 'off')
            set(handles.popRT_Hybrid, 'visible', 'off')
            set(handles.zotxt_Hybrid, 'visible', 'off')
            set(handles.tttxt_Hybrid, 'visible', 'off')
            set(handles.rttxt_Hybrid, 'visible', 'off')
            set(handles.pb_infoZO_Hybrid, 'visible', 'off')
            set(handles.pb_infoTT_Hybrid, 'visible', 'off')
            set(handles.pb_infoRT_Hybrid, 'visible', 'off')
            set(handles.txtTimeStep_Hybrid, 'visible', 'off')
            set(handles.tstxt_Hybrid, 'visible', 'off')

            set(handles.txtTO_Hybrid, 'visible', 'off')
            set(handles.txtTO2_Hybrid, 'visible', 'off')
            set(handles.txtIO_Hybrid, 'visible', 'off')
            set(handles.txtIO2_Hybrid, 'visible', 'off')
            set(handles.txtEO_Hybrid, 'visible', 'off')
            set(handles.txtEO2_Hybrid, 'visible', 'off')
            set(handles.txtCP_Hybrid, 'visible', 'off')

            set(handles.rbConLin_Hybrid, 'visible', 'off')
            set(handles.txtConLinTO_Hybrid, 'visible', 'off')
            set(handles.pbUpConLinTO_Hybrid, 'visible', 'off')
            set(handles.pbDownConLinTO_Hybrid, 'visible', 'off')
            set(handles.txtIO_Hybrid_lin, 'visible', 'off')
            set(handles.txtEO_Hybrid_lin, 'visible', 'off')
            set(handles.pbUpIO_Hybrid_lin, 'visible', 'off')
            set(handles.pbDownIO_Hybrid_lin, 'visible', 'off')
            set(handles.pbUpEO_Hybrid_lin, 'visible', 'off')
            set(handles.pbDownEO_Hybrid_lin, 'visible', 'off')
            set(handles.pb_infoCPIO_Hybrid_lin, 'visible', 'off')
            set(handles.pb_infoCPEE_Hybrid_lin, 'visible', 'off')
            set(handles.txtIO2_Hybrid_lin, 'visible', 'off')
            set(handles.txtEO2_Hybrid_lin, 'visible', 'off')

            set(handles.pb_infoCL_Hybrid, 'visible', 'off')
            set(handles.pb_infoCLTO_Hybrid, 'visible', 'off')
            set(handles.pb_infoCP_Hybrid, 'visible', 'off')
            set(handles.pb_infoCPTO_Hybrid, 'visible', 'off')
            set(handles.pb_infoCPIO_Hybrid, 'visible', 'off')
            set(handles.pb_infoCPEE_Hybrid, 'visible', 'off')
            set(handles.txtConPolTO_Hybrid, 'visible', 'off')
            set(handles.pbUpConPolTO_Hybrid, 'visible', 'off')
            set(handles.pbDownConPolTO_Hybrid, 'visible', 'off')
            set(handles.pbUpIO_Hybrid, 'visible', 'off')
            set(handles.pbDownIO_Hybrid, 'visible', 'off')
            set(handles.pbUpEO_Hybrid, 'visible', 'off')
            set(handles.pbDownEO_Hybrid, 'visible', 'off')
        else
            set(handles.popLA_Hybrid, 'visible', 'off')
            set(handles.txtA_Hybrid, 'visible', 'off')
            set(handles.pbInfoLA_Hybrid, 'visible', 'off')
            set(handles.OptionsError_Hybrid, 'visible', 'off')
            set(handles.txtOptionError_Hybrid, 'visible', 'off')
            set(handles.txtOptionErrorUp_Hybrid, 'visible', 'off')
            set(handles.txtOptionErrorDown_Hybrid, 'visible', 'off')
            
            set(handles.txtZOrder_Hybrid, 'visible', 'on')
            set(handles.pbUpZOrder_Hybrid, 'visible', 'on')
            set(handles.pbDownZOrder_Hybrid, 'visible', 'on')
            set(handles.txtTaylor_Hybrid, 'visible', 'on')
            set(handles.pbUpTaylor_Hybrid, 'visible', 'on')
            set(handles.pbDownTaylor_Hybrid, 'visible', 'on')
            set(handles.popRT_Hybrid, 'visible', 'on')
            set(handles.zotxt_Hybrid, 'visible', 'on')
            set(handles.tttxt_Hybrid, 'visible', 'on')
            set(handles.rttxt_Hybrid, 'visible', 'on')
            set(handles.pb_infoZO_Hybrid, 'visible', 'on')
            set(handles.pb_infoTT_Hybrid, 'visible', 'on')
            set(handles.pb_infoRT_Hybrid, 'visible', 'on')
            set(handles.txtTimeStep_Hybrid, 'visible', 'on')
            set(handles.tstxt_Hybrid, 'visible', 'on')
            set(handles.txtTO_Hybrid, 'visible', 'on')
            set(handles.txtTO2_Hybrid, 'visible', 'on')
            set(handles.txtIO_Hybrid, 'visible', 'on')
            set(handles.txtIO2_Hybrid, 'visible', 'on')
            set(handles.txtEO_Hybrid, 'visible', 'on')
            set(handles.txtEO2_Hybrid, 'visible', 'on')
            set(handles.txtCP_Hybrid, 'visible', 'on')
            set(handles.rbConLin_Hybrid, 'visible', 'on')
            set(handles.txtConLinTO_Hybrid, 'visible', 'on')
            set(handles.pbUpConLinTO_Hybrid, 'visible', 'on')
            set(handles.pbDownConLinTO_Hybrid, 'visible', 'on')
            set(handles.txtIO_Hybrid_lin, 'visible', 'on')
            set(handles.txtEO_Hybrid_lin, 'visible', 'on')
            set(handles.pbUpIO_Hybrid_lin, 'visible', 'on')
            set(handles.pbDownIO_Hybrid_lin, 'visible', 'on')
            set(handles.pbUpEO_Hybrid_lin, 'visible', 'on')
            set(handles.pbDownEO_Hybrid_lin, 'visible', 'on')
            set(handles.pb_infoCPIO_Hybrid_lin, 'visible', 'on')
            set(handles.pb_infoCPEE_Hybrid_lin, 'visible', 'on')
            set(handles.txtIO2_Hybrid_lin, 'visible', 'on')
            set(handles.txtEO2_Hybrid_lin, 'visible', 'on')
            set(handles.pb_infoCL_Hybrid, 'visible', 'on')
            set(handles.pb_infoCLTO_Hybrid, 'visible', 'on')
            set(handles.pb_infoCP_Hybrid, 'visible', 'on')
            set(handles.pb_infoCPTO_Hybrid, 'visible', 'on')
            set(handles.pb_infoCPIO_Hybrid, 'visible', 'on')
            set(handles.pb_infoCPEE_Hybrid, 'visible', 'on')
            set(handles.txtConPolTO_Hybrid, 'visible', 'on')
            set(handles.pbUpConPolTO_Hybrid, 'visible', 'on')
            set(handles.pbDownConPolTO_Hybrid, 'visible', 'on')
            set(handles.pbUpIO_Hybrid, 'visible', 'on')
            set(handles.pbDownIO_Hybrid, 'visible', 'on')
            set(handles.pbUpEO_Hybrid, 'visible', 'on')
            set(handles.pbDownEO_Hybrid, 'visible', 'on')
            
        end
        
     end
else
    uiwait(msgbox('The loaded SpaceEx file is not hybrid', 'Error', 'error', 'modal'))
    return
end


% --- Executes on selection change in popUTransLoc_hybrid.
function popUTransLoc_hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popUTransLoc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUTransLoc_hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUTransLoc_hybrid


% --- Executes during object creation, after setting all properties.
function popUTransLoc_hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUTransLoc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu85.
function popupmenu85_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu85 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu85


% --- Executes during object creation, after setting all properties.
function popupmenu85_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtUTransLoc_hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtUTransLoc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtUTransLoc_hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtUTransLoc_hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtUTransLoc_hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtUTransLoc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit96_Callback(hObject, eventdata, handles)
% hObject    handle to edit96 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit96 as text
%        str2double(get(hObject,'String')) returns contents of edit96 as a double


% --- Executes during object creation, after setting all properties.
function edit96_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit96 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_infoZO_Hybrid.
function pb_infoZO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoZO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_zonotopeOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoTT_Hybrid.
function pb_infoTT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoTT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_taylorTerms.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoRT_Hybrid.
function pb_infoRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_reductionTechnique.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pbInfoLA_Hybrid.
function pbInfoLA_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbInfoLA_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_linAlg.png';
infoBox({[path_im, im_ZO]});
uiwait;


function txtTimeStep_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtTimeStep_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimeStep_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtTimeStep_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtTimeStep_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimeStep_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton45.
function radiobutton45_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton45



function edit98_Callback(hObject, eventdata, handles)
% hObject    handle to edit98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit98 as text
%        str2double(get(hObject,'String')) returns contents of edit98 as a double


% --- Executes during object creation, after setting all properties.
function edit98_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit98 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox61.
function checkbox61_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox61



function edit99_Callback(hObject, eventdata, handles)
% hObject    handle to edit99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit99 as text
%        str2double(get(hObject,'String')) returns contents of edit99 as a double


% --- Executes during object creation, after setting all properties.
function edit99_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox62.
function checkbox62_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox62



function edit100_Callback(hObject, eventdata, handles)
% hObject    handle to edit100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit100 as text
%        str2double(get(hObject,'String')) returns contents of edit100 as a double


% --- Executes during object creation, after setting all properties.
function edit100_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox63.
function checkbox63_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox63



function edit101_Callback(hObject, eventdata, handles)
% hObject    handle to edit101 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit101 as text
%        str2double(get(hObject,'String')) returns contents of edit101 as a double


% --- Executes during object creation, after setting all properties.
function edit101_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit101 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox64.
function checkbox64_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox64



function edit102_Callback(hObject, eventdata, handles)
% hObject    handle to edit102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit102 as text
%        str2double(get(hObject,'String')) returns contents of edit102 as a double


% --- Executes during object creation, after setting all properties.
function edit102_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox65.
function checkbox65_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox65



function edit103_Callback(hObject, eventdata, handles)
% hObject    handle to edit103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit103 as text
%        str2double(get(hObject,'String')) returns contents of edit103 as a double


% --- Executes during object creation, after setting all properties.
function edit103_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox66.
function checkbox66_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox66


% --- Executes on selection change in popupmenu118.
function popupmenu118_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu118 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu118 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu118


% --- Executes during object creation, after setting all properties.
function popupmenu118_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu118 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu119.
function popupmenu119_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu119 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu119 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu119


% --- Executes during object creation, after setting all properties.
function popupmenu119_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu119 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu120.
function popupmenu120_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu120 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu120 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu120


% --- Executes during object creation, after setting all properties.
function popupmenu120_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu120 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu121.
function popupmenu121_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu121 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu121 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu121


% --- Executes during object creation, after setting all properties.
function popupmenu121_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu121 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu122.
function popupmenu122_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu122 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu122


% --- Executes during object creation, after setting all properties.
function popupmenu122_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu123.
function popupmenu123_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu123 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu123 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu123


% --- Executes during object creation, after setting all properties.
function popupmenu123_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu123 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in txtStartLoc_Hybrid.
function popStartLoc_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtStartLoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns txtStartLoc_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from txtStartLoc_Hybrid


% --- Executes during object creation, after setting all properties.
function txtStartLoc_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtStartLoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in txtFinalLoc_Hybrid.
function popFinalLoc_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtFinalLoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns txtFinalLoc_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from txtFinalLoc_Hybrid


% --- Executes during object creation, after setting all properties.
function txtFinalLoc_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFinalLoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popULoc_Hybrid.
function popULoc_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popULoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popULoc_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popULoc_Hybrid
contents = cellstr(get(hObject,'String'));
loc = contents{get(hObject,'Value')};
set(handles.pbU_Hybrid,'String',['Compose Input for Loc ',num2str(loc)])

% --- Executes during object creation, after setting all properties.
function popULoc_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popULoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popMethod_Hybrid.
function popMethod_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popMethod_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject, 'String'));
method = contents{get(hObject,'Value')};

if strcmp(method, 'polytope') || strcmp(method,'conZonotope') || strcmp(method,'zonoGirard') || strcmp(method,'nondetGuard') 
    set(handles.pb_info_enclose_Hybrid, 'Enable', 'on')
    set(handles.checkbox_Hybrid, 'Enable', 'on')
    set(handles.checkpca_Hybrid, 'Enable', 'on')
    set(handles.checkflow_Hybrid, 'Enable', 'on')   
else
    set(handles.pb_info_enclose_Hybrid, 'Enable', 'off')
    set(handles.checkbox_Hybrid, 'Enable', 'off')
    set(handles.checkpca_Hybrid, 'Enable', 'off')
    set(handles.checkflow_Hybrid, 'Enable', 'off')
end

if strcmp(method, 'hyperplaneMap')|| strcmp(method,'conZonotope')
    set(handles.txtGO_Hybrid, 'Enable', 'on')
    set(handles.GOtxt_Hybrid, 'Enable', 'on')
else
    set(handles.txtGO_Hybrid, 'Enable', 'off')
    set(handles.GOtxt_Hybrid, 'Enable', 'off')
end
    
% Hints: contents = cellstr(get(hObject,'String')) returns popMethod_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popMethod_Hybrid


% --- Executes during object creation, after setting all properties.
function popMethod_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popMethod_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu132.
function popupmenu132_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu132 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu132 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu132


% --- Executes during object creation, after setting all properties.
function popupmenu132_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu132 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Hybrid System -------------------------------------------------------

% --- Callback Functions - Hybrid -----------------------------------------


function popStartLoc_hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtstartloc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtstartloc_hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtstartloc_hybrid as a double
start_loc = get(hObject, 'String');
start_loc = str2double(start_loc);
handles.start_loc = start_loc;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popStartLoc_hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtstartloc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popEndLoc_hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popEndLoc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popEndLoc_hybrid as text
%        str2double(get(hObject,'String')) returns contents of popEndLoc_hybrid as a double
end_loc = get(hObject, 'String');
end_loc = str2double(end_loc);
handles.end_loc = end_loc;
guidata(hObject, handles)



% --- Executes during object creation, after setting all properties.
function popEndLoc_hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popEndLoc_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLocR0_hybrid.
function popLocR0_hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popLocR0_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popLocR0_hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLocR0_hybrid

% --- Executes during object creation, after setting all properties.
function popLocR0_hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLocR0_hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popDim1Plot_Hybrid.
function popDim1Plot_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popDim1Plot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popDim1Plot_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popDim1Plot_Hybrid


% --- Executes during object creation, after setting all properties.
function popDim1Plot_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDim1Plot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popDim2Plot_Hybrid.
function popDim2Plot_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popDim2Plot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popDim2Plot_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popDim2Plot_Hybrid


% --- Executes during object creation, after setting all properties.
function popDim2Plot_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDim2Plot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushAddPlot_Hybrid.
function pushAddPlot_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pushAddPlot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listPlots_Hybrid, 'Value', 1);

contents_dim1 = cellstr(get(handles.popDim1Plot_Hybrid, 'String'));
dim1 = contents_dim1{get(handles.popDim1Plot_Hybrid, 'Value')};
dim1 = dim1(~isspace(dim1));

contents_dim2 = cellstr(get(handles.popDim2Plot_Hybrid, 'String'));
dim2 = contents_dim2{get(handles.popDim2Plot_Hybrid, 'Value')};
dim2 = dim2(~isspace(dim2));

if isempty(dim1) || isempty(dim2)
    return
end

if dim1 == 'time'
    dim = sprintf('[%s, %i]', dim1, str2double(dim2));
else
    dim = sprintf('[%i, %i]', str2double(dim1), str2double(dim2));
end
handles.hybrid_plots = [handles.hybrid_plots, dim];
guidata(hObject, handles);

set(handles.listPlots_Hybrid, 'String', handles.hybrid_plots);

function listPlots_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to listPlots_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listPlots_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listPlots_Hybrid


% --- Executes during object creation, after setting all properties.
function  listPlots_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listPlots_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkReachPlot_Hybrid.
function checkReachPlot_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkReachPlot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in checkInitialPlot_Hybrid.
function checkInitialPlot_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkInitialPlot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in checkSimulationPlot_Hybrid.
function checkSimulationPlot_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkSimulationPlot_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushDeletePlot_Nonlinear.
function pushDeletePlot_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pushDeletePlot_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.listPlots_Hybrid,'String'));
value = get(handles.listPlots_Hybrid, 'Value');

if isempty(contents)
    return
end

contents(value,:) = [];

handles.hybrid_plots = contents';
guidata(hObject, handles)
set(handles.listPlots_Hybrid, 'Value', 1);
set(handles.listPlots_Hybrid, 'String', contents);

% --- Executes on selection change in popColorSimulation_Hybrid.
function popColorSimulation_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popColorSimulation_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorSimulation_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorSimulation_Hybrid


% --- Executes during object creation, after setting all properties.
function popColorSimulation_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorSimulation_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLineSimulation_Hybrid.
function popLineSimulation_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popLineSimulation_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popLineSimulation_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLineSimulation_Hybrid


% --- Executes during object creation, after setting all properties.
function popLineSimulation_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLineSimulation_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popColorInitial_Hybrid.
function popColorInitial_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popColorInitial_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorInitial_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorInitial_Hybrid


% --- Executes during object creation, after setting all properties.
function popColorInitial_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorInitial_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popEdgeColorInitial_Hybrid.
function popEdgeColorInitial_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popEdgeColorInitial_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popEdgeColorInitial_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popEdgeColorInitial_Hybrid


% --- Executes during object creation, after setting all properties.
function popEdgeColorInitial_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popEdgeColorInitial_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popColorReach_Hybrid.
function popColorReach_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popColorReach_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popColorReach_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColorReach_Hybrid


% --- Executes during object creation, after setting all properties.
function popColorReach_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColorReach_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popEdgeColorReach_Hybrid.
function popEdgeColorReach_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popEdgeColorReach_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popEdgeColorReach_Hybrid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popEdgeColorReach_Hybrid


% --- Executes during object creation, after setting all properties.
function popEdgeColorReach_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popEdgeColorReach_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in rbSim_Hybrid.
function rbSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to rbSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.rbRSim_Hybrid, 'Value', 0)
    set(handles.rbSimRRT_Hybrid, 'Value', 0)

    set(handles.checkSimulationPlot_Hybrid, 'Enable', 'on')
    set(handles.popColorSimulation_Hybrid, 'Enable', 'on')
    set(handles.popLineSimulation_Hybrid, 'Enable', 'on')
    set(handles.text196, 'Enable', 'on')
    set(handles.text197, 'Enable', 'on')

    set(handles.txtX0_Sim_Hybrid, 'Enable', 'on')
    set(handles.txtU_Sim_Hybrid, 'Enable', 'on')
    set(handles.popX0_Sim_Hybrid, 'Enable', 'on')
    set(handles.popU_Sim_Hybrid, 'Enable', 'on')

    set(handles.txtNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtFV_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtFIV_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtNCI_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.rbYesEPS_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Hybrid, 'Enable', 'off')

else
    
    set(handles.checkSimulationPlot_Hybrid, 'Enable', 'off')
    set(handles.popColorSimulation_Hybrid, 'Enable', 'off')
    set(handles.popLineSimulation_Hybrid, 'Enable', 'off')
    set(handles.text196, 'Enable', 'off')
    set(handles.text197, 'Enable', 'off')

    set(handles.txtX0_Sim_Hybrid, 'Enable', 'off')
    set(handles.txtU_Sim_Hybrid, 'Enable', 'off')
    set(handles.popX0_Sim_Hybrid, 'Enable', 'off')
    set(handles.popU_Sim_Hybrid, 'Enable', 'off')
    
end



function txtX0_Sim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtX0_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x0 = get(hObject, 'String');
handles.x0_hybrid = x0;
guidata(hObject, handles);
set(handles.popX0_Sim_Hybrid, 'Value', 1)

% --- Executes during object creation, after setting all properties.
function txtX0_Sim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtX0_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtU_Sim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtU_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
u = get(hObject, 'String');
handles.u_hybrid = u;
guidata(hObject, handles);
set(handles.popU_Sim_Hybrid, 'Value', 1)

% --- Executes during object creation, after setting all properties.
function txtU_Sim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtU_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbRSim_Hybrid.
function rbRSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to rbRSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.rbSim_Hybrid, 'Value', 0)
    set(handles.rbSimRRT_Hybrid, 'Value', 0)

    set(handles.checkSimulationPlot_Hybrid, 'Enable', 'on')
    set(handles.popColorSimulation_Hybrid, 'Enable', 'on')
    set(handles.popLineSimulation_Hybrid, 'Enable', 'on')
    set(handles.text196, 'Enable', 'on')
    set(handles.text197, 'Enable', 'on')

    set(handles.txtX0_Sim_Hybrid, 'Enable', 'off')
    set(handles.txtU_Sim_Hybrid, 'Enable', 'off')
    set(handles.popX0_Sim_Hybrid, 'Enable', 'off')
    set(handles.popU_Sim_Hybrid, 'Enable', 'off')
    
    set(handles.txtNoP_RSim_Hybrid, 'Enable', 'on')
    set(handles.pbUpNoP_RSim_Hybrid, 'Enable', 'on')
    set(handles.pbDownNoP_RSim_Hybrid, 'Enable', 'on')
    set(handles.txtFV_RSim_Hybrid, 'Enable', 'on')
    set(handles.txtFIV_RSim_Hybrid, 'Enable', 'on')
    set(handles.txtNCI_RSim_Hybrid, 'Enable', 'on')
    set(handles.pbUpNCI_RSim_Hybrid, 'Enable', 'on')
    set(handles.pbDownNCI_RSim_Hybrid, 'Enable', 'on')
    
    set(handles.txtNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.rbYesEPS_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Hybrid, 'Enable', 'off')
    
else
    
    set(handles.checkSimulationPlot_Hybrid, 'Enable', 'off')
    set(handles.popColorSimulation_Hybrid, 'Enable', 'off')
    set(handles.popLineSimulation_Hybrid, 'Enable', 'off')
    set(handles.text196, 'Enable', 'off')
    set(handles.text197, 'Enable', 'off')

    set(handles.txtNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtFV_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtFIV_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtNCI_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Hybrid, 'Enable', 'off')
    
end



function txtNoP_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtNoP_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtNoP_RSim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNoP_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownNoP_RSim_Hybrid.
function pbDownNoP_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNoP_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_RSim_Hybrid, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNoP_RSim_Hybrid, 'String', value);
end

% --- Executes on button press in pbUpNoP_RSim_Hybrid.
function pbUpNoP_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNoP_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_RSim_Hybrid, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNoP_RSim_Hybrid, 'String', value);


function txtFV_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtFV_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));
if value > 1
    value = 1;
elseif value < 0
    value = 0;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtFV_RSim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFV_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFIV_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtFIV_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));
if value > 1
    value = 1;
elseif value < 0
    value = 0;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtFIV_RSim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFIV_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtNCI_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtNCI_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtNCI_RSim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNCI_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbSimRRT_Hybrid.
function rbSimRRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to rbSimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.rbRSim_Hybrid, 'Value', 0)
    set(handles.rbSim_Hybrid, 'Value', 0)

    set(handles.checkSimulationPlot_Hybrid, 'Enable', 'on')
    set(handles.popColorSimulation_Hybrid, 'Enable', 'on')
    set(handles.popLineSimulation_Hybrid, 'Enable', 'on')
    set(handles.text196, 'Enable', 'on')
    set(handles.text197, 'Enable', 'on')

    set(handles.txtX0_Sim_Hybrid, 'Enable', 'off')
    set(handles.txtU_Sim_Hybrid, 'Enable', 'off')
    set(handles.popX0_Sim_Hybrid, 'Enable', 'off')
    set(handles.popU_Sim_Hybrid, 'Enable', 'off')
    set(handles.txtNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbUpNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbDownNoP_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtFV_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtFIV_RSim_Hybrid, 'Enable', 'off')
    set(handles.txtNCI_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbUpNCI_RSim_Hybrid, 'Enable', 'off')
    set(handles.pbDownNCI_RSim_Hybrid, 'Enable', 'off')
    
    set(handles.txtNoP_SimRRT_Hybrid, 'Enable', 'on')
    set(handles.pbUpNoP_SimRRT_Hybrid, 'Enable', 'on')
    set(handles.pbDownNoP_SimRRT_Hybrid, 'Enable', 'on')
    set(handles.rbNoEPS_SimRRT_Hybrid, 'Enable', 'on')
    set(handles.txtSF_SimRRT_Hybrid, 'Enable', 'on')
    
else
    set(handles.checkSimulationPlot_Hybrid, 'Enable', 'off')
    set(handles.popColorSimulation_Hybrid, 'Enable', 'off')
    set(handles.popLineSimulation_Hybrid, 'Enable', 'off')
    set(handles.text196, 'Enable', 'off')
    set(handles.text197, 'Enable', 'off')
    
    set(handles.txtNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.pbUpNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.pbDownNoP_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.rbNoEPS_SimRRT_Hybrid, 'Enable', 'off')
    set(handles.txtSF_SimRRT_Hybrid, 'Enable', 'off')
end



function txtNoP_SimRRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtNoP_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtNoP_SimRRT_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNoP_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownNoP_SimRRT_Hybrid.
function pbDownNoP_SimRRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNoP_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_SimRRT_Hybrid, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNoP_SimRRT_Hybrid, 'String', value);
end

% --- Executes on button press in pbUpNoP_SimRRT_Hybrid.
function pbUpNoP_SimRRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNoP_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNoP_SimRRT_Hybrid, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNoP_SimRRT_Hybrid, 'String', value);

% --- Executes on button press in pbDownNCI_RSim_Hybrid.
function pbDownNCI_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownNCI_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNCI_RSim_Hybrid, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtNCI_RSim_Hybrid, 'String', value);
end


% --- Executes on button press in pbUpNCI_RSim_Hybrid.
function pbUpNCI_RSim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpNCI_RSim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtNCI_RSim_Hybrid, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtNCI_RSim_Hybrid, 'String', value);

% --- Executes on button press in rbYesEPS_SimRRT_Hybrid.
function rbYesEPS_SimRRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to rbYesEPS_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'off')
set(handles.rbNoEPS_SimRRT_Hybrid, 'Enable', 'on')
set(handles.rbNoEPS_SimRRT_Hybrid, 'Value', 0)


% --- Executes on button press in rbNoEPS_SimRRT_Hybrid.
function rbNoEPS_SimRRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to rbNoEPS_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'off')
set(handles.rbYesEPS_SimRRT_Hybrid, 'Enable', 'on')
set(handles.rbYesEPS_SimRRT_Hybrid, 'Value', 0)



function txtSF_SimRRT_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtSF_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject, 'String'));
if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtSF_SimRRT_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSF_SimRRT_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popX0_Sim_Hybrid.
function popX0_Sim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popX0_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
x0 = contents{get(hObject,'Value')};
x0_value = evalin('base', x0);
set(handles.txtX0_Sim_Hybrid,'String', mat2str(x0_value));
handles.x0_hybrid = x0;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popX0_Sim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popX0_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popU_Sim_Hybrid.
function popU_Sim_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to popU_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
u = contents{get(hObject,'Value')};
u_value = evalin('base', u);
set(handles.txtU_Sim_Hybrid,'String', mat2str(u_value));
handles.u_hybrid = u;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popU_Sim_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popU_Sim_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Hybrid.
function checkbox_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Hybrid


% --- Executes on button press in checkpca_Hybrid.
function checkpca_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkpca_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkpca_Hybrid


% --- Executes on button press in checkflow_Hybrid.
function checkflow_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to checkflow_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkflow_Hybrid


% --- Executes on button press in txtLA_Hybrid.
function radiobutton51_Callback(hObject, eventdata, handles)
% hObject    handle to txtLA_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of txtLA_Hybrid


% --- Executes on button press in rbConLin_Hybrid.
function rbConLin_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to rbConLin_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Enable', 'off')
set(handles.txtConLinTO_Hybrid, 'Enable', 'on')
set(handles.pbUpConLinTO_Hybrid, 'Enable', 'on')
set(handles.pbDownConLinTO_Hybrid, 'Enable', 'on')
set(handles.txtIO_Hybrid_lin, 'Enable', 'on')
set(handles.txtEO_Hybrid_lin, 'Enable', 'on')
set(handles.pbUpIO_Hybrid_lin, 'Enable', 'on')
set(handles.pbDownIO_Hybrid_lin, 'Enable', 'on')
set(handles.pbUpEO_Hybrid_lin, 'Enable', 'on')
set(handles.pbDownEO_Hybrid_lin, 'Enable', 'on')
set(handles.pb_infoCPIO_Hybrid_lin, 'Enable', 'on')
set(handles.pb_infoCPEE_Hybrid_lin, 'Enable', 'on')
set(handles.txtIO2_Hybrid_lin, 'Enable', 'on')
set(handles.txtEO2_Hybrid_lin, 'Enable', 'on')

set(handles.txtCP_Hybrid, 'Enable', 'on')
set(handles.txtCP_Hybrid, 'Value', 0)
set(handles.txtConPolTO_Hybrid, 'Enable', 'off')
set(handles.pbUpConPolTO_Hybrid, 'Enable', 'off')
set(handles.pbDownConPolTO_Hybrid, 'Enable', 'off')
set(handles.txtIO_Hybrid, 'Enable', 'off')
set(handles.pbUpIO_Hybrid, 'Enable', 'off')
set(handles.pbDownIO_Hybrid, 'Enable', 'off')
set(handles.txtEO_Hybrid, 'Enable', 'off')
set(handles.pbUpEO_Hybrid, 'Enable', 'off')
set(handles.pbDownEO_Hybrid, 'Enable', 'off')

value = int16(str2double(get(handles.txtConLinTO_Hybrid, 'String')));

if value < 3
    set(handles.txtIO_Hybrid_lin, 'Enable', 'off')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'off')
else
    set(handles.txtIO_Hybrid_lin, 'Enable', 'on')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'on')
end



function txtConLinTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtConLinTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtConLinTO_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtConLinTO_Hybrid as a double
value = int16(str2double(get(hObject, 'String')));

if value < 3
    set(handles.txtIO_Hybrid_lin, 'Enable', 'off')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'off')
else
    set(handles.txtIO_Hybrid_lin, 'Enable', 'on')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'on')
end

if value <= 3
    value = num2str(value);
    set(hObject, 'String', value);
end


% --- Executes during object creation, after setting all properties.
function txtConLinTO_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtConLinTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownConLinTO_Hybrid.
function pbDownConLinTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownConLinTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConLinTO_Hybrid, 'String');
value = str2double(value) - 1;

if value < 3
    set(handles.txtIO_Hybrid_lin, 'Enable', 'off')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'off')
else
    set(handles.txtIO_Hybrid_lin, 'Enable', 'on')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'on')
end

if value >=2
    value = num2str(value);
    set(handles.txtConLinTO_Hybrid, 'String', value);
end

% --- Executes on button press in pbUpConLinTO_Hybrid.
function pbUpConLinTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpConLinTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConLinTO_Hybrid, 'String');
value = str2double(value) + 1;

if value < 3
    set(handles.txtIO_Hybrid_lin, 'Enable', 'off')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'off')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'off')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'off')
else
    set(handles.txtIO_Hybrid_lin, 'Enable', 'on')
    set(handles.txtEO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpIO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownIO_Hybrid_lin, 'Enable', 'on')
    
    set(handles.pbUpEO_Hybrid_lin, 'Enable', 'on')
    set(handles.pbDownEO_Hybrid_lin, 'Enable', 'on')
end

if value <= 3
    value = num2str(value);
    set(handles.txtConLinTO_Hybrid, 'String', value);
end

% --- Executes on button press in txtCP_Hybrid.
function txtCP_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtCP_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Enable', 'off')
set(handles.txtConLinTO_Hybrid, 'Enable', 'off')
set(handles.pbUpConLinTO_Hybrid, 'Enable', 'off')
set(handles.pbDownConLinTO_Hybrid, 'Enable', 'off')
set(handles.txtIO_Hybrid_lin, 'Enable', 'off')
set(handles.txtEO_Hybrid_lin, 'Enable', 'off')
set(handles.pbUpIO_Hybrid_lin, 'Enable', 'off')
set(handles.pbDownIO_Hybrid_lin, 'Enable', 'off')
set(handles.pbUpEO_Hybrid_lin, 'Enable', 'off')
set(handles.pbDownEO_Hybrid_lin, 'Enable', 'off')

set(handles.rbConLin_Hybrid, 'Enable', 'on')
set(handles.rbConLin_Hybrid, 'Value', 0)
set(handles.txtConPolTO_Hybrid, 'Enable', 'on')
set(handles.pbUpConPolTO_Hybrid, 'Enable', 'on')
set(handles.pbDownConPolTO_Hybrid, 'Enable', 'on')
set(handles.txtIO_Hybrid, 'Enable', 'on')
set(handles.pbUpIO_Hybrid, 'Enable', 'on')
set(handles.pbDownIO_Hybrid, 'Enable', 'on')
set(handles.txtEO_Hybrid, 'Enable', 'on')
set(handles.pbUpEO_Hybrid, 'Enable', 'on')
set(handles.pbDownEO_Hybrid, 'Enable', 'on')



function txtConPolTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtConPolTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtConPolTO_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtConPolTO_Hybrid as a double
value = int16(str2double(get(hObject, 'String')));

if value < 3
    value = 3;
end

if value > 4
    value = 4;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtConPolTO_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtConPolTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownConPolTO_Hybrid.
function pbDownConPolTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownConPolTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConPolTO_Hybrid, 'String');
value = str2double(value) - 1;

if value >=3
    value = num2str(value);
    set(handles.txtConPolTO_Hybrid, 'String', value);
end

% --- Executes on button press in pbUpConPolTO_Hybrid.
function pbUpConPolTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpConPolTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtConPolTO_Hybrid, 'String');
value = str2double(value) + 1;

if value <= 4
    value = num2str(value);
    set(handles.txtConPolTO_Hybrid, 'String', value);
end


function txtIO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtIO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIO_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtIO_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtIO_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownIO_Hybrid.
function pbDownIO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownIO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Hybrid, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtIO_Hybrid, 'String', value);
end

% --- Executes on button press in pbUpIO_Hybrid.
function pbUpIO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpIO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Hybrid, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtIO_Hybrid, 'String', value);



function txtEO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtEO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEO_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtEO_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtEO_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownEO_Hybrid.
function pbDownEO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownEO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtEO_Hybrid, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtEO_Hybrid, 'String', value);
end

% --- Executes on button press in pbUpEO_Hybrid.
function pbUpEO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpEO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtEO_Hybrid, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtEO_Hybrid, 'String', value);

% --- Executes on button press in pb_infoCL_Hybrid.
function pb_infoCL_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCL_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_consLinearization.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCP_Hybrid.
function pb_infoCP_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCP_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_consPolynomialization.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCLTO_Hybrid.
function pb_infoCLTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCLTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_tensorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCPTO_Hybrid.
function pb_infoCPTO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPTO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_tensorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCPIO_Hybrid.
function pb_infoCPIO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPIO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_intermediateOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_infoCPEE_Hybrid.
function pb_infoCPEE_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPEE_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_errorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in txtNS_Hybrid.
function radiobutton54_Callback(hObject, eventdata, handles)
% hObject    handle to txtNS_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of txtNS_Hybrid



function txtGO_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtGO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGO_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtGO_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtGO_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGO_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFinalLoc_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtFinalLoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFinalLoc_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtFinalLoc_Hybrid as a double


% --- Executes on button press in pb_web1.
function pb_web1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_web1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
url = ' http://spaceex.imag.fr/download-6' ;
web (url)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pb_web1.
function pb_web1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pb_web1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_web2.
function pb_web2_Callback(hObject, eventdata, handles)
% hObject    handle to pb_web2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
url = ' http://github.com/nikos-kekatos/SL2SX' ;
web (url)


% --- Executes on button press in pb_info_enclose_Hybrid.
function pb_info_enclose_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_info_enclose_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_enclose.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pb_input_set_Linear.
function pb_input_set_Linear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_input_set_Linear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_inputSet.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_input_set_Nonlinear.
function pb_input_set_Nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to pb_input_set_Nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_inputSet.png';
infoBox({[path_im, im_ZO]});
uiwait;

% --- Executes on button press in pb_input_set_Hybrid.
function pb_input_set_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_input_set_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_inputSet.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pb_infoCPIO_Nonlinear_lin.
function pb_infoCPIO_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPIO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_intermediateOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pb_infoCPEE_Nonlinear_lin.
function pb_infoCPEE_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPEE_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_errorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;


function txtIO_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to txtIO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIO_Nonlinear_lin as text
%        str2double(get(hObject,'String')) returns contents of txtIO_Nonlinear_lin as a double
value = int16(str2double(get(hObject, 'String')));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);

% --- Executes during object creation, after setting all properties.
function txtIO_Nonlinear_lin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownIO_Nonlinear_lin.
function pbDownIO_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownIO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Nonlinear_lin, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtIO_Nonlinear_lin, 'String', value);
end


% --- Executes on button press in pbUpIO_Nonlinear_lin.
function pbUpIO_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpIO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Nonlinear_lin, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtIO_Nonlinear_lin, 'String', value);



function txtEO_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to txtEO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEO_Nonlinear_lin as text
%        str2double(get(hObject,'String')) returns contents of txtEO_Nonlinear_lin as a double
value = int16(str2double(get(hObject, 'String')));

if value < 1
    value = 1;
end

value = num2str(value);
set(hObject, 'String', value);


% --- Executes during object creation, after setting all properties.
function txtEO_Nonlinear_lin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownEO_Nonlinear_lin.
function pbDownEO_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownEO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtEO_Nonlinear_lin, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtEO_Nonlinear_lin, 'String', value);
end


% --- Executes on button press in pbUpEO_Nonlinear_lin.
function pbUpEO_Nonlinear_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpEO_Nonlinear_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = get(handles.txtEO_Nonlinear_lin, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtEO_Nonlinear_lin, 'String', value);



function txtIO_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to txtIO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIO_Hybrid_lin as text
%        str2double(get(hObject,'String')) returns contents of txtIO_Hybrid_lin as a double


% --- Executes during object creation, after setting all properties.
function txtIO_Hybrid_lin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownIO_Hybrid_lin.
function pbDownIO_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownIO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Hybrid_lin, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtIO_Hybrid_lin, 'String', value);
end


% --- Executes on button press in pbUpIO_Hybrid_lin.
function pbUpIO_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpIO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtIO_Hybrid_lin, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtIO_Hybrid_lin, 'String', value);


function txtEO_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to txtEO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEO_Hybrid_lin as text
%        str2double(get(hObject,'String')) returns contents of txtEO_Hybrid_lin as a double


% --- Executes during object creation, after setting all properties.
function txtEO_Hybrid_lin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDownEO_Hybrid_lin.
function pbDownEO_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbDownEO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtEO_Hybrid_lin, 'String');
value = str2double(value) - 1;

if value > 0
    value = num2str(value);
    set(handles.txtEO_Hybrid_lin, 'String', value);
end


% --- Executes on button press in pbUpEO_Hybrid_lin.
function pbUpEO_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pbUpEO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtEO_Hybrid_lin, 'String');
value = str2double(value) + 1;
value = num2str(value);
set(handles.txtEO_Hybrid_lin, 'String', value);

% --- Executes on button press in pb_infoCPIO_Hybrid_lin.
function pb_infoCPIO_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPIO_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_intermediateOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;


% --- Executes on button press in pb_infoCPEE_Hybrid_lin.
function pb_infoCPEE_Hybrid_lin_Callback(hObject, eventdata, handles)
% hObject    handle to pb_infoCPEE_Hybrid_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_im= [CORAROOT, filesep, 'app', filesep, 'images', filesep];
im_ZO = 'Info_errorOrder.png';
infoBox({[path_im, im_ZO]});
uiwait;



function txtOptionError_Callback(hObject, eventdata, handles)
% hObject    handle to txtOptionError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOptionError as text
%        str2double(get(hObject,'String')) returns contents of txtOptionError as a double


% --- Executes during object creation, after setting all properties.
function txtOptionError_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOptionError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in txtOptionErrorDown.
function txtOptionErrorDown_Callback(hObject, eventdata, handles)
% hObject    handle to txtOptionErrorDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtOptionError, 'String');
value = str2double(value) - 0.01;

if value > 0
    set(handles.txtOptionError, 'String', value);
end

% --- Executes on button press in txtOptionErrorUp.
function txtOptionErrorUp_Callback(hObject, eventdata, handles)
% hObject    handle to txtOptionErrorUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtOptionError, 'String');
value = str2double(value) + 0.01;
set(handles.txtOptionError, 'String', value);



% --- Executes on button press in refresh.
function refresh_Callback(hObject, eventdata, handles)
% hObject    handle to refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
WS_vars = evalin('base', 'who');
WS_vars = [{''}; WS_vars];

%linear
set(handles.popUTrans_Linear, 'String', WS_vars)
set(handles.popUTrans_Linear, 'Value', 1)

set(handles.popX0_Sim_Linear, 'String', WS_vars)
set(handles.popX0_Sim_Linear, 'Value', 1)
set(handles.popU_Sim_Linear, 'String', WS_vars)
set(handles.popU_Sim_Linear, 'Value', 1)

set(handles.popA_Linear, 'String', WS_vars)
set(handles.popA_Linear, 'Value', 1)
set(handles.popB_Linear, 'String', WS_vars)
set(handles.popB_Linear, 'Value', 1)
set(handles.popC_Linear, 'String', WS_vars)
set(handles.popC_Linear, 'Value', 1)
set(handles.popD_Linear, 'String', WS_vars)
set(handles.popD_Linear, 'Value', 1)
set(handles.popE_Linear, 'String', WS_vars)
set(handles.popE_Linear, 'Value', 1)
set(handles.popF_Linear, 'String', WS_vars)
set(handles.popF_Linear, 'Value', 1)

%nonlinear
set(handles.popUTrans_Nonlinear, 'String', WS_vars)
set(handles.popUTrans_Nonlinear, 'Value', 1)

set(handles.popX0_Sim_Nonlinear, 'String', WS_vars)
set(handles.popX0_Sim_Nonlinear, 'Value', 1)
set(handles.popU_Sim_Nonlinear, 'String', WS_vars)
set(handles.popU_Sim_Nonlinear, 'Value', 1)

%hybrid
set(handles.popX0_Sim_Hybrid, 'String', WS_vars)
set(handles.popX0_Sim_Hybrid, 'Value', 1)
set(handles.popU_Sim_Hybrid, 'String', WS_vars)
set(handles.popU_Sim_Hybrid, 'Value', 1)



function txtOptionError_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtOptionError_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOptionError_Hybrid as text
%        str2double(get(hObject,'String')) returns contents of txtOptionError_Hybrid as a double


% --- Executes during object creation, after setting all properties.
function txtOptionError_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOptionError_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in txtOptionErrorDown_Hybrid.
function txtOptionErrorDown_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtOptionErrorDown_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtOptionError_Hybrid, 'String');
value = str2double(value) - 0.01;

if value > 0
    set(handles.txtOptionError_Hybrid, 'String', value);
end

% --- Executes on button press in txtOptionErrorUp_Hybrid.
function txtOptionErrorUp_Hybrid_Callback(hObject, eventdata, handles)
% hObject    handle to txtOptionErrorUp_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.txtOptionError_Hybrid, 'String');
value = str2double(value) + 0.01;
set(handles.txtOptionError_Hybrid, 'String', value);


% --- Executes during object creation, after setting all properties.
function popStartLoc_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popStartLoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popFinalLoc_Hybrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popFinalLoc_Hybrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function text248_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text248 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
