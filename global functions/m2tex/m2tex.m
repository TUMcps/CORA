function varargout = m2tex(varargin)
% This file "m2tex.m" translates a normal m-file into a LaTeX-file, so that
% the layout of both files will look the same.
% 
% Features:     - recognizes all keywords, strings and comments
%               - recognizes all indents and tabs
%               - recognizes cell titles
%               - writes all recognized objects correct in a tex-file
%               - uses original Matlab colors
%               - the font looks almost the same as in Matlab Editor
%               - tex-file is saved using fontencoding "UTF-8", so that
%                 German Umlauts will be written correct
%               - option for numbered code lines ('num' or 'no_num')
%               - recognizes a linebreak, but only with a leading space
%                 character, i.e. " ..." and only once in a line
%               - option for input- and output-file (*.m respectively 
%                 *.tex-file), to specify as input argument
%               - execution via GUI or command line
%               - option ('bwc') for output in black/white for non-colored code
% 
% Application:  There are SEVERAL WAYS TO EXECUTE the program.
%               1. Run the file by pressing F5 oder click the "Run" button 
%                  in the Matlab Editor, while the m2tex.m file is open or
%                  type "m2tex" in the Command Window.
%                  ==> The program will open the GUI and the rest goes from 
%                  there.
%               2. Call the program from the Command Window or another 
%                  m-file with at least one option. The order of the 
%                  option is not important. For example:
% 
%                  m2tex('num','C:\Programme\MATLAB\work\testfile.m','C:\latex_documents\testfile.tex')
%                  m2tex('testfile.m','testoutputfile.tex','no_num','bwc')
% 
%               - If you insert only a filename, m2tex tries to read the
%                 file from your actual current directory.
%               - If you specify no target directory and filename, m2tex 
%                 takes the same path and filename as stated for the 
%                 source m-file.
%               - If you have not specified the source or the option for 
%                 numbered code lines, the program will still ask for it.
% 
% Input:        One m-file.
% Output:       A tex-file will be written to look like the m-file after it
%               has been set to pdf.
% 
% Integration:  Include the output.tex file via the following commands:
%               - \include{output.tex} or
%               - \input{output.tex}
%               You also have to include the following package into the
%               header of your main LaTeX document:
%               - \usepackage{color}
%               - \usepackage{fancyvrb}
%
% created on:   09.06.2009 by USL with Matlab 7.5.0 (R2007b)
% Version:      2.1
% finished:     17.06.2009 (1.0)
% last change:  29.09.2009
%
% If there are any suggestions, problems or ideas for missing functions or
% how to improve the program,
% please contact me at Matlab Central, where you've got the file...
%
% VERSION History:
%   - 1.0   17.06.09    -> translating m-code to tex
%   - 1.1   19.06.09    -> added option 'num' for numbered code lines
%   - 1.2   26.06.09    -> recognizes now a linebreak, but only with a
%                          leading space character, i.e. " ..."
%   - 2.0   06.08.09    -> added option for source and target path and or
%                          filename
%                       -> added a GUI to the command line application
%                       -> set paths, filenames and other options via GUI
%   - 2.1   29.09.09    -> added option 'bwc' for code in black & white
%
% Comments from the writer of this file: (I'm working on it...)
% - recognizes " ...", but only one in a line without problems and only
% with a leading space character
% - no such string as "$" or "#" is allowed to be in the whole m-file in any way
%       --> $,# = variable, which can be replaced in a line, if nessecary
% - only optimized (numbers+their indents) for ?\normalsize?
%
%
%
%
% % % % m2tex M-file for m2tex.fig
% % % %      m2tex, by itself, creates a new m2tex or raises the existing
% % % %      singleton*.
% % % %
% % % %      H = m2tex returns the handle to a new m2tex or the handle to
% % % %      the existing singleton*.
% % % %
% % % %      m2tex('Property','Value',...) creates a new m2tex using the
% % % %      given property value pairs. Unrecognized properties are passed via
% % % %      varargin to m2tex_OpeningFcn.  This calling syntax produces a
% % % %      warning when there is an existing singleton*.
% % % %
% % % %      m2tex('CALLBACK') and m2tex('CALLBACK',hObject,...) call the
% % % %      local function named CALLBACK in m2tex.M with the given input
% % % %      arguments.
% % % %
% % % %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% % % %      instance to run (singleton)".
% % % %
% % % % See also: GUIDE, GUIDATA, GUIHANDLES
% % % 
% % % % Edit the above text to modify the response to help m2tex
% % % 
% % % % Last Modified by GUIDE v2.5 17-Sep-2009 11:06:26

%% check the input options
if ~isempty(varargin)
    if length(varargin) > 1
        % input has to be a real input option for the CLT
        if ~ishandle(varargin{2})
            % call the the command line tool (CLT)
            input = '''';
            for i = 1:length(varargin)
                input = [input varargin{i},''',''']; %#ok<AGROW>
            end
            input = input(1:end-2);
            eval(['m2tex_CLT(',input,')']);
            % finish the execution of this m-file now
            return
        end
    else
        % call the the command line tool (CLT)
        eval(['m2tex_CLT(''',varargin{1},''')']);
        % finish the execution of this m-file now
        return
    end
else
    eval('m2tex_CLT(''test'')');
end

%% normal gui code starts here
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @m2tex_OpeningFcn, ...
                   'gui_OutputFcn',  @m2tex_OutputFcn, ...
                   'gui_LayoutFcn',  @m2tex_LayoutFcn, ...
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

% --- Executes just before m2tex is made visible.
function m2tex_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
% % % % % if ~isempty(varargin)
% % % % %     condition = varargin{1};
% % % % %     clear varargin
% % % % % else
% % % % %     condition = 'm2tex_CLT.m_un-active';
% % % % % end
% % % % % set(handles.text1,'UserData',condition)
% Choose default command line output for m2tex
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes m2tex wait for user response (see UIRESUME)
% uiwait(handles.figure1);
uiwait
% --- Outputs from this function are returned to the command line.

% --- Executes during object creation, after setting all properties.
function button_source_Callback(hObject, eventdata, handles)
% hObject    handle to button_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('*.m','Select your m-file');
if FileName ~= 0
    set(handles.textfield_source, 'String', [PathName,FileName]);
end
% --- Executes on button press in button_targetfile.
function button_targetfile_Callback(hObject, eventdata, handles)
% hObject    handle to button_targetfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
trail.src = get(handles.textfield_source, 'String');
idx = findstr(trail.src,'\');
if isempty(idx)
    idx = findstr(trail.src,'/');
end
FileName_tex = [trail.src(idx(end)+1:end-2),'.tex'];
[FileName,PathName,FilterIndex] = uiputfile('*.tex','Select the output destination',FileName_tex);
if FileName ~= 0
    if ~strcmp(FileName(end-3:end),'.tex')
        FileName = [FileName,'.tex'];
    end
    set(handles.textfield_targetfile, 'String', PathName);
    set(handles.textfield_targetname, 'String', FileName);
end
% --- Executes on button press in reset_savedpath.
function reset_savedpath_Callback(hObject, eventdata, handles)
% hObject    handle to reset_savedpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist('m2tex_prefs.mat','file') == 2
    load m2tex_prefs.mat
    set(handles.textfield_source, 'String', trail.src);
    set(handles.textfield_targetfile, 'String', trail.targetpath);
    set(handles.textfield_targetname, 'String', [trail.targetname]);
end
% --- Executes on button press in reset_workdir.
function reset_workdir_Callback(hObject, eventdata, handles)
% hObject    handle to reset_workdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
trail.src        = get(handles.textfield_source,'string');
if strcmp(trail.src(end-1:end),'.m')
    idx = findstr(trail.src,'\');
    if isempty(idx)
        idx = findstr(trail.src,'/');
    end
    FileName = [trail.src(idx(end)+1:end-2),'.tex'];
    PathName = trail.src(1:idx(end));
    set(handles.textfield_targetfile, 'String',PathName);
    set(handles.textfield_targetname, 'String',FileName);
end
% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
message_error{1} = 'No m-file selected!';
message_error{2} = 'No tex-file selected!';
message_error{3} = 'No target path selected!';
message_error{4} = '';
trail.src        = get(handles.textfield_source,'string');
trail.targetpath = get(handles.textfield_targetfile,'string');
trail.targetname = get(handles.textfield_targetname,'string');
if get(handles.checkbox_numb_codelines,'Value')
    trail.numb = 'num';
else
    trail.numb = 'no_num';
end
if get(handles.checkbox_bwc,'Value')
    trail.bwc = 'bwc';
else
    trail.bwc = 'color';
end
cond = [strcmp(trail.src(end-1:end),'.m') ...
        strcmp(trail.targetname(end-3:end),'.tex') ...
        ~strcmp(trail.targetpath,'nothing is set')];
if cond(1) == 0, idx.a = 1;else idx.a = 4;end
if cond(2) == 0, idx.b = 2;else idx.b = 4;end
if cond(3) == 0, idx.c = 3;else idx.c = 4;end
if all(cond)
    if get(handles.checkbox_saveoptions,'Value')
        save('m2tex_prefs.mat','trail','-mat')
        disp('Source and Target Path and Name are saved.')
    end
    uiresume;
    disp('m2tex_gui properly executed.')
    % call the CLT
    eval(['m2tex_CLT(''',trail.numb,''',''',trail.src,''',''',...
    [trail.targetpath,trail.targetname],''',''',trail.bwc,''')'])
else
    errordlg(sprintf([message_error{idx.a},'\n',message_error{idx.b} ...
        ,'\n',message_error{idx.c}]),'File Error');
end

%% empty object functions
% --- Executes on button press in checkbox_bwc.
function checkbox_bwc_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_bwc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of checkbox_bwc
function checkbox_saveoptions_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_saveoptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of checkbox_saveoptions
function checkbox_numb_codelines_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_numb_codelines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of
% checkbox_numb_codelines

% --- Executes when selected object is changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in save_options.
function save_options_Callback(hObject, eventdata, handles)
% hObject    handle to save_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of save_options
function textfield_targetname_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_targetname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of textfield_targetname as text
%        str2double(get(hObject,'String')) returns contents of textfield_targetname as a double
% --- Executes during object creation, after setting all properties.
function textfield_targetname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_targetname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function textfield_targetfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_targetfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function textfield_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function textfield_targetfile_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_targetfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of textfield_targetfile as text
%        str2double(get(hObject,'String')) returns contents of
%        textfield_targetfile as a double
% --- Executes during object creation, after setting all properties.
% --- Executes on button press in button_source.
function textfield_source_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of textfield_source as text
%        str2double(get(hObject,'String')) returns contents of textfield_source as a double
function varargout = m2tex_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Get default command line output from handles structure
% varargout{1} = handles.output;

%% Graphical Part
% --- Creates and returns a handle to the GUI figure. 
function h1 = m2tex_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', 2, ...
    'pushbutton', 11, ...
    'text', 19, ...
    'edit', 7, ...
    'frame', 4, ...
    'radiobutton', 19, ...
    'uipanel', 13, ...
    'checkbox', 3), ...
    'override', 1, ...
    'release', 13, ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', 'C:\Dokumente und Einstellungen\Admin\Eigene Dateien\My Dropbox\A_m-files\mcode2tex\Studie\m2tex.m', ...
    'lastFilename', 'C:\Dokumente und Einstellungen\Admin\Eigene Dateien\My Dropbox\A_m-files\mcode2tex\m2tex_gui.fig');
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = [];
appdata.initTags = struct(...
    'handle', [], ...
    'tag', 'figure1');

h1 = figure(...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.925490196078431 0.913725490196078 0.847058823529412],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'DockControls','off',...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','m2tex_gui',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[14.840602904 20.974356673],...
'PaperType','A5',...
'Position',[544 485 500 325],...
'Resize','off',...
'UseHG2','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[],...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'calculate';

h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','m2tex(''calculate_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'ListboxTop',0,...
'Position',[69.8 1.30769230769231 26.2 5.92307692307692],...
'String','Run',...
'Tag','calculate',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel1';

h3 = uipanel(...
'Parent',h1,...
'Units','characters',...
'Title','Path Preferences',...
'Tag','uipanel1',...
'UserData',[],...
'Clipping','on',...
'Position',[3.2 8.69230769230769 92.8 11.3076923076923],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text1';

h4 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[2 8.38461538461539 30.6 1.30769230769231],...
'String','Source Path (m-file)',...
'Style','text',...
'Tag','text1',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text2';

h5 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[2 4.84615384615385 26.6 1.15384615384615],...
'String','Target Path (tex-file)',...
'Style','text',...
'Tag','text2',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'textfield_source';

h6 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','m2tex(''textfield_source_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[2 7.15384615384615 88.6 1.15384615384615],...
'String','nothing is set',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'm2tex(''textfield_source_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','textfield_source',...
'UserData',[]);

appdata = [];
appdata.lastValidTag = 'textfield_targetfile';

h7 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','m2tex(''textfield_targetfile_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[2 3.61538461538462 88.6 1.15384615384615],...
'String','nothing is set',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'm2tex(''textfield_targetfile_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','textfield_targetfile',...
'UserData',[]);

appdata = [];
appdata.lastValidTag = 'button_source';

h8 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'Callback','m2tex(''button_source_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[75.8333333333335 8.5576923076923 15 1.33333333333333],...
'String','Browse',...
'Tag','button_source',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'button_targetfile';

h9 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'Callback','m2tex(''button_targetfile_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[75.8333333333335 4.89102564102564 15 1.33333333333333],...
'String','Browse',...
'Tag','button_targetfile',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'textfield_targetname';

h10 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','m2tex(''textfield_targetname_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[19 1.30769230769231 35.6 1.15384615384615],...
'String','nothing is set',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'm2tex(''textfield_targetname_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','textfield_targetname',...
'UserData',[]);

appdata = [];
appdata.lastValidTag = 'text15';

h11 = uicontrol(...
'Parent',h3,...
'Units','characters',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[2 1.3525641025641 16.2 1.15384615384615],...
'String','Target Filename',...
'Style','text',...
'Tag','text15',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'unitgroup';

h12 = uibuttongroup(...
'Parent',h1,...
'Units','characters',...
'Title','Options',...
'Tag','unitgroup',...
'Clipping','on',...
'Position',[3.16666666666667 1.1474358974359 28.6666666666667 6.75],...
'SelectedObject',[],...
'SelectionChangeFcn','test_gui_m2tex(''unitgroup_SelectionChangeFcn'',gcbo,[],guidata(gcbo))',...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'checkbox_numb_codelines';

h13 = uicontrol(...
'Parent',h12,...
'Units','characters',...
'Callback','m2tex(''checkbox_numb_codelines_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[1.83333333333333 2.10897435897436 25 1.91666666666667],...
'String','Numbered Code Lines',...
'Style','checkbox',...
'Value',1,...
'Tag','checkbox_numb_codelines',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'checkbox_saveoptions';

h14 = uicontrol(...
'Parent',h12,...
'Units','characters',...
'Callback','m2tex(''checkbox_saveoptions_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[1.83333333333333 0.442307692307692 21.6666666666667 1.75],...
'String','Save Paths',...
'Style','checkbox',...
'Value',1,...
'Tag','checkbox_saveoptions',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'checkbox_bwc';

h15 = uicontrol(...
'Parent',h12,...
'Units','characters',...
'Callback','m2tex(''checkbox_bwc_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[1.86666666666667 3.75641025641026 25 1.91666666666667],...
'String','All in Black - No Color',...
'Style','checkbox',...
'Tag','checkbox_bwc',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel12';

h16 = uipanel(...
'Parent',h1,...
'Units','characters',...
'Title','Restore Target Path & Name to',...
'Tag','uipanel12',...
'Clipping','on',...
'Position',[33.2 1.15384615384615 34.8 6.84615384615385],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'reset_savedpath';

h17 = uicontrol(...
'Parent',h16,...
'Units','characters',...
'Callback','m2tex(''reset_savedpath_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'ListboxTop',0,...
'Position',[7.43333333333334 3.42948717948718 18.5 1.75],...
'String','Saved Ones',...
'Tag','reset_savedpath',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'reset_workdir';

h18 = uicontrol(...
'Parent',h16,...
'Units','characters',...
'Callback','m2tex(''reset_workdir_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'ListboxTop',0,...
'Position',[7.43333333333334 0.929487179487177 18.5 1.75],...
'String','Source',...
'Tag','reset_workdir',...
'UserData',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel11';

h19 = uipanel(...
'Parent',h1,...
'Units','characters',...
'Title',blanks(0),...
'Tag','uipanel11',...
'Clipping','on',...
'Position',[15 20.8461538461538 70 3.92307692307692],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text18';

h20 = uicontrol(...
'Parent',h19,...
'Units','characters',...
'HorizontalAlignment','left',...
'Position',[8.2 0.282051282051283 26.8333333333333 0.916666666666667],...
'String','Version: 2.1',...
'Style','text',...
'Tag','text18',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text17';

h21 = uicontrol(...
'Parent',h19,...
'Units','characters',...
'HorizontalAlignment','right',...
'Position',[48.2 0.198717948717949 14.1666666666667 1],...
'String',' 2009 USL',...
'Style','text',...
'Tag','text17',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text16';

h22 = uicontrol(...
'Parent',h19,...
'Units','characters',...
'FontSize',19.5,...
'FontWeight','bold',...
'Position',[0.4 1.07692307692308 68.6 2.61538461538462],...
'String','"m2tex.m" Preferences',...
'Style','text',...
'Tag','text16',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );


hsingleton = h1;

% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   eval(createfcn);
end

% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
    'gui_Singleton'
    'gui_OpeningFcn'
    'gui_OutputFcn'
    'gui_LayoutFcn'
    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('MATLAB:gui_mainfcn:FieldNotFound', 'Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % M2TEX
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % M2TEX(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallbak(gui_State, varargin{:})
    % M2TEX('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % M2TEX(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = true;
end

if ~gui_Create
    % In design time, we need to mark all components possibly created in
    % the coming callback evaluation as non-serializable. This way, they
    % will not be brought into GUIDE and not be saved in the figure file
    % when running/saving the GUI from GUIDE.
    designEval = false;
    if (numargin>1 && ishghandle(varargin{2}))
        fig = varargin{2};
        while ~isempty(fig) && ~isa(handle(fig),'figure')
            fig = get(fig,'parent');
        end
        
        designEval = isappdata(0,'CreatingGUIDEFigure') || isprop(fig,'__GUIDEFigure');
    end
        
    if designEval
        beforeChildren = findall(fig);
    end
    
    % evaluate the callback now
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else       
        feval(varargin{:});
    end
    
    % Set serializable of objects created in the above callback to off in
    % design time. Need to check whether figure handle is still valid in
    % case the figure is deleted during the callback dispatching.
    if designEval && ishandle(fig)
        set(setdiff(findall(fig),beforeChildren), 'Serializable','off');
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end

    % Check user passing 'visible' P/V pair first so that its value can be
    % used by oepnfig to prevent flickering
    gui_Visible = 'auto';
    gui_VisibleInput = '';
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        % Recognize 'visible' P/V pair
        len1 = min(length('visible'),length(varargin{index}));
        len2 = min(length('off'),length(varargin{index+1}));
        if ischar(varargin{index+1}) && strncmpi(varargin{index},'visible',len1) && len2 > 1
            if strncmpi(varargin{index+1},'off',len2)
                gui_Visible = 'invisible';
                gui_VisibleInput = 'off';
            elseif strncmpi(varargin{index+1},'on',len2)
                gui_Visible = 'visible';
                gui_VisibleInput = 'on';
            end
        end
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.

    
    % Do feval on layout code in m-file if it exists
    gui_Exported = ~isempty(gui_State.gui_LayoutFcn);
    % this application data is used to indicate the running mode of a GUIDE
    % GUI to distinguish it from the design mode of the GUI in GUIDE. it is
    % only used by actxproxy at this time.   
    setappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]),1);
    if gui_Exported
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % make figure invisible here so that the visibility of figure is
        % consistent in OpeningFcn in the exported GUI case
        if isempty(gui_VisibleInput)
            gui_VisibleInput = get(gui_hFigure,'Visible');
        end
        set(gui_hFigure,'Visible','off')

        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen');
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        end
    end
    if isappdata(0, genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]))
        rmappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]));
    end

    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    % Singleton setting in the GUI M-file takes priority if different
    gui_Options.singleton = gui_State.gui_Singleton;

    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end

    % Apply input P/V pairs other than 'visible'
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        len1 = min(length('visible'),length(varargin{index}));
        if ~strncmpi(varargin{index},'visible',len1)
            try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
        end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    if isscalar(gui_hFigure) && ishandle(gui_hFigure)
        % Handle the default callbacks of predefined toolbar tools in this
        % GUI, if any
        guidemfile('restoreToolbarToolPredefinedCallback',gui_hFigure); 
        
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        % Call openfig again to pick up the saved visibility or apply the
        % one passed in from the P/V pairs
        if ~gui_Exported
            gui_hFigure = local_openfig(gui_State.gui_Name, 'reuse',gui_Visible);
        elseif ~isempty(gui_VisibleInput)
            set(gui_hFigure,'Visible',gui_VisibleInput);
        end
        if strcmpi(get(gui_hFigure, 'Visible'), 'on')
            figure(gui_hFigure);
            
            if gui_Options.singleton
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        if isappdata(gui_hFigure,'InGUIInitialization')
            rmappdata(gui_hFigure,'InGUIInitialization');
        end

        % If handle visibility is set to 'callback', turn it on until
        % finished with OutputFcn
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end

    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end

    if isscalar(gui_hFigure) && ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end

function gui_hFigure = local_openfig(name, singleton, visible)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
if nargin('openfig') == 2
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
else
    gui_hFigure = openfig(name, singleton, visible);
end

function result = local_isInvokeActiveXCallback(gui_State, varargin)

try
    result = ispc && iscom(varargin{1}) ...
             && isequal(varargin{1},gcbo);
catch
    result = false;
end

function result = local_isInvokeHGCallbak(gui_State, varargin)

try
    fhandle = functions(gui_State.gui_Callback);
    result = ~isempty(findstr(gui_State.gui_Name,fhandle.file)) || ...
             (ischar(varargin{1}) ...
             && isequal(ishandle(varargin{2}), 1) ...
             && ~isempty(strfind(varargin{1},[get(varargin{2}, 'Tag'), '_'])));
catch
    result = false;
end

%% Command Line Tool (CLT)
function m2tex_CLT(varargin)

% clc
dbstop if error % dbclear if error
%% check the input options
[option] = input_options(varargin);
%% after clarification of the options and paths...
if option.FilterIndex == 0
    disp('The process was aborted from the user.')
elseif option.FilterIndex == 1
    fid = fopen([option.srcpath,option.srcname]);
    if fid ~= -1
        tex_body = '';
        %--- for debug purpose only ---------------------------------------
%         cnt_skip_lines = 58;    % ? % how many lines to skip
%         for t = 1 : cnt_skip_lines
%             line = fgetl(fid);      %#ok<NASGU> % skip lines
%         end
        %------------------------------------------------------------------
        l_anz = 0;lib_ex = [0 0];
        while feof(fid) == 0        % read until end of file is reached
            % as of this loop, read line after line
            line  = fgetl(fid);
            % check existence of the 3 'colors'
            [idx] = check_colors(line);
            % verify validation
            [idx,lib_ex] = check_points(line,idx,lib_ex);
            % execute conversion
            [line] = write_tex(line,idx,option);
            % numbered code lines
            l_anz = l_anz + 1;
            [prae_line,dist] = numb_line(option.numb,l_anz);
            % add line to already edited text
            tex_body = sprintf('\r %s',tex_body,[prae_line,line(1),dist,line(2:end)]);
            tex_body = tex_body(2:end);
        end
        status1 = fclose(fid);       % close file
        [header,footer] = write_header;
        tex_ges = sprintf('%s \n',header,['\noindent',tex_body],footer);
        tex_ges = tex_ges(1:end-2);
        % % %         clc
        % % %         tex_ges
        % save finished code as tex-file
        fid = fopen([option.tgtpath,option.tgtname],'w','native','UTF-8');
        if fid ~= -1
            count = fprintf(fid,'%s',tex_ges');
            disp('Program successful executed.')
        else
            disp([option.tgtpath,option.tgtname])
            disp('File couldn''t be opened with "fopen".')
        end
        status2 = fclose(fid);       % close file
    else
        disp([option.srcpath,option.srcname])
        disp('File couldn''t be opened with "fopen".')
    end
    dbclear if error
end
function [idx] = check_colors(line)
% This subroutine looks and finds all possible hints, which could cause a
% change of color in the m-code.

% initialisation of the struct variables
idx.com     = [];   % location of comment
idx.com2    = [];   % location of double-comment aka cell titles
idx.str     = [];   % location of strings
idx.str_d   = [];   % length
idx.key     = [];   % location of keywords
idx.key_d   = [];   % length
idx.fsp     = [];   % location of space characters
idx.fsp_d   = [];   % length
idx.fs2     = [];   % location of "second" word strings
idx.lib     = [];   % location of a linebreak

% find all percentage signs
idx.com = findstr('%',line);
if length(line) >= 2, idx.com2= findstr('%%',line(1:2));end
% find all apostrophe signs
idx.str = findstr('''',line);
% if there are only space characters, make line really empty
if all(isspace(line))
    line = [];
end
if ~isempty(line)   % only if line is not empty
    word = textscan(line,'%s');
    wort = unique(word{:});
    word = word{:};
    % find all 'keywords'
    for i=1:length(wort)
        if iskeyword(wort{i})
            idx.key(length(idx.key)+(1:length(findstr(wort{i},line))))     = findstr(wort{i},line);
            idx.key_d(length(idx.key_d)+(1:length(findstr(wort{i},line)))) = length(wort{i});
            if strcmp(wort{i},'elseif')
                line = strrep(line,'elseif','abcdef');
            end
        end
    end
    % find "second" words
    if length(line) > 3 && length(word) >= 2
        if ~iskeyword(word{1}) && ~strcmp(word{2}(1),'=')
            index = findstr(word{2},line);
%             idx.fs2(length(idx.fs2)+(1:length(index(1))))     = index(1);
            idx.fs2 = index(1);
        end
    end
    % find a linebreak " ..."
    idx.lib = findstr(' ...',line);
end
% find all space characters
idx.fsp = findstr(' ',line);
idx.fsp_d = diff(idx.fsp);

if ~isempty(idx.lib)
    % treat everything after linebreak like a comment
    idx.com = [idx.com idx.lib+4];
    % treat the " ..." like a keyword
    idx.key = [idx.key idx.lib];
    idx.key_d = [idx.key_d idx.lib./idx.lib*4];
end
function [idx,lib_ex] = check_points(line,idx,lib_ex)
% This subroutine checks all found hints of the subroutine
% "check_colors", wether they are valid.
%% if double-comments exist
if ~isempty(idx.com2)
    idx.com     = [];
    idx.str     = [];
    idx.key     = [];
    idx.fs2     = [];
    idx.lib     = [];
end
%% Group I
if ~isempty(idx.str) && ~isempty(idx.com)
    % clears "%" inside a string existing of pairs of apostrophs
    for i = 1:2:2*fix(length(idx.str)/2)
        for j = length(idx.com):-1:1
            index = find((idx.com(j) > idx.str(i) & idx.com(j) < idx.str(i+1)));
            idx.com(index) = [];
        end
    end
    % of valid "%" only the first will be needed
    if ~isempty(idx.com)
        idx.com = idx.com(1);
        idx.str(find(idx.com < idx.str)) = [];
        % clear all keywords after first "%"
        index = find(idx.com < idx.key);
        idx.key(index) = [];
        idx.key_d(index) = [];
    end
    % clear keywords, which are inside a string
    for i = 1:2:2*fix(length(idx.str)/2)
        index = find((idx.key > idx.str(i) & idx.key < idx.str(i+1)));
        idx.key(index) = [];
        idx.key_d(index) = [];
    end
end
%% Group II
if isempty(idx.str) && ~isempty(idx.com) && ~isempty(idx.key)
    % remember only the first "%" sign
    idx.com = min(idx.com);
    % clear all keywords after first "%"
    index = find(idx.com < idx.key);
    idx.key(index) = [];
    idx.key_d(index) = [];
end
%% Group IIa
if ~isempty(idx.str)
    % clear strings, which are really a transpose
    index_es = findstr(line,'=');
    % if one cond. is 1, then "'" is beginning of a string
    % if you find between "=" and "'" a left bracket
    cond1 = all(isspace(line(index_es+1:idx.str(1)-1)));
    cond2 = any(findstr(line(index_es+1:idx.str(1)-1),'('));
    cond3 = any(findstr(line(index_es+1:idx.str(1)-1),'{'));
    cond4 = any(findstr(line(index_es+1:idx.str(1)-1),'['));
    % if any condition is true, then delete idx.str(i)
    if ~any([cond1 cond2 cond3 cond4])
        idx.str(1) = [];
    end
    % if you find left next to "'" a right bracket, number or a letter
    % iteration for all "'", from the beginning
    % 1. delete all apostrophs in line2, which are not strings
    line2 = line;
    apost = strfind(line2,'''');
    for i = 1:length(idx.str)
        apost(find(apost == idx.str(i))) = [];
    end
    line2(apost) = '$';
    % 2. for all apostrophs left, make the following test
    idxstr_length = length(idx.str);
    for i = 1:idxstr_length
        if ~isempty(idx.str)
            if idx.str(1) > 1
                [token, remain] = strtok(line2,'''');
                cond2(1) = strcmp(token(end),')');
                cond2(2) = strcmp(token(end),']');
                cond2(3) = strcmp(token(end),'}');
                cond2(4) = ~isempty(str2num(token(end)));
                cond2(5) = isletter(token(end));
                % if any condition is true, then delete idx.str(i)
                if any(cond2)
                    line2(idx.str(1)) = '$';
                    idx.str(1) = [];
                    idxstr_length = length(idx.str);
                end
            end
        end
    end
end
%% Group III
if ~isempty(idx.str) && isempty(idx.com) && ~isempty(idx.key)
    % clear keywords, which are inside a string
    for i = 1:2:2*fix(length(idx.str)/2)
        index = find((idx.key > idx.str(i) & idx.key < idx.str(i+1)));
        idx.key(index) = [];
        idx.key_d(index) = [];
    end
end
%% last Group: "second" words and other stuff...
% clear the keyword "end", if it's right after ":", for example: "1:end"
if ~isempty(idx.key)
    index = idx.key(find(idx.key-1 >0 ))-1;
    index = index(line(index) == ':')+1;
    if ~isempty(index)
        index = find(index == idx.key);
        idx.key(index) = [];
        idx.key_d(index) = [];
    end
end
% this part does what? again?
if length(idx.str) >= 4
    idx.str([find(diff(idx.str) == 1) find(diff(idx.str) == 1)+1]) = [];
end
% ?? was war das hier? war fr abstnde gedacht, nicht mehr ntig!!!
if ~isempty(idx.fsp)
    index = find(diff(idx.fsp)~=1);
    idx.fsp([index length(idx.fsp)])   = [];
    idx.fsp_d(index) = [];
    if isempty(idx.fsp)
        idx.fsp = [];
        idx.fsp_d = [];
    end
end
% clear "second" words on the right side of "%"
if ~isempty(idx.com) && ~isempty(idx.fs2)
    idx.fs2(find(idx.com(1) <= idx.fs2)) = [];
end
if ~isempty(idx.key) && ~isempty(idx.fs2)
    % clear keywords, which are behind a "second" word, except a linebreak
    index = find(idx.key >= idx.fs2 & ~strcmp(' ...',line(idx.key:idx.key+3)));
    idx.key(index)   = [];
    idx.key_d(index) = [];
end
% clear "second" words, which are inside a string
if ~isempty(idx.str) && ~isempty(idx.fs2)
    % fs2's rausschmeissen, die innerhalb von strings sind % NOCH NTIG ??
    for i = 1:2:2*fix(length(idx.str)/2)
        index = find(idx.fs2 > idx.str(i) & idx.fs2 < idx.str(i+1));
        idx.fs2(index) = [];
    end
end
% if linebreak exists in previous line
if (lib_ex(1) == 1) && (lib_ex(2) == 1)
    idx.fs2 = 1;
elseif (lib_ex(1) == 1) && (lib_ex(2) == 0)
    idx.fs2 = [];
end
% set conditions for next line if linebreak exists
lib_ex = [0 0];
if ~isempty(idx.key) && ~isempty(idx.lib)
    if any(idx.key == idx.lib)
        lib_ex(1) = 1; % linebreak  exists in this line
    end
    if ~isempty(idx.fs2)
        lib_ex(2) = 1; % "2nd" word exists
    end
end
function [line] = write_tex(line,idx,option)
%% color definitions for black/white or colored code
str_col.key_front = '$\color{mblue}$';    % blue
str_col.key_rear  = '$\color{black}$';    % black
if strcmp(option.bwc,'bwc')
    str_col.key_front = '$#';             % change to boldface
    str_col.key_rear  = '#$';             % return to not boldface
end
%% read the whole strings, which are colored
text.com = [];
text.str = [];
text.key = [];
text.fs2 = [];
if ~isempty(idx.com)
    text.com{1} = line(idx.com:end);
end
if ~isempty(idx.str)
    for i = 1:2:2*fix(length(idx.str)/2)
        text.str{i} = line(idx.str(i):idx.str(i+1));
    end
end
if ~isempty(idx.key)
    for i = 1:length(idx.key)
        text.key{i} = line(idx.key(i):idx.key(i)+idx.key_d(i)-1);
    end
end
if ~isempty(idx.fs2)
    text.fs2{1} = line(idx.fs2:end);
    if ~isempty(idx.com)
        text.fs2{1} = line(idx.fs2:idx.com-1);
    end
    if ~isempty(idx.lib)
        text.fs2{1} = line(idx.fs2:idx.lib-1);
    end
end
%% replace text-strings in line
% DOUBLE_COMMENT
if ~isempty(idx.com2)
    line = ['$\color{mgreen}#',line,'#\color{black}$'];
end
% KEYWORDS
if ~isempty(text.key)
    keys = sort(idx.key);
    for i = 1:length(keys)
        key_text_sorted{i} = text.key{find(keys(i) == idx.key)};
    end
    text.key = key_text_sorted;
    for i = length(text.key):-1:1
        n = [];j = length(keys);
        while isempty(n)
            n = find(findstr(text.key{i},line) == keys(j));
            j = max(j - 1,1);
        end
        text_rep = [str_col.key_front,text.key{i},str_col.key_rear];
        if strcmp(' ...',text.key{i})
            % regexprep won't work in case of looking for ' ...', why???
            % so only one ' ...' in line is allowed
            line = strrep(line,text.key{i},text_rep);
        else
            line = regexprep(line,text.key{i},text_rep,n);
        end
    end
    if ~isempty(findstr('$color',line))
        line = strrep(line,'$color','$\color');
    end
end
% COMMENTS
if ~isempty(text.com)
    text_rep = ['$\color{mgreen}$',text.com{1},'$\color{black}$'];
    line = strrep(line,text.com{1},text_rep);
end
% STRINGS
if ~isempty(text.str)
    for i = 1:length(text.str)
        text_rep = ['$\color{mred}$',text.str{i},'$\color{black}$'];
        if ~isempty(text.str{i})
            line = strrep(line,text.str{i},text_rep);
        end
    end
end
% "SECOND" WORDS, are strings too
if ~isempty(text.fs2)
    text_rep = ['$\color{mred}$',text.fs2{1},'$\color{black}$'];
    line = strrep(line,text.fs2{1},text_rep);
end
line = ['$',line,'$\\'];
%% option for black/white code activated
if strcmp(option.bwc,'bwc')
    line = strrep(line,'\color{mred}','\color{mdarkgrey}');
    line = strrep(line,'\color{mgreen}','\color{mgrey}');
end
function [header,footer] = write_header
header = sprintf('%s \n'...
    ,'% This file was automatically created from the m-file'...
    ,'% "m2tex.m" written by USL.'...
    ,'% The fontencoding in this file is UTF-8.'...
    ,'% '...
    ,'% You will need to include the following two packages in'...
    ,'% your LaTeX-Main-File.'...
    ,'% '...
    ,'% \usepackage{color}'...
    ,'% \usepackage{fancyvrb}'...
    ,'% '...
    ,'% It is advised to use the following option for Inputenc'...
    ,'% \usepackage[utf8]{inputenc}'...
    ,'% '...
    ,' '...
    ,'% definition of matlab colors:'...
    ,'\definecolor{mblue}{rgb}{0,0,1}'...
    ,'\definecolor{mgreen}{rgb}{0.13333,0.5451,0.13333}'...
    ,'\definecolor{mred}{rgb}{0.62745,0.12549,0.94118}'...
    ,'\definecolor{mgrey}{rgb}{0.5,0.5,0.5}'...
    ,'\definecolor{mdarkgrey}{rgb}{0.25,0.25,0.25}'...
    ,' '...
    ,'\DefineShortVerb[fontfamily=courier,fontseries=m]{\$}'...
    ,'\DefineShortVerb[fontfamily=courier,fontseries=b]{\#}'...
    ,' '...
    );
header = header(1:end-2);

footer = sprintf('%s \n'...
    ,' '...
    ,'\UndefineShortVerb{\$}'...
    ,'\UndefineShortVerb{\#}'...
    );
footer = footer(1:end-2);
function [prae_line,dist_line] = numb_line(numb,l_anz)
if strcmp(numb,'num')
    format bank
    % distances only optimized for numbers in \scriptsize and code in
    % \normalsize
    offset = -3.2; % = 0: nicht eingerckt, sonst code auf line (-3.25)
    if l_anz > 0,   n_dist = ['\hspace*{',num2str(offset+1.6),'em}'];end
    if l_anz > 9,   n_dist = ['\hspace*{',num2str(offset+1.2),'em}'];end
    if l_anz > 99,  n_dist = ['\hspace*{',num2str(offset+0.8),'em}'];end
    if l_anz > 999, n_dist = ['\hspace*{',num2str(offset+0.4),'em}'];end
    if l_anz > 9999,n_dist = ['\hspace*{',num2str(offset+0.0),'em}'];end
    % tiny, scriptsize, footnotesize, small, normalsize
    prae_line = [n_dist,'{\scriptsize ',num2str(l_anz),'}'];
    dist_line = '  ';
    format short
else
    prae_line = '';
    dist_line = '';
end
function [option] = input_options(varargin)
trail = varargin{1};
%% checking which input option there are
for i = 1:length(trail)
    if strcmp(trail{i}(end-2:end),'num')
        spot.num(i) = i;
    else
        spot.num(i) = 0;
    end
    if strcmp(trail{i}(end-1:end),'.m')
        spot.m(i) = i;
    else
        spot.m(i) = 0;
    end
    if ~isempty(findstr(trail{i},'.tex'))
        spot.tex(i) = i;
    else
        spot.tex(i) = 0;
    end
    if strcmp(trail{i}(1:3),'bwc')
        spot.bwc(i) = i;
    else
        spot.bwc(i) = 0;
    end
end
% makes sure there are no double extensions or options
if length(find(spot.num)) == 1
    spot.num = spot.num(find(spot.num));
else
    spot.num = [];
end
if length(find(spot.m)) == 1
    spot.m = spot.m(find(spot.m));
else
    spot.m = [];
end
if length(find(spot.tex)) == 1
    spot.tex = spot.tex(find(spot.tex));
else
    spot.tex = [];
end
if length(find(spot.bwc)) == 1
    spot.bwc = spot.bwc(find(spot.bwc));
else
    spot.bwc = [];
end
%% about numbered code lines
if ~isempty(spot.num)
    option.numb = trail{spot.num};
    if ~strcmp(option.numb,'num') && ~strcmp(option.numb,'no_num')
        option.numb = 'no_num';
        disp('The code will not be numbered, because the option is wrong.')
    end
else
    % no option for numbered lines stated, ask for that
    button = questdlg('Do you want numbered code lines?','Options','Yes','No','Yes');
    if strcmp(button,'Yes')
        option.numb = 'num';
    else
        option.numb = 'no_num';
    end
end
%% about source path and name
if ~isempty(spot.m)
    idx = findstr(trail{spot.m},'\');
    if isempty(idx)
        idx = findstr(trail{spot.m},'/');
    end
    if ~isempty(idx)
        option.srcname = trail{spot.m}(idx(end)+1:end);
        option.srcpath = trail{spot.m}(1:idx(end));
    else
        option.srcname = trail{spot.m};
        option.srcpath = [cd,'\'];
    end
    option.FilterIndex = 1;
else
    % menu to choose m-file for conversion
    [option.srcname,option.srcpath,option.FilterIndex] = uigetfile('*.m','Select the m-file');
end
%% about target path and name
if ~isempty(spot.tex)
    idx = findstr(trail{spot.tex},'\');
    if isempty(idx)
        idx = findstr(trail{spot.tex},'/');
    end
    if ~isempty(idx)
        option.tgtname = trail{spot.tex}(idx(end)+1:end);
        option.tgtpath = trail{spot.tex}(1:idx(end));
    else
        option.tgtname = trail{spot.tex};
        option.tgtpath = [cd,'\'];
    end
    option.FilterIndex = 1;
else
    option.tgtname = [option.srcname(1:end-2),'.tex'];
    option.tgtpath = option.srcpath;
    option.FilterIndex = 1;
    % menu to choose destination and name for tex-file to save
%     FileName_tex = [trail{spot.m}(idx(end)+1:end-2),'.tex'];
%     [option.tgtname,option.tgtpath,option.FilterIndex] = uiputfile('*.m','Select the m-file',FileName_tex);
end
%% about black white code
if ~isempty(spot.bwc)
    option.bwc = trail{spot.bwc};
else
    option.bwc = 'color';
end

