function varargout = TaskStimSelector(varargin)
% TASKSTIMSELECTOR MATLAB code for TaskStimSelector.fig
%      TASKSTIMSELECTOR, by itself, creates a new TASKSTIMSELECTOR or raises the existing
%      singleton*.
%
%      H = TASKSTIMSELECTOR returns the handle to a new TASKSTIMSELECTOR or the handle to
%      the existing singleton*.
%
%      TASKSTIMSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TASKSTIMSELECTOR.M with the given input arguments.
%
%      TASKSTIMSELECTOR('Property','Value',...) creates a new TASKSTIMSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TaskStimSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TaskStimSelector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TaskStimSelector

% Last Modified by GUIDE v2.5 15-May-2023 18:27:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TaskStimSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @TaskStimSelector_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
    %gui_State.gui_Callback = varargin{1};
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TaskStimSelector is made visible.
function TaskStimSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TaskStimSelector (see VARARGIN)

% Choose default command line output for TaskStimSelector
handles.output = hObject;

if ~isempty(varargin) 
end

%keyboard
handles.basedir = varargin{2};
handles.filelist = varargin{3};

set(handles.task_dir_edit, 'String', handles.basedir)

if ~isempty(handles.filelist)
    if ~iscell(handles.filelist) 
        handles.filelist = {handles.filelist}; 
    end
    set(handles.task_file_list, 'String', handles.filelist); 
    handles.filelist_full = fullfile(handles.basedir, handles.filelist); 
    if ~iscell(handles.filelist_full) 
        handles.filelist_full = {handles.filelist_full}; 
    end 
else
    handles.filelist = []; 
    handles.filelist_full = []; 
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TaskStimSelector wait for user response (see UIRESUME)
% uiwait(handles.TaskStimSelectorGUI);


% --- Outputs from this function are returned to the command line.
function varargout = TaskStimSelector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in task_file_list.
function task_file_list_Callback(hObject, eventdata, handles)
% hObject    handle to task_file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns task_file_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from task_file_list


% --- Executes during object creation, after setting all properties.
function task_file_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to task_file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add_push.
function add_push_Callback(hObject, eventdata, handles)
% hObject    handle to add_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('Add task file', '*.mat', handles.basedir);
if ~isequal(filename,0)
    tflist = get(handles.task_file_list, 'String');
    tflist = [tflist; {filename}];
    set(handles.task_file_list, 'String', tflist)
    if ~isempty(handles.filelist_full)
        handles.filelist_full = [handles.filelist_full; {fullfile(pathname, filename)}];
        handles.filelist = [handles.filelist; {filename}];
    else
        handles.filelist_full = {fullfile(pathname, filename)};
        handles.filelist = {filename};
    end
end
guidata(hObject, handles);


% --- Executes on button press in remove_push.
function remove_push_Callback(hObject, eventdata, handles)
% hObject    handle to remove_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rmidx = get(handles.task_file_list, 'Value'); 
handles.filelist(rmidx) = [];
handles.filelist_full(rmidx) = [];
set(handles.task_file_list, 'String', handles.filelist);
guidata(hObject, handles); 


function task_dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to task_dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of task_dir_edit as text
%        str2double(get(hObject,'String')) returns contents of task_dir_edit as a double


% --- Executes during object creation, after setting all properties.
function task_dir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to task_dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
