function varargout = GUI_ellipsoid(varargin)
% GUI_ELLIPSOID MATLAB code for GUI_ellipsoid.fig
%      GUI_ELLIPSOID, by itself, creates a new GUI_ELLIPSOID or raises the existing
%      singleton*.
%
%      H = GUI_ELLIPSOID returns the handle to a new GUI_ELLIPSOID or the handle to
%      the existing singleton*.
%
%      GUI_ELLIPSOID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ELLIPSOID.M with the given input arguments.
%
%      GUI_ELLIPSOID('Property','Value',...) creates a new GUI_ELLIPSOID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ellipsoid_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ellipsoid_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ellipsoid

% Last Modified by GUIDE v2.5 17-Sep-2018 19:12:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ellipsoid_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ellipsoid_OutputFcn, ...
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


% --- Executes just before GUI_ellipsoid is made visible.
function GUI_ellipsoid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ellipsoid (see VARARGIN)

% Choose default command line output for GUI_ellipsoid
handles.output = hObject;

handles.X0 = [0;0];
run(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ellipsoid wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ellipsoid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Q = eval(get(handles.edit4, 'String'));
handles.R = str2double(get(handles.edit5, 'String'));

[x10,x20] = ginput(1);
handles.X0 = [x10; x20];

run(hObject, eventdata, handles)



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
    run(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
    run(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
    run(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Other Functions

function run(hObject, eventdata, handles)
    handles.Q = eval(get(handles.edit4, 'String'));
    handles.R = str2double(get(handles.edit5, 'String'));
    
    A = [0 1; 0 0]; B = [0; 1]; C = eye(2); D = 0;
    sys = ss(A,B,C,D);

    % Discrete Plant
    Ts = 0.01;
    sysd = c2d(sys,Ts,'zoh');
    Ad = sysd.a; Bd = sysd.b; Cd = sysd.c; Dd = sysd.d;

    handles.K = dlqr(Ad, Bd, handles.Q, handles.R);
    Acl = Ad-Bd*handles.K;
    P = dlyap(Acl', handles.Q');

    nonlcon = @(X) nonlconstraints (X, Acl);
    x0 = reshape(handles.Q,4,1);

    options = struct('MaxFunctionEvaluations', 1000, 'MaxIterations', 1000);

    fun = @(X) objFunction (X, Acl);
    x = fmincon(fun,x0,[],[],[],[],[],[],nonlcon,options);
    Qnew = reshape(x,2,2);
    Pnew = dlyap(Acl', Qnew');
    axes(handles.axes1); cla; hold on;
    rectangle('Position', [-1 -1 2 2]); 
    Pnew = inv(Pnew);
    E = ellipsoid([0; 0], Pnew); 
    plot(E)
    axis([-1.5 1.5 -1.5 1.5]);
    axis equal
    
    X0 = handles.X0;
    
    t = 0:Ts:10;
    K = handles.K;
    sys_cl = ss(Ad-Bd*K,B,C,D,Ts);
    [X, t] = lsim(sys_cl,zeros(size(t)), t, X0);

    plot(X(:,1), X(:,2), 'LineWidth',2)

    X0 = round(X0,2);
    axes(handles.axes3); plot(t,X(:,1)); title(['Initial Condition: ', mat2str(X0)])
    axes(handles.axes4); plot(t,X(:,2)); xlabel('Time [s]')

% end

function out = objFunction (X, Acl)
    Q = reshape(X, 2,2);
    P = dlyap(Acl', Q');
    out = trace(P);
% end

function [c, ceq] = nonlconstraints (X, Acl)
    c = [];
    Q = reshape(X, 2,2);

    P = dlyap(Acl', Q');
    gama1 = sqrt([1 0]*inv(P)*[1 0]');
    gama3 = sqrt([0 1]*inv(P)*[0 1]');
        
    ceq = X(2)-X(3);
    c = [c; -Q(1,1); -Q(1,1)*Q(2,2)+Q(2,1)*Q(1,2); gama1-1; gama3-1];   
% end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

run(hObject, eventdata, handles);
