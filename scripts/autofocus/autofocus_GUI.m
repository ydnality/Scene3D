function varargout = autofocus_GUI(varargin)
% AUTOFOCUS_GUI MATLAB code for autofocus_GUI.fig
%      AUTOFOCUS_GUI, by itself, creates a new AUTOFOCUS_GUI or raises the existing
%      singleton*.
%
%      H = AUTOFOCUS_GUI returns the handle to a new AUTOFOCUS_GUI or the handle to
%      the existing singleton*.
%
%      AUTOFOCUS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOFOCUS_GUI.M with the given input arguments.
%
%      AUTOFOCUS_GUI('Property','Value',...) creates a new AUTOFOCUS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before autofocus_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to autofocus_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help autofocus_GUI

% Last Modified by GUIDE v2.5 25-Mar-2013 00:16:43
 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @autofocus_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @autofocus_GUI_OutputFcn, ...
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



% --- Executes just before autofocus_GUI is made visible.
function autofocus_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to autofocus_GUI (see VARARGIN)

% Choose default command line output for autofocus_GUI
handles.output = hObject;
fprintf('Opening Fcn\n');

% --- Initialize the scene and show the image.
%Render cones scene
sceneName = 'cones autofocus';
oi = s3dRenderScene('../../autofocus/cones/autofocusTest.pbrt', .050);
oi = oiSet(oi, 'name', sceneName);
% vcAddAndSelectObject(oi);
% oiWindow;

rgbImage = oiGet(oi, 'rgb image');
 
% Show the image on axes1.
axes(handles.axes1);
imHandle = imshow(rgbImage);

set(imHandle,'ButtonDownFcn',@image_ButtonDownFcn);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes autofocus_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = autofocus_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
[x,y] = ginput(1);
xmin = x-10; xmax = x+10;
ymin = y-10; ymax = y+10;
%fprintf('x:%d , y:%d\n',x,y);
%fprintf('x:%d , y:%d\n',xres,yres);
if x<0 || y<0 || x>xres || y>yres
    fprintf('ERROR: Out of image!\n');
    return;
elseif 0<=x && x<10 || 0<=y && y<10 || xres-10<x && x<=xres || yres-10<y && y<=yres
    fprintf('ERROR: Can not get portion image!\n');
    return;
else
    smallImage = rgbImage(xmin:xmax,ymin:ymax,:);
    %imshow(smallImage);
    tempVar = var(smallImage(:));
end
fprintf('axes buttondown fcn\n');

% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('figure button down\n');

function image_ButtonDownFcn(vargin)

[x,y] = ginput(1);
xmin = x-10; xmax = x+10;
ymin = y-10; ymax = y+10;
%fprintf('x:%d , y:%d\n',x,y);
%fprintf('x:%d , y:%d\n',xres,yres);
if x<0 || y<0 || x>xres || y>yres
    fprintf('ERROR: Out of image!\n');
    return;
elseif 0<=x && x<10 || 0<=y && y<10 || xres-10<x && x<=xres || yres-10<y && y<=yres
    fprintf('ERROR: Can not get portion image!\n');
    return;
else
    smallImage = rgbImage(xmin:xmax,ymin:ymax,:);
    %imshow(smallImage);
    tempVar = var(smallImage(:));
end
fprintf('image Button Down Fcn');
