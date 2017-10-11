function varargout = Exploratory_Data_Analysis_GUI_v02(varargin)
% exploratory_data_analysis_gui_v02 MATLAB code for Exploratory_Data_Analysis_GUI_v02.fig
%      exploratory_data_analysis_gui_v02, by itself, creates a new exploratory_data_analysis_gui_v02 or raises the existing
%      singleton*.
%
%      H = exploratory_data_analysis_gui_v02 returns the handle to a new exploratory_data_analysis_gui_v02 or the handle to
%      the existing singleton*.
%
%      exploratory_data_analysis_gui_v02('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in exploratory_data_analysis_gui_v02.M with the given input arguments.
%
%      exploratory_data_analysis_gui_v02('Property','Value',...) creates a new exploratory_data_analysis_gui_v02 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Exploratory_Data_Analysis_GUI_v02_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Exploratory_Data_Analysis_GUI_v02_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Exploratory_Data_Analysis_GUI_v02

% Last Modified by GUIDE v2.5 01-Sep-2017 18:14:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Exploratory_Data_Analysis_GUI_v02_OpeningFcn, ...
                   'gui_OutputFcn',  @Exploratory_Data_Analysis_GUI_v02_OutputFcn, ...
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


% --- Executes just before Exploratory_Data_Analysis_GUI_v02 is made visible.
function Exploratory_Data_Analysis_GUI_v02_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Exploratory_Data_Analysis_GUI_v02 (see VARARGIN)



%% Test required dependencies 
if ~exist('plot_google_map','file')
    errordlg('plot_google_map.m function from Matlab File Exchange is required to be in the working paths','plot_google_map error')
end

if ~exist('dynamicDateTicks','file')
    errordlg('dynamicDateTicks.m function from Matlab File Exchange is required to be in the working paths','dynamicDateTicks error')
end

if ~exist('diverging_map','file')
    errordlg('diverging_map.m function must be in the working paths. It is available for download @ http://www.kennethmoreland.com/color-maps/','diverging color maps error')
end

if ~exist('export_fig','file')
    errordlg('export_fig.m function from Matlab File Exchange is required to be in the working paths','export_fig error')
end

if ~exist('create_Drive_Pass_Index','file')
    errordlg('create_Drive_Pass_Index.m function must be in the working path. This is availble at https://github.com/kpmessier/SV-Library/tree/master/SV_Library2.0/Basic_Functions','drive pass error')
end

dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});



% Choose default command line output for Exploratory_Data_Analysis_GUI_v02
handles.output = hObject;

 % create the listener for the slider
 handles.sliderListener1 = addlistener(handles.slider4,'ContinuousValueChange', ...
                                      @(hFigure,eventdata) slider_minhour_ContValCallback(...
                                        hObject,eventdata));
                                    
 handles.sliderListener2 = addlistener(handles.slider5,'ContinuousValueChange', ...
                                      @(hFigure,eventdata) slider_maxhour_ContValCallback(...
                                        hObject,eventdata)); 
                                    
 handles.sliderListener3 = addlistener(handles.slider1,'ContinuousValueChange', ...
                                      @(hFigure,eventdata) slider_day_ContValCallback(...
                                        hObject,eventdata));                                      
 % update handles structure
 guidata(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Exploratory_Data_Analysis_GUI_v02 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function [txt] = myupdatefcn(~,event_obj,handles)


UserData = handles.figure1.UserData;
tdata = UserData.time;

ytemp1 = UserData.y;
if ~isfield(UserData,'Index2')
    Index2=true(length(ytemp1),1);
else
    Index2 = UserData.Index2;
end
if ~isfield(UserData,'date1')
    Index_Date1 = true(length(ytemp1),1);
else
    Index_Date1 = tdata >= datenum(UserData.date1);
end
if ~isfield(UserData,'date2')
    Index_Date2 = true(length(ytemp1),1);
else
    Index_Date2 = tdata <= datenum(UserData.date2);
end

idx_reduce = UserData.index_reduce;

% uidx_reduce = unique(idx_reduce);
% n_days = accumarray(idx_reduce,floor(tdata),[],@(x)length(unique(x)),NaN);
% min_days =  handles.Min_Days_Text.Value;
% idx_days = n_days >= min_days;
% idx_remove = ~ismember(idx_reduce,uidx_reduce(idx_days));


min_hour_index = hour(tdata)>= handles.Min_Slider_Text.Value;
max_hour_index = hour(tdata)<= handles.Max_Slider_Text.Value;
idx_hour = min_hour_index & max_hour_index;

Index_Final =  Index2 & Index_Date1 & Index_Date2 & idx_hour;
ytemp2 = ytemp1(Index_Final);

pos = event_obj.Position;
Xcoord = event_obj.Target.XData';
Ycoord = event_obj.Target.YData';
% Xcoord_all = UserData.xcoord;
% Ycoord_all = UserData.ycoord;
ull = UserData.ull;
[~,idx1] = ismember(pos,[Xcoord Ycoord],'rows');
[~,idx2] = ismember(pos,ull,'rows');
Ydata = event_obj.Target.CData(idx1);
idx_data = idx_reduce == idx2;
if isfield(UserData,'DrivePassIndex')
    ytemp = ytemp2(idx_data); 
    driveIndex = UserData.DrivePassIndex; driveIndex = driveIndex(idx_data);
    driveFunc = UserData.driveFunc;   
    yplot = grpstats(ytemp,driveIndex,driveFunc);
else
    yplot = ytemp2(idx_data);
end
tdates =  datestr(unique(floor(UserData.time(idx_data))));
set(handles.DataListBox,'Value',1)
set(handles.DataListBox,'String',tdates);
if ~isfield(UserData,'mapFunc')
    mapFunc = @(x)median(x,'omitnan');
else
    mapFunc=UserData.mapFunc;
end
yFunc = mapFunc(yplot);
yinterp = linspace(UserData.yrange(1),UserData.yrange(2),size(UserData.cmap,1));
idxy = knnsearch(yinterp',yFunc);

hist_color = UserData.cmap(idxy,:);
h=histogram(handles.CursorFig,yplot,'BinMethod','fd');
set(h,'FaceColor',hist_color);
set(h,'FaceAlpha',0.9);
grid(handles.CursorFig,'on')
ylabel(handles.CursorFig,'Count')
if isfield(UserData,'units')
    xlabel(handles.CursorFig,UserData.units)
end



Time_1hz = UserData.time(idx_data);
TimeDay = floor(Time_1hz);
ud = unique(TimeDay);
k=1;
yp = UserData.y(idx_data);
for i = 1:length(ud)
    idx = TimeDay==ud(i);
    tplot(k:k+sum(idx),1)=[yp(idx); NaN];
    filled_time(k:k+sum(idx),1)=[Time_1hz(idx); NaN];
    k=length(tplot)+1;
end

if isfield(UserData,'TimeSeries_t') && isfield(UserData,'TimeSeries_y')
    set(handles.TimeSeriesFig,'NextPlot','replace')
  plot(handles.TimeSeriesFig,UserData.TimeSeries_t,UserData.TimeSeries_y,'-b');
end
   

set(handles.TimeSeriesFig,'NextPlot','add')
plot(handles.TimeSeriesFig,filled_time,tplot,'-r','linew',3);
if isfield(UserData,'Time_YLim')
    ylim(handles.TimeSeriesFig,UserData.Time_YLim);
end
dynamicDateTicks(handles.TimeSeriesFig);

% zoom_obj = zoom(handles.ColorMapFig);
% set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});

if isfield(UserData,'DrivePassIndex')
    txt = {['Long: ',num2str(pos(1))], ['Lat:  ',num2str(pos(2))],...
        ['Val:  ',num2str(Ydata)],['# of 1-Hz obs: ',num2str(sum(idx_data))],...
        ['# of Drive Passes: ',num2str(length(yplot))],['# of Days: ',num2str(size(tdates,1))]};    
    
else
    txt = {['Long: ',num2str(pos(1))], ['Lat:  ',num2str(pos(2))],...
        ['Val:  ',num2str(Ydata)],['# of obs: ',num2str(sum(idx_data))],...
        ['# of Days: ',num2str(size(tdates,1))]};
end


function zoomUpdateFcn(~,event_obj,handles)
set(handles.figure1, 'pointer', 'watch')

UserData = get(handles.output,'UserData');

if ~isfield(UserData,'y')
    errordlg('Must input YData','YData Input Error')
else
    y = UserData.y;
end
if ~isfield(UserData,'time')
    errordlg('Time is a flat circle','Time Error Input')
else
    tdata = UserData.time;
end
if ~isfield(UserData,'units')
    units='';
else
    units=UserData.units;
end
if ~isfield(UserData,'xcoord')
    error('Must input X Coordinate')
else
    xcoord = UserData.xcoord;
end
if ~isfield(UserData,'ycoord')
    error('Must input Y Coordinate')
else
    ycoord = UserData.ycoord;
end
% if ~isfield(UserData,'timeFunc')
%     timeFunc =  @(x)median(x,'omitnan');
% else
%     timeFunc=UserData.timeFunc;
% end
if ~isfield(UserData,'TimeTitle')
    titlestr='';
else
    titlestr = UserData.TimeTitle;
end
if ~isfield(UserData,'Index1')
    Index1=true(length(y),1);
else
    Index1 = UserData.Index1;
end
if ~isfield(UserData,'Index2')
    Index2=true(length(y),1);
else
    Index2 = UserData.Index2;
end
if ~isfield(UserData,'date1')
    Index_Date1 = true(length(y),1);
else
    Index_Date1 = tdata >= datenum(UserData.date1);
end
if ~isfield(UserData,'date2')
    Index_Date2 = true(length(y),1);
else
    Index_Date2 = tdata <= datenum(UserData.date2);
end

min_hour_index = hour(tdata)>= handles.Min_Slider_Text.Value;
max_hour_index = hour(tdata)<= handles.Max_Slider_Text.Value;
idx_hour = min_hour_index & max_hour_index;

% Window axis index
XLim = handles.ColorMapFig.XLim;
YLim = handles.ColorMapFig.YLim;
ll = [xcoord ycoord];
[ull,~,index_ax] = unique(ll,'rows');
[in,on] = inpolygon(ull(:,1),ull(:,2),XLim,YLim);
idx_axis = in | on;
Index_Axis = idx_axis(index_ax);
%
Index_Final =  Index2 & Index_Date1 & Index_Date2 & Index_Axis & idx_hour;

tdata = tdata(Index_Final);
y = y(Index_Final);

TimeDay = floor(tdata);
ud = unique(TimeDay);
k=1;
for i = 1:length(ud)
    idx = TimeDay==ud(i);
    yplot(k:k+sum(idx),1)=[y(idx); NaN];
    filled_time(k:k+sum(idx),1)=[tdata(idx); NaN];
    k=length(yplot)+1;
end

axes(handles.TimeSeriesFig);
hold off
plot(filled_time,yplot,'-b');
ht=title(titlestr);
set(ht,'BackgroundColor',[0.94 0.94 0.94])
ylabel(units)
if isfield(UserData,'Time_YLim')
    ylim(UserData.Time_YLim)
end
dynamicDateTicks;
UserData.TimeSeries_t = filled_time;
UserData.TimeSeries_y = yplot;

axes(handles.ColorMapFig);
plot_google_map('MapType','terrain','ShowLabels',false)

zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});
set(handles.output,'UserData',UserData);
set(handles.figure1, 'pointer', 'arrow');


% --- Outputs from this function are returned to the command line.
function varargout = Exploratory_Data_Analysis_GUI_v02_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Enter_YData_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

y_new =  get(hObject,'String');
set(hObject,'ForegroundColor','k')
UserData = get(handles.output,'UserData');
try
    UserData.y = evalin('base',y_new);
catch
    errordlg(sprintf('%s does not exist as a variable in the workspace!',y_new),'Y Input error')
end
set(handles.output,'UserData',UserData);
% --- Executes during object creation, after setting all properties.
function Enter_YData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Enter_XCoord_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

xcoord_new =  get(hObject,'String');
set(hObject,'ForegroundColor','k')
UserData = get(handles.output,'UserData');
try
    UserData.xcoord = evalin('base',xcoord_new);
catch
    errordlg(sprintf('%s does not exist as a variable in the workspace!',xcoord_new),'X Coordinate error')
end
set(handles.output,'UserData',UserData);

% --- Executes during object creation, after setting all properties.
function Enter_XCoord_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Enter_YCoord_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
ycoord_new =  get(hObject,'String');
UserData = get(handles.output,'UserData');
set(hObject,'ForegroundColor','k')
try
    UserData.ycoord = evalin('base',ycoord_new);
    set(hObject,'String',ycoord_new);
catch
    errordlg(sprintf('%s does not exist as a variable in the workspace!',ycoord_new),'Y coordinate error')
end
set(handles.output,'UserData',UserData);

% --- Executes during object creation, after setting all properties.
function Enter_YCoord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function Enter_Units_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
units =  get(hObject,'String');
UserData = get(handles.output,'UserData');

UserData.units = units;
set(handles.output,'UserData',UserData);


% --- Executes during object creation, after setting all properties.
function Enter_Units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Time_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Time_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Time_Edit as text
%        str2double(get(hObject,'String')) returns contents of Time_Edit as a double
 date_str = get(hObject,'String');
 set(hObject,'ForegroundColor','k')
 UserData = get(handles.output,'UserData');
try
    UserData.time = evalin('base',date_str);
    set(hObject,'String',date_str);
    
catch
    errordlg(sprintf('Time is a flat circle \n %s does not exist as a variable in the workspace!',date_str),'Time Input error')
end
set(handles.output,'UserData',UserData);

% --- Executes during object creation, after setting all properties.
function Time_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Time_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function Mapping_Functions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
contents = cellstr(get(hObject,'String'));
func = contents{get(hObject,'Value')};
if strcmp(func,'Mean')
    mapFunc = @(x)mean(x,'omitnan');
elseif strcmp(func,'Median')
    mapFunc = @(x)median(x,'omitnan');
elseif strcmp(func,'5th Percentile')
    mapFunc = @(x)prctile(x,5);
elseif strcmp(func,'25th Percentile')
    mapFunc = @(x)prctile(x,25);
elseif strcmp(func,'75th Percentile')
    mapFunc = @(x)prctile(x,75);
elseif strcmp(func,'95th Percentile')
    mapFunc = @(x)prctile(x,95);
elseif strcmp(func,'Min')
    mapFunc = @(x)min(x,'omitnan');
elseif strcmp(func,'Max')
    mapFunc = @(x)max(x,'omitnan');
elseif strcmp(func,'Number of Unique')
    mapFunc = @(x)length(unique(x));
elseif strcmp(func,'Variance')
    mapFunc = @(x)var(x,'omitnan');
elseif strcmp(func,'Standard Deviation')
    mapFunc = @(x)std(x,'omitnan');
elseif strcmp(func,'IQR')
    mapFunc = @(x)iqr(x);
end
 UserData = get(handles.output,'UserData');
 UserData.mapFunc = mapFunc;
set(handles.output,'UserData',UserData);

% --- Executes during object creation, after setting all properties.
function Mapping_Functions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function CreateColorMap_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1, 'pointer', 'watch')
drawnow;

UserData = get(handles.output,'UserData');

if ~isfield(UserData,'y')
    error('Must input YData')
else
    y = UserData.y;
end
if ~isfield(UserData,'xcoord')
    error('Must input X Coordinate')
else
    xcoord = UserData.xcoord;
end
if ~isfield(UserData,'ycoord')
    error('Must input Y Coordinate')
else
    ycoord = UserData.ycoord;
end
if ~isfield(UserData,'time')
    error('Time is a flat circle')
else
    tdata = UserData.time;
end
if ~isfield(UserData,'veh_id')
    error('Vehicle Index is a required field!')
end

if ~isfield(UserData,'units')
    units='';
else
    units=UserData.units;
end
if ~isfield(UserData,'mapFunc')
    mapFunc = @(x)median(x,'omitnan');
else
    mapFunc=UserData.mapFunc;
end
if ~isfield(UserData,'MarkerSize')
    MarkerSize = 8;
else
    MarkerSize=UserData.MarkerSize;
end
if ~isfield(UserData,'shape')
    shape = 'o';
else
    shape=UserData.shape;
end
if ~isfield(UserData,'yrange')
    yrange = [0 prctile(y,95)];
    UserData.yrange=yrange;
else
    yrange=UserData.yrange;
end
if ~isfield(UserData,'cmap')
    cmap = flipud(hot(256));
    cmap(1:20,:)=[];
    UserData.cmap = cmap;
else
    cmap = UserData.cmap;
end
if ~isfield(UserData,'MapTitle')
    titlestr='';
else
    titlestr = UserData.MapTitle;
end

if ~isfield(UserData,'Index2')
    Index2=true(length(y),1);
else
    Index2 = UserData.Index2;
end
if ~isfield(UserData,'date1')
    Index_Date1 = true(length(y),1);
else
    Index_Date1 = tdata >= datenum(UserData.date1);
end
if ~isfield(UserData,'date2')
    Index_Date2 = true(length(y),1);
else
    Index_Date2 = tdata <= datenum(UserData.date2);
end

min_hour_index = hour(tdata)>= handles.Min_Slider_Text.Value;
max_hour_index = hour(tdata)<= handles.Max_Slider_Text.Value;
idx_hour = min_hour_index & max_hour_index;



Index_Final =  Index2 & Index_Date1 & Index_Date2 & idx_hour;
y = y(Index_Final);

ll = [xcoord ycoord];
ll = ll(Index_Final,:);

[ull,~,index] = unique(ll,'rows');


if isfield(UserData,'DrivePassIndex')
    index_pass = UserData.DrivePassIndex; index_pass = index_pass(Index_Final);
    Pass_Mean = accumarray(index_pass,y,[],@nanmean,NaN);
    Pass_idx = accumarray(index_pass,index,[],@(x)x(1),NaN);
    try
    y_reduce = accumarray(Pass_idx,Pass_Mean,[],mapFunc,NaN);
    catch 
        idx_nan = isnan(Pass_idx); 
        y_reduce = accumarray(Pass_idx(~idx_nan),Pass_Mean(~idx_nan),[],mapFunc,NaN);
    end
else
    try 
        y_reduce = accumarray(index,y,[],mapFunc,NaN);
    catch
        idx_nan = isnan(index);
        y_reduce = accumarray(index(~idx_nan),y(~idx_nan),[],mapFunc,NaN);
    end
end

% Get min-days for plotting index 
min_days =  handles.Min_Days_Text.Value;
n_days = accumarray(index,floor(tdata(Index_Final)),[],@(x)length(unique(x)),NaN);
idx_days = n_days >= min_days;

axes(handles.ColorMapFig);
XLim = get(handles.ColorMapFig,'XLim');
YLim = get(handles.ColorMapFig,'YLim');
hold off
scatter(ull(idx_days,1),ull(idx_days,2),MarkerSize,y_reduce(idx_days),'filled',shape);
% Keep current zoom axes
if XLim(1,1) == 0 && XLim(1,2) == 1 && YLim(1,1) == 0 && YLim(1,2) == 1
else
    axis([XLim YLim])
end
title(titlestr)
colormap(cmap);
hc=colorbar;
title(hc,units)
set(hc,'FontName','Arial')
caxis(yrange)
plot_google_map('MapType','terrain','ShowLabels',false)

zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});

% plot_google_map('MapType','terrain','ShowLabels',false,'APIKey','AIzaSyCNCMyVWdDqocbEjUyWIMtQLz6xfHgRd-k')
UserData.ull = ull;
UserData.index_reduce = index;
UserData.index_st = Index_Final;
set(handles.output,'UserData',UserData);
set(handles.figure1, 'pointer', 'arrow');
    % --- Executes during object creation, after setting all properties.
function CreateColorMap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes during object creation, after setting all properties.
function ColorMapFig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColorMapFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ColorMapFig







function ColorScaleMin_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

cax_min = str2double(get(hObject,'String'));
UserData = get(handles.output,'UserData');
UserData.yrange(1) = cax_min;
set(handles.output,'UserData',UserData);

% --- Executes during object creation, after setting all properties.
function ColorScaleMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ColorScaleMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
cax_max = str2double(get(hObject,'String'));
UserData = get(handles.output,'UserData');
UserData.yrange(2) = cax_max;
set(handles.output,'UserData',UserData);

% --- Executes during object creation, after setting all properties.
function ColorScaleMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------


% --- Executes on selection change in cmap.
function cmap_Callback(hObject, eventdata, handles)
% hObject    handle to cmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap

contents =  cellstr(get(hObject,'String'));
c1=contents{get(hObject,'Value')};

UserData = get(handles.output,'UserData');

if strcmp(c1,'Hot')
    cmap = flipud(hot(256));
    cmap(1:20,:)=[];
elseif strcmp(c1,'Blue-to-Red')
    cmap = diverging_map(0:1/256:1,[0.230, 0.299, 0.754],[0.706, 0.016, 0.150]);
elseif strcmp(c1,'Grayscale')
    cmap = flipud(gray(256));
elseif strcmp(c1,'Autumn')
    cmap = flipud(autumn(256));
elseif strcmp(c1,'hsv')
    cmap = hsv(256);
end
UserData.cmap = cmap;
colormap(cmap);
set(handles.output,'UserData',UserData);

zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});

% --- Executes during object creation, after setting all properties.
function cmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Index1_Callback(hObject, eventdata, handles)
% hObject    handle to Index1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Index1 as text
%        str2double(get(hObject,'String')) returns contents of Index1 as a double
index1 =  get(hObject,'String');
set(hObject,'ForegroundColor','k')
UserData = get(handles.output,'UserData');
try
    UserData.veh_id = evalin('base',index1);
catch
    error('%s does not exist as a variable in the workspace!',index1)
end
set(handles.output,'UserData',UserData);
zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});
% --- Executes during object creation, after setting all properties.
function Index1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Index1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Index2_Callback(hObject, eventdata, handles)
% hObject    handle to Index2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Index2 as text
%        str2double(get(hObject,'String')) returns contents of Index2 as a double
index2 =  get(hObject,'String');
set(hObject,'ForegroundColor','k')
UserData = get(handles.output,'UserData');
try
    UserData.Index2 = evalin('base',index2);
catch
    error('%s does not exist as a variable in the workspace!',index2)
end
set(handles.output,'UserData',UserData);
zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});

% --- Executes during object creation, after setting all properties.
function Index2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Index2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TimeYLim_Callback(hObject, eventdata, handles)
% hObject    handle to TimeYLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeYLim as text
%        str2double(get(hObject,'String')) returns contents of TimeYLim as a double
Time_YLim =  get(hObject,'String');
YLim_Dbl = str2num(Time_YLim);%#ok
if size(YLim_Dbl,1) ~= 1 && size(YLim_Dbl,2) ~= 2
    error('Invalid Time Series Y Limits. Time Series YLim must be a 1 x 2 numeric matrix')
end
UserData = get(handles.output,'UserData');
UserData.Time_YLim = YLim_Dbl;
set(hObject,'ForegroundColor','k')
set(handles.output,'UserData',UserData);zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});

% --- Executes during object creation, after setting all properties.
function TimeYLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeYLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DateRange1_Callback(hObject, eventdata, handles)
% hObject    handle to DateRange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DateRange1 as text
%        str2double(get(hObject,'String')) returns contents of DateRange1 as a double
 date1 = datestr(get(hObject,'String'));
 set(hObject,'String',date1);
UserData = get(handles.output,'UserData');
UserData.date1 = date1;
set(hObject,'ForegroundColor','k')
set(handles.output,'UserData',UserData);
zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});

% --- Executes during object creation, after setting all properties.
function DateRange1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DateRange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DateRange2_Callback(hObject, eventdata, handles)
% hObject    handle to DateRange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DateRange2 as text
%        str2double(get(hObject,'String')) returns contents of DateRange2 as a double

 date2 = datestr(get(hObject,'String'));
 set(hObject,'String',date2);
UserData = get(handles.output,'UserData');
UserData.date2 = date2;
set(hObject,'ForegroundColor','k')
set(handles.output,'UserData',UserData);
zoom_obj = zoom(handles.ColorMapFig);
set(zoom_obj,'ActionPostCallback',{@zoomUpdateFcn,handles})
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles});
% --- Executes during object creation, after setting all properties.
function DateRange2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DateRange2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Map_Title_Callback(hObject, eventdata, handles)
% hObject    handle to Map_Title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Map_Title as text
%        str2double(get(hObject,'String')) returns contents of Map_Title as a double
MapTitle =  get(hObject,'String');
UserData = get(handles.output,'UserData');

UserData.MapTitle = MapTitle;
set(handles.output,'UserData',UserData);


% --- Executes during object creation, after setting all properties.
function Map_Title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Map_Title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TimeSeriesTitle_Callback(hObject, eventdata, handles)
% hObject    handle to TimeSeriesTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeSeriesTitle as text
%        str2double(get(hObject,'String')) returns contents of TimeSeriesTitle as a double

TimeTitle =  get(hObject,'String');
UserData = get(handles.output,'UserData');

UserData.TimeTitle = TimeTitle;
set(handles.output,'UserData',UserData);

% --- Executes during object creation, after setting all properties.
function TimeSeriesTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeSeriesTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportColorMap.
function ExportColorMap_Callback(hObject, eventdata, handles)
% hObject    handle to ExportColorMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_title = inputdlg('Please input the title of your figure');
fig_format = inputdlg(sprintf('Please specify the format.\n Options include fig, png, jpg.\n Do NOT include the leading dot(.) \n .Use eps and .pdf only if your systems supports it.'));
export_fig(fig_format{1},handles.ColorMapFig,fig_title{1})


% ------------------------------------------------------------------


% --------------------------------------------------------------------
% function DataCursor_ClickedCallback(hObject, eventdata, handles)
% % hObject    handle to DataCursor (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%  pos = get(eventdata,'Position');
% 

% % --------------------------------------------------------------------
% function DataCursor_OnCallback(hObject, eventdata, handles)
% % hObject    handle to DataCursor (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% h = datacursormode;
% datacursormode on
% h.DisplayStyle = 'datatip';
% h.SnapToData = 'on';


% -------------------------------------------------------------------


% --- Executes on selection change in DataListBox.
function hObject=DataListBox_Callback(hObject, eventdata, handles)
% hObject    handle to DataListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DataListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DataListBox
% --- Executes during object creation, after setting all properties.
function DataListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu6.
function reduceByDrivePass_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6

contents =  cellstr(get(hObject,'String'));
c1=contents{get(hObject,'Value')};

UserData = get(handles.output,'UserData');

if strcmp(c1,'None')
    if isfield(UserData,'DrivePassIndex')
        UserData = rmfield(UserData,'DrivePassIndex');
    end
elseif strcmp(c1,'Mean')
    driveFunc = @nanmean;
    UserData.driveFunc = driveFunc;
    xcoord = UserData.xcoord;
    ycoord = UserData.ycoord;
    tdata = UserData.time;
    veh_id = UserData.veh_id;
    index = create_Drive_Pass_Index(xcoord,ycoord,tdata,veh_id);
    UserData.DrivePassIndex = index;
elseif strcmp(c1,'Median')
    driveFunc = @nanmedian;
    UserData.driveFunc = driveFunc;
    xcoord = UserData.xcoord;
    ycoord = UserData.ycoord;
    tdata = UserData.time;
    veh_id = UserData.veh_id;
    index = create_Drive_Pass_Index(xcoord,ycoord,tdata,veh_id);
    UserData.DrivePassIndex = index;
end


set(handles.output,'UserData',UserData);


% --- Executes during object creation, after setting all properties.
function reduceByDrivePass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

 function slider_minhour_ContValCallback(hFigure,eventdata)
 % test it out - get the handles object and write the current value
 % to the edit box
 handles = guidata(hFigure);
 sliderValue = get(handles.slider4,'Value');
 set(handles.Min_Slider_Text,'String',[num2str(sliderValue),':','00']);
set(handles.Min_Slider_Text,'Value',sliderValue)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


 function slider_maxhour_ContValCallback(hFigure,eventdata)
 % test it out - get the handles object and write the current value
 % to the edit box
 handles = guidata(hFigure);
 sliderValue = get(handles.slider5,'Value');
 set(handles.Max_Slider_Text,'String',[num2str(sliderValue),':00']);
set(handles.Max_Slider_Text,'Value',sliderValue)


 function slider_day_ContValCallback(hFigure,eventdata)
 % test it out - get the handles object and write the current value
 % to the edit box
 handles = guidata(hFigure);
 sliderValue = get(handles.slider1,'Value');
 set(handles.Min_Days_Text,'String',num2str(sliderValue));
set(handles.Min_Days_Text,'Value',sliderValue)




% --- Executes during object creation, after setting all properties.
function Max_Slider_Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_Slider_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
