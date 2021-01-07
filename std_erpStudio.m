% History
% 02/19/2019 Makoto. Mean (SD) display supported.
% 12/29/2016 Makoto. Created.

function varargout = std_erpStudio(varargin)
% STD_ERPSTUDIO MATLAB code for std_erpStudio.fig
%      STD_ERPSTUDIO, by itself, creates a new STD_ERPSTUDIO or raises the existing
%      singleton*.
%
%      H = STD_ERPSTUDIO returns the handle to a new STD_ERPSTUDIO or the handle to
%      the existing singleton*.
%
%      STD_ERPSTUDIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STD_ERPSTUDIO.M with the given input arguments.
%
%      STD_ERPSTUDIO('Property','Value',...) creates a new STD_ERPSTUDIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before std_erpStudio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to std_erpStudio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help std_erpStudio

% Last Modified by GUIDE v2.5 29-Dec-2016 13:45:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @std_erpStudio_OpeningFcn, ...
                   'gui_OutputFcn',  @std_erpStudio_OutputFcn, ...
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


% --- Executes just before std_erpStudio is made visible.
function std_erpStudio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to std_erpStudio (see VARARGIN)

% Import STUDY and ALLEEG.
STUDY  = evalin('base','STUDY');
ALLEEG = evalin('base','ALLEEG');
STUDY = std_readerp(STUDY, ALLEEG, 'design', STUDY.currentdesign, 'clusters', 2:length(STUDY.cluster));
assignin('base', 'STUDY', STUDY);

% Obtain condition names.
variableLabel11 = STUDY.design(STUDY.currentdesign).variable(1,1).label;
variableValue11 = STUDY.design(STUDY.currentdesign).variable(1,1).value;
variableLabel12 = STUDY.design(STUDY.currentdesign).variable(1,2).label;
variableValue12 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
design11 = STUDY.design(STUDY.currentdesign).variable(1,1).value;
design12 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
design21 = STUDY.design(STUDY.currentdesign).variable(1,1).value;
design22 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
design31 = STUDY.design(STUDY.currentdesign).variable(1,1).value;
design32 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
design41 = STUDY.design(STUDY.currentdesign).variable(1,1).value;
design42 = STUDY.design(STUDY.currentdesign).variable(1,2).value;

% Convert names into non-cell format: variable(1,1)
studyDesignVariable11 = STUDY.design(STUDY.currentdesign).variable(1,1).value;
for n = 1:length(studyDesignVariable11)
    if iscell(studyDesignVariable11{n})
        tmpCell = studyDesignVariable11{n};
        tmpString = '';
        for m = 1:length(tmpCell)
            tmpString = [tmpString '_' tmpCell{m}];
        end
        tmpString = tmpString(2:end);
        studyDesignVariable11{n} = tmpString;
    end
end
% Convert names into non-cell format: variable(1,2)
studyDesignVariable12 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
for n = 1:length(studyDesignVariable12)
    if iscell(studyDesignVariable12{n})
        tmpCell = studyDesignVariable12{n};
        tmpString = '';
        for m = 1:length(tmpCell)
            tmpString = [tmpString '_' tmpCell{m}];
        end
        tmpString = tmpString(2:end);
        studyDesignVariable12{n} = tmpString;
    end
end
design11 = studyDesignVariable11;
design12 = studyDesignVariable12;
design21 = studyDesignVariable11;
design22 = studyDesignVariable12;
design31 = studyDesignVariable11;
design32 = studyDesignVariable12;
design41 = studyDesignVariable11;
design42 = studyDesignVariable12;

% Set up popupmenus 11
if any(cell2mat(design11))
    set(handles.design11, 'string', [{'(none)'} design11]');
else
    set(handles.design11, 'string', 'Not present', 'enable', 'off');
end

% Set up popupmenus 12
if any(cell2mat(design12))
    set(handles.design12, 'string', [{'(none)'} design12]');
else
    set(handles.design12, 'string', 'Not present', 'enable', 'off');
end

% Set up popupmenus 21
if any(cell2mat(design21))
    set(handles.design21, 'string', [{'(none)'} design21]');
else
    set(handles.design21, 'string', 'Not present', 'enable', 'off');
end

% Set up popupmenus 22
if any(cell2mat(design22))
    set(handles.design22, 'string', [{'(none)'} design22]');
else
    set(handles.design22, 'string', 'Not present', 'enable', 'off');
end

% Set up popupmenus 31
if any(cell2mat(design31))
    set(handles.design31, 'string', [{'(none)'} design31]');
else
    set(handles.design31, 'string', 'Not present', 'enable', 'off');
end

% Set up popupmenus 32
if any(cell2mat(design32))
    set(handles.design32, 'string', [{'(none)'} design32]');
else
    set(handles.design32, 'string', 'Not present', 'enable', 'off');
end

% Set up popupmenus 41
if any(cell2mat(design41))
    set(handles.design41, 'string', [{'(none)'} design41]');
else
    set(handles.design41, 'string', 'Not present', 'enable', 'off');
end

% Set up popupmenus 42
if any(cell2mat(design42))
    set(handles.design42, 'string', [{'(none)'} design42]');
else
    set(handles.design42, 'string', 'Not present', 'enable', 'off');
end

% Enter default baseline period
set(handles.baselineEdit, 'String', num2str([min(STUDY.cluster(1,3).erptimes) 0]))

% Change the figure title
set(gcf, 'Name', 'std_erpStudio()')

% Choose default command line output for std_erpStudio
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes std_erpStudio wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = std_erpStudio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in clusterList.
function clusterList_Callback(hObject, eventdata, handles)
% hObject    handle to clusterList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clusterList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clusterList

%%%%%%%%%%%%%%%%%%%
%%% Input check %%%
%%%%%%%%%%%%%%%%%%%
design11Idx = get(handles.design11, 'value')-1;
design12Idx = get(handles.design12, 'value')-1;
design21Idx = get(handles.design21, 'value')-1;
design22Idx = get(handles.design22, 'value')-1;
design31Idx = get(handles.design31, 'value')-1;
design32Idx = get(handles.design32, 'value')-1;
design41Idx = get(handles.design41, 'value')-1;
design42Idx = get(handles.design42, 'value')-1;
% If nothing is chosen
if ~any([design11Idx design12Idx design21Idx design22Idx design31Idx design32Idx design41Idx design42Idx])
    error('Define either A and C, or all A, B, C, and D.')
end
% If only A and B (or C and D) are used
if ~any([design11Idx design12Idx design21Idx design22Idx]) | ~any([design31Idx design32Idx design41Idx design42Idx])
    error('Use A and C instead.')
end



%%%%%%%%%%%%%%%%%%%%%%
%%% Show scalp map %%%
%%%%%%%%%%%%%%%%%%%%%%
STUDY  = evalin('base','STUDY');
ALLEEG = evalin('base','ALLEEG');
currentCluster = get(handles.clusterList,'value')+1;
axes(handles.topoCentroidPlot);
std_topoplot(STUDY,ALLEEG,'clusters',currentCluster,'figure','off');



%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show dipole plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%
cla(handles.axialAxes,    'reset')
cla(handles.sagittalAxes, 'reset')
cla(handles.coronalAxes,  'reset')
axes(handles.axialAxes)
std_dipplot(STUDY,ALLEEG, 'clusters', currentCluster, 'mode', 'apart', 'figure', 'off');
view(0,90)
axes(handles.sagittalAxes)
std_dipplot(STUDY,ALLEEG, 'clusters', currentCluster, 'mode', 'apart', 'figure', 'off');
view(90,0)
axes(handles.coronalAxes)
std_dipplot(STUDY,ALLEEG, 'clusters', currentCluster, 'mode', 'apart', 'figure', 'off');
view(0,0)
set(gcf, 'color', [0.66 0.76 1]);
latency = STUDY.cluster(currentCluster).erptimes;
data = cell(4,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract data from STUDY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data{1,1}
if nnz(design11Idx) | nnz(design12Idx)
    if design11Idx == 0
        design11Idx = 1;
    end
    if design12Idx == 0
        design12Idx = 1;
    end
    data{1,1} = STUDY.cluster(currentCluster).erpdata{design11Idx, design12Idx};
else
    data{1,1} = zeros(length(latency),1);
end
% Extract data{2,1}
if nnz(design21Idx) | nnz(design22Idx)
    if design21Idx == 0
        design21Idx = 1;
    end
    if design22Idx == 0
        design22Idx = 1;
    end
    data{2,1} = STUDY.cluster(currentCluster).erpdata{design21Idx, design22Idx};
else
    data{2,1} = zeros(length(latency),1);
end
% Extract data{3,1}
if nnz(design31Idx) | nnz(design32Idx)
    if design31Idx == 0
        design31Idx = 1;
    end
    if design32Idx == 0
        design32Idx = 1;
    end
    data{3,1} = STUDY.cluster(currentCluster).erpdata{design31Idx, design32Idx};
else
    data{3,1} = zeros(length(latency),1);
end
% Extract data{4,1}
if nnz(design41Idx) | nnz(design42Idx)
    if design41Idx == 0
        design41Idx = 1;
    end
    if design42Idx == 0
        design42Idx = 1;
    end
    data{4,1} = STUDY.cluster(currentCluster).erpdata{design41Idx, design42Idx};
else
    data{4,1} = zeros(length(latency),1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low-pass filter the signal %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upperPassbandEdge = str2double(get(handles.lowPassFilterEdit, 'String')); % Hz
if any(upperPassbandEdge)
    for n = 1:4
        tmpData            = double(data{n,1})';
        unfoldededData     = [tmpData(:,end:-1:2) tmpData tmpData(:,end-1:-1:1)];
        tmpFiltData.data   = unfoldededData;
        tmpFiltData.srate  = ALLEEG(1,1).srate;
        tmpFiltData.trials = 1;
        tmpFiltData.event  = [];
        tmpFiltData.pnts   = size(tmpFiltData.data, 2);
        filtorder          = pop_firwsord('hamming', tmpFiltData.srate, upperPassbandEdge/3); % The last argument represents transition band width.
        tmpFiltData_done   = pop_eegfiltnew(tmpFiltData, 0, upperPassbandEdge, filtorder);
        filteredData       = tmpFiltData_done.data(:, size(tmpData,2):size(tmpData,2)*2-1);
        data{n,1}          = filteredData';
    end
else
    disp('No low-pass filter applied.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subtract baseline %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:4
    tmpData = data{n,1};
    baselineLatency = str2num(get(handles.baselineEdit, 'String'));
    baselineIdx = find(latency >= baselineLatency(1) & latency <= baselineLatency(2));
    baselineCorrected = bsxfun(@minus, tmpData, mean(tmpData(baselineIdx,:),1));
    data{n,1} = baselineCorrected;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare plot titles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
design11Items  = cellstr(      get(handles.design11, 'String'));
design11String = design11Items{get(handles.design11, 'Value')};
design12Items  = cellstr(      get(handles.design12, 'String'));
design12String = design12Items{get(handles.design12, 'Value')};
design21Items  = cellstr(      get(handles.design21, 'String'));
design21String = design21Items{get(handles.design21, 'Value')};
design22Items  = cellstr(      get(handles.design22, 'String'));
design22String = design22Items{get(handles.design22, 'Value')};
design31Items  = cellstr(      get(handles.design31, 'String'));
design31String = design31Items{get(handles.design31, 'Value')};
design32Items  = cellstr(      get(handles.design32, 'String'));
design32String = design32Items{get(handles.design32, 'Value')};
design41Items  = cellstr(      get(handles.design41, 'String'));
design41String = design41Items{get(handles.design41, 'Value')};
design42Items  = cellstr(      get(handles.design42, 'String'));
design42String = design42Items{get(handles.design42, 'Value')};
if (strcmp(design11String, '(none)') | strcmp(design11String, 'Not present')) & (strcmp(design12String, '(none)') | strcmp(design12String, 'Not present'))
    titleA = '';
elseif strcmp(design11String, '(none)') | strcmp(design11String, 'Not present')
    titleA = design12String;
elseif strcmp(design12String, '(none)') | strcmp(design12String, 'Not present')
    titleA = design11String;
else
    titleA = [design11String '-' design12String];
end
if (strcmp(design21String, '(none)') | strcmp(design21String, 'Not present')) & (strcmp(design22String, '(none)') | strcmp(design22String, 'Not present'))
    titleB = '';
elseif strcmp(design21String, '(none)') | strcmp(design21String, 'Not present')
    titleB = design22String;
elseif strcmp(design22String, '(none)') | strcmp(design22String, 'Not present')
    titleB = design21String;
else
    titleB = [design21String '-' design22String];
end
if (strcmp(design31String, '(none)') | strcmp(design31String, 'Not present')) & (strcmp(design32String, '(none)') | strcmp(design32String, 'Not present'))
    titleC = '';
elseif strcmp(design31String, '(none)') | strcmp(design31String, 'Not present')
    titleC = design32String;
elseif strcmp(design32String, '(none)') | strcmp(design32String, 'Not present')
    titleC = design31String;
else
    titleC = [design31String '-' design32String];
end
if (strcmp(design41String, '(none)') | strcmp(design41String, 'Not present')) & (strcmp(design42String, '(none)') | strcmp(design42String, 'Not present'))
    titleD = '';
elseif strcmp(design41String, '(none)') | strcmp(design41String, 'Not present')
    titleD = design42String;
elseif strcmp(design42String, '(none)') | strcmp(design42String, 'Not present')
    titleD = design41String;
else
    titleD = [design41String '-' design42String];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain pairing status %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
designMatrix = [design11Idx design12Idx; design21Idx design22Idx; design31Idx design32Idx; design41Idx design42Idx];
pairingCheckNeeds1 = ~(designMatrix(1,1)==designMatrix(2,1));
pairingCheckNeeds2 = ~(designMatrix(1,2)==designMatrix(2,2));
pairingCheckNeeds3 = ~(designMatrix(3,1)==designMatrix(4,1));
pairingCheckNeeds4 = ~(designMatrix(3,2)==designMatrix(4,2));
pairingCheckNeeds5 = ~(designMatrix(1,1)==designMatrix(3,1));
pairingCheckNeeds6 = ~(designMatrix(1,2)==designMatrix(3,2));

pairingStatus = ones(1,6); % [AB's var1, var2 CD's var1 var2] 
if pairingCheckNeeds1
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,1).pairing, 'off')
        pairingStatus(1,1) = 0; % A vs. B var 1, paired or not
    end
end
if pairingCheckNeeds2
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,2).pairing, 'off')
        pairingStatus(1,2) = 0; % A vs. B var 2, paired or not
    end
end
if pairingCheckNeeds3
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,1).pairing, 'off')
        pairingStatus(1,3) = 0; % C vs. D var 1, paired or not
    end
end
if pairingCheckNeeds4
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,2).pairing, 'off')
        pairingStatus(1,4) = 0; % C vs. D var 2, paired or not
    end
end
if pairingCheckNeeds5
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,1).pairing, 'off')
        pairingStatus(1,5) = 0; % A vs. C var 1, paired or not
    end
end
if pairingCheckNeeds6
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,2).pairing, 'off')
        pairingStatus(1,6) = 0; % A vs. C var 2, paired or not
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate means and errors %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data1Mean = mean(data{1,1},2);
data1Err  = std(data{1,1},0,2)/sqrt(size(data{1,1},2));
data2Mean = mean(data{2,1},2);
data2Err  = std(data{2,1},0,2)/sqrt(size(data{2,1},2));
data3Mean = mean(data{3,1},2);
data3Err  = std(data{3,1},0,2)/sqrt(size(data{3,1},2));
data4Mean = mean(data{4,1},2);
data4Err  = std(data{4,1},0,2)/sqrt(size(data{4,1},2));
diffMean  = (data1Mean-data2Mean)-(data3Mean-data4Mean);
% diffErr will be defined for A-C only condition.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract statsWindow latency in millisecond %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
statsWindowLatency = str2num(get(handles.statsWindowEdit, 'String'));
statsWindowIdx     = find(latency >= statsWindowLatency(1) & latency <= statsWindowLatency(2));
statsWindowColor   = [0.95 0.95 0.85];
fillTimes  = latency(statsWindowIdx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select how to build surrogate statistics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.permOrBootstrapPopupmenu, 'Value') == 1
    surroMethod = 'perm';
else
    surroMethod = 'bootstrap';
end



%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform 2x2 ANOVA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if any(data2Mean) | any(data4Mean)
    
    % Calculate A-B Error
    if nnz([pairingStatus(1) pairingStatus(2)])==2
        diffErr1 = std(data{1,1}-data{2,1},0,2)/sqrt(size(data{1,1}-data{2,1},2));
    else
        % From ttest2_cell() line 90: inhomogenous variance from Howel, 2009, "Statistical Methods for Psychology"
        diffErr1 = sqrt(var(data{1,1},0,2)./size(data{1,1},2) + var(data{2,1},0,2)./size(data{2,1},2));
    end
    
    % Calculate C-D Error
    if nnz([pairingStatus(3) pairingStatus(4)])==2
        diffErr2 = std(data{3,1}-data{4,1},0,2)/sqrt(size(data{3,1}-data{4,1},2));
    else
        % From ttest2_cell() line 90: inhomogenous variance from Howel, 2009, "Statistical Methods for Psychology"
        diffErr2 = sqrt(var(data{3,1},0,2)./size(data{3,1},2) + var(data{4,1},0,2)./size(data{4,1},2));
    end
    
    % Calculate (A-B)-(C-D) Error
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,1).pairing, 'on') & strcmp(STUDY.design(STUDY.currentdesign).variable(1,2).pairing, 'on')
        diffErr3 = std((data{1,1}-data{2,1})-(data{3,1}-data{4,1}),0,2)/sqrt(size((data{1,1}-data{2,1})-(data{3,1}-data{4,1}),2));
    else
        % I don't know how to calculate this... perform geometic mean for time being.
        diffErr3 = sqrt(sqrt(var(data{1,1},0,2)./size(data{1,1},2) + var(data{2,1},0,2)./size(data{2,1},2)) .* sqrt(var(data{3,1},0,2)./size(data{3,1},2) + var(data{4,1},0,2)./size(data{4,1},2)));
    end
    
    % Compute min and max
    upperEnv1 = data1Mean-data2Mean+diffErr1;
    upperEnv2 = data3Mean-data4Mean+diffErr2;
    upperEnv3 = diffMean           +diffErr3;
    lowerEnv1 = data1Mean-data2Mean-diffErr1;
    lowerEnv2 = data3Mean-data4Mean-diffErr2;
    lowerEnv3 = diffMean           -diffErr3;
    plotMinMax = [min([lowerEnv1' lowerEnv2' lowerEnv3']) max([upperEnv1' upperEnv2' upperEnv3'])];
    plotMinMax = [plotMinMax(1)-sum(abs(plotMinMax))/20 plotMinMax(2)+sum(abs(plotMinMax))/20];
    
    % Condition 1
    axes(handles.condition1Axes);
    linePlotHandle = plot(latency, data1Mean-data2Mean, 'Color', [0,0,1], 'LineWidth', 1);
    hold on
    fill([latency latency(end:-1:1)], [upperEnv1' lowerEnv1(end:-1:1)'], [0.8 0.9 1], 'Linestyle', 'none')
    title([titleA ' - ' titleB], 'interpreter', 'none', 'fontsize', 12)
    xlim([latency(1) latency(end)]); ylim(plotMinMax);
    line(xlim, [0 0], 'Color', [0.75 0.75 0.75])
    line([0 0], ylim, 'Color', [0.75 0.75 0.75])
    fillHandle = fill([fillTimes fillTimes(end:-1:1)], [repmat(plotMinMax(1), length(fillTimes)) repmat(plotMinMax(2), length(fillTimes))], statsWindowColor);
    set(fillHandle, 'EdgeColor', statsWindowColor);
    uistack(fillHandle, 'bottom');
    uistack(linePlotHandle, 'top');
    hold off

    % Condition 2 
    axes(handles.condition2Axes);
    linePlotHandle = plot(latency, data3Mean-data4Mean, 'Color', [1,0,0], 'LineWidth', 1)
    hold on
    fill([latency latency(end:-1:1)], [upperEnv2' lowerEnv2(end:-1:1)'], [1 0.8 0.9], 'Linestyle', 'none')
    title([titleC ' - ' titleD], 'interpreter', 'none', 'fontsize', 12)
    xlim([latency(1) latency(end)]); ylim(plotMinMax);
    line(xlim, [0 0], 'Color', [0.75 0.75 0.75])
    line([0 0], ylim, 'Color', [0.75 0.75 0.75])
    fillHandle = fill([fillTimes fillTimes(end:-1:1)], [repmat(plotMinMax(1), length(fillTimes)) repmat(plotMinMax(2), length(fillTimes))], statsWindowColor);
    set(fillHandle, 'EdgeColor', statsWindowColor);
    uistack(fillHandle, 'bottom');
    uistack(linePlotHandle, 'top');
    hold off
    
    % Condition 3
    axes(handles.differenceAxes);
    linePlotHandle = plot(latency, diffMean, 'Color', [0 0.75 0.25], 'LineWidth', 1)
    hold on
    fill([latency latency(end:-1:1)], [upperEnv3' lowerEnv3(end:-1:1)'], [0.8 1 0.9], 'Linestyle', 'none')
    title(sprintf('(%s - %s) -\n(%s - %s)', titleA, titleB, titleC, titleD), 'interpreter', 'none', 'fontsize', 12)
    xlim([latency(1) latency(end)]); ylim(plotMinMax);
    line(xlim, [0 0], 'Color', [0.75 0.75 0.75])
    line([0 0], ylim, 'Color', [0.75 0.75 0.75])
    fillHandle = fill([fillTimes fillTimes(end:-1:1)], [repmat(plotMinMax(1), length(fillTimes)) repmat(plotMinMax(2), length(fillTimes))], statsWindowColor);
    set(fillHandle, 'EdgeColor', statsWindowColor);
    uistack(fillHandle, 'bottom');
    uistack(linePlotHandle, 'top');
    hold off

    % Prepare 2x2 data
    data1 = data{1,1};
    data1Window = data1(statsWindowIdx,:);
    data2 = data{2,1};
    data2Window = data2(statsWindowIdx,:);
    data3 = data{3,1};
    data3Window = data3(statsWindowIdx,:);
    data4 = data{4,1};
    data4Window = data4(statsWindowIdx,:);
    if     get(handles.meanOrPeak, 'value') == 1
        data1AmpForStats = mean(data1Window);
        data2AmpForStats = mean(data2Window);
        data3AmpForStats = mean(data3Window);
        data4AmpForStats = mean(data4Window);
        
    else
        % Find positive peak
        if     get(handles.meanOrPeak, 'value') == 2
            [data1AmpForStats, latency1Idx] = max(data1Window);
            [data2AmpForStats, latency2Idx] = max(data2Window);
            [data3AmpForStats, latency3Idx] = max(data3Window);
            [data4AmpForStats, latency4Idx] = max(data4Window);
            
        % Find negative peak 
        elseif get(handles.meanOrPeak, 'value') == 3
            [data1AmpForStats, latency1Idx] = min(data1Window);
            [data2AmpForStats, latency2Idx] = min(data2Window);
            [data3AmpForStats, latency3Idx] = min(data3Window);
            [data4AmpForStats, latency4Idx] = min(data4Window);
        end
        
        data1LatencyForStats = latency(statsWindowIdx(latency1Idx));
        data2LatencyForStats = latency(statsWindowIdx(latency2Idx));
        data3LatencyForStats = latency(statsWindowIdx(latency3Idx));
        data4LatencyForStats = latency(statsWindowIdx(latency4Idx));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Simple effect A-B %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % For simple effect, perform t-test
    if nnz([pairingStatus(1) pairingStatus(2)])==2 % paired test
        [F_amp, df_amp, pval_amp, surro_amp] = statcond({data1AmpForStats data2AmpForStats},         'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'on');
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data1LatencyForStats data2LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'on');
        end
    else % unpaired test
        [F_amp, df_amp, pval_amp, surog_amp] = statcond({data1AmpForStats data2AmpForStats},         'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'off');
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data1LatencyForStats data2LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'off');
        end
    end
    
    % Generate a report.
    if get(handles.meanOrPeak, 'value') == 1
        str11 = sprintf('Pertumation Test on Mean Ampitude\n');
        str12 = sprintf('     Simple Effect A-B: A = %.3f (%.3f), B = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data1AmpForStats), std(data1AmpForStats), mean(data2AmpForStats), std(data2AmpForStats), df_amp, F_amp, pval_amp);
    else
        str11 = sprintf('Pertumation Test on Peak Ampitude\n');
        str12 = sprintf('     Simple Effect A-B: A = %.3f (%.3f), B = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data1AmpForStats), std(data1AmpForStats), mean(data2AmpForStats), std(data2AmpForStats), df_amp, F_amp, pval_amp);
        str13 = sprintf('Pertumation Test on Peak latency\n');
        str14 = sprintf('     Simple Effect A-B: A = %.3f (%.3f), B = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data1LatencyForStats), std(data1LatencyForStats), mean(data2LatencyForStats), std(data2LatencyForStats), df_lat, F_lat, pval_lat);
%         str11 = sprintf('Pertumation Test on Peak Ampitude\n');
%         str12 = sprintf('     Simple Effect A-B: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n', mean(data1AmpForStats)-mean(data2AmpForStats), df_amp, F_amp, pval_amp);
%         str13 = sprintf('Pertumation Test on Peak latency\n');
%         str14 = sprintf('     Simple Effect A-B: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n', mean(data1LatencyForStats)-mean(data2LatencyForStats), df_lat, F_lat, pval_lat);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Simple effect C-D %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % For simple effect, perform t-test
    if nnz([pairingStatus(3) pairingStatus(4)])==2 % paired test
        [F_amp, df_amp, pval_amp, surro_amp] = statcond({data3AmpForStats data4AmpForStats},         'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'on');
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data3LatencyForStats data4LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'on');
        end
    else % unpaired test
        [F_amp, df_amp, pval_amp, surog_amp] = statcond({data3AmpForStats data4AmpForStats},         'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'off');
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data3LatencyForStats data4LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'off');
        end
    end
    
    % Generate a report.
    if get(handles.meanOrPeak, 'value') == 1
        str21 = sprintf('     Simple Effect C-D: C = %.3f (%.3f), D = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data3AmpForStats), std(data3AmpForStats), mean(data4AmpForStats), std(data4AmpForStats), df_amp, F_amp, pval_amp);
       %str21 = sprintf('     Simple Effect C-D: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n',                  mean(data3AmpForStats)-mean(data4AmpForStats), df_amp, F_amp, pval_amp);
    else
        str21 = sprintf('     Simple Effect C-D: C = %.3f (%.3f), D = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data3AmpForStats), std(data3AmpForStats), mean(data4AmpForStats), std(data4AmpForStats), df_amp, F_amp, pval_amp);
        str22 = sprintf('     Simple Effect C-D: C = %.3f (%.3f), D = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data3LatencyForStats), std(data3LatencyForStats), mean(data4LatencyForStats), std(data4LatencyForStats), df_lat, F_lat, pval_lat);
       %str21 = sprintf('     Simple Effect C-D: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n', mean(data3AmpForStats)-mean(data4AmpForStats), df_amp, F_amp, pval_amp);
       %str22 = sprintf('     Simple Effect C-D: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n', mean(data3LatencyForStats)-mean(data4LatencyForStats), df_lat, F_lat, pval_lat);
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Interaction (A-B)-(C-D) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(STUDY.design(STUDY.currentdesign).variable(1,1).pairing, 'on') & strcmp(STUDY.design(STUDY.currentdesign).variable(1,2).pairing, 'on')
            [F_amp, df_amp, pval_amp, surro_amp] = statcond({data1AmpForStats data2AmpForStats; data3AmpForStats data4AmpForStats},                 'method', surroMethod, 'naccu', 100000, 'paired', 'on'); 
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data1LatencyForStats data2LatencyForStats; data3LatencyForStats data4LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'paired', 'on');
        end
    else
            [F_amp, df_amp, pval_amp, surro_amp] = statcond({data1AmpForStats data2AmpForStats; data3AmpForStats data4AmpForStats},                 'method', surroMethod, 'naccu', 100000, 'paired', 'off'); 
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data1LatencyForStats data2LatencyForStats; data3LatencyForStats data4LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'paired', 'off');
        end
    end
    
    % Generate a report.
    if get(handles.meanOrPeak, 'value') == 1
        % From: https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes,_equal_variance
        equalVarSd12 = sqrt(((length(data1AmpForStats)-1)*var(data1AmpForStats)+(length(data2AmpForStats)-1)*var(data2AmpForStats))/(length(data1AmpForStats)+length(data2AmpForStats)-2));
        equalVarSd34 = sqrt(((length(data3AmpForStats)-1)*var(data3AmpForStats)+(length(data4AmpForStats)-1)*var(data4AmpForStats))/(length(data3AmpForStats)+length(data4AmpForStats)-2));
        str31 = sprintf('     Interaction: A-B = %.3f (%.3f), C-D = %.3f (%.3f), F(%.0f,%.0f) = %.2f, p = %.4f', mean(data1AmpForStats)-mean(data2AmpForStats), equalVarSd12, mean(data3AmpForStats)-mean(data4AmpForStats), equalVarSd34, df_amp{1,3}(1,1), df_amp{1,3}(1,2), F_amp{1,3}, pval_amp{1,3});
        %str31 = sprintf('     Interaction: M_diff = %.2f, F(%.0f,%.0f) = %.2f, p = %.4f', (mean(data1AmpForStats)-mean(data2AmpForStats))-(mean(data3AmpForStats)-mean(data4AmpForStats)), df_amp{1,3}(1,1), df_amp{1,3}(1,2), F_amp{1,3}, pval_amp{1,3});
        str = [str11 str12 str21 str31];
    else
        equalVarSd12 = sqrt(((length(data1AmpForStats)-1)*var(data1AmpForStats)+(length(data2AmpForStats)-1)*var(data2AmpForStats))/(length(data1AmpForStats)+length(data2AmpForStats)-2));
        equalVarSd34 = sqrt(((length(data3AmpForStats)-1)*var(data3AmpForStats)+(length(data4AmpForStats)-1)*var(data4AmpForStats))/(length(data3AmpForStats)+length(data4AmpForStats)-2));
        str31 = sprintf('     Interaction: A-B = %.3f (%.3f), C-D = %.3f (%.3f), F(%.0f,%.0f) = %.2f, p = %.4f\n', mean(data1AmpForStats)-mean(data2AmpForStats), equalVarSd12, mean(data3AmpForStats)-mean(data4AmpForStats), equalVarSd34, df_amp{1,3}(1,1), df_amp{1,3}(1,2), F_amp{1,3}, pval_amp{1,3});
        
        equalVarSd12 = sqrt(((length(data1LatencyForStats)-1)*var(data1LatencyForStats)+(length(data2LatencyForStats)-1)*var(data2LatencyForStats))/(length(data1LatencyForStats)+length(data2LatencyForStats)-2));
        equalVarSd34 = sqrt(((length(data3LatencyForStats)-1)*var(data3LatencyForStats)+(length(data4LatencyForStats)-1)*var(data4LatencyForStats))/(length(data3LatencyForStats)+length(data4LatencyForStats)-2));
        str32 = sprintf('     Interaction: A-B = %.3f (%.3f), C-D = %.3f (%.3f), F(%.0f,%.0f) = %.2f, p = %.4f\n', mean(data1LatencyForStats)-mean(data2LatencyForStats), equalVarSd12, mean(data3LatencyForStats)-mean(data4LatencyForStats), equalVarSd34, df_lat{1,3}(1,1), df_lat{1,3}(1,2), F_lat{1,3}, pval_lat{1,3});

        str = [str11 str12 str21 str31 str13 str14 str22 str32];
        
%         str31 = sprintf('     Interaction: M_diff = %.2f, F(%.0f,%.0f) = %.2f, p = %.4f\n', (mean(data1AmpForStats)-mean(data2AmpForStats))-(mean(data3AmpForStats)-mean(data4AmpForStats)), df_amp{1,3}(1,1), df_amp{1,3}(1,2), F_amp{1,3}, pval_amp{1,3});
%         str32 = sprintf('     Interaction: M_diff = %.2f, F(%.0f,%.0f) = %.2f, p = %.4f',   (mean(data1LatencyForStats)-mean(data2LatencyForStats)) - (mean(data3LatencyForStats)-mean(data4LatencyForStats)), df_lat{1,3}(1,1), df_lat{1,3}(1,2), F_lat{1,3}, pval_lat{1,3});
%         str = [str11 str12 str21 str31 str13 str14 str22 str32];
    end

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform paired/two-sample t-test %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    if nnz([pairingStatus(5) pairingStatus(6)])==2
        diffErr = std(data{1,1}-data{3,1},0,2)/sqrt(size(data{1,1}-data{3,1},2));
    else
        % From ttest2_cell() line 90: inhomogenous variance from Howel, 2009, "Statistical Methods for Psychology"
        diffErr = sqrt(var(data{1,1},0,2)./size(data{1,1},2) + var(data{3,1},0,2)./size(data{3,1},2));
    end
    
    % Compute max and min
    plotMinMax = [min([data1Mean-data1Err; data3Mean-data3Err; diffMean-diffErr]) max([data1Mean+data1Err; data3Mean+data3Err; diffMean+diffErr])];
    plotMinMax = [plotMinMax(1)-sum(abs(plotMinMax))/20 plotMinMax(2)+sum(abs(plotMinMax))/20];
    
    % Condition 1
    axes(handles.condition1Axes);
    
    linePlotHandle = plot(latency, data1Mean, 'Color', [0,0,1], 'LineWidth', 2);
    hold on
    upperEnv = data1Mean+data1Err;
    lowerEnv = data1Mean-data1Err;
    fill([latency latency(end:-1:1)], [upperEnv' lowerEnv(end:-1:1)'], [0.8 0.9 1], 'Linestyle', 'none')
    title(titleA, 'interpreter', 'none', 'fontsize', 12)
    xlim([latency(1) latency(end)]); ylim(plotMinMax);
    line(xlim, [0 0], 'Color', [0.75 0.75 0.75])
    line([0 0], ylim, 'Color', [0.75 0.75 0.75])
    fillHandle = fill([fillTimes fillTimes(end:-1:1)], [repmat(plotMinMax(1), length(fillTimes)) repmat(plotMinMax(2), length(fillTimes))], statsWindowColor);
    set(fillHandle, 'EdgeColor', statsWindowColor);
    uistack(fillHandle, 'bottom');
    uistack(linePlotHandle, 'top');
    hold off
    
    % Condition 2
    axes(handles.condition2Axes);
    linePlotHandle = plot(latency, data3Mean, 'Color', [1,0,0], 'LineWidth', 2);
    hold on
    upperEnv = data3Mean+data3Err;
    lowerEnv = data3Mean-data3Err;
    fill([latency latency(end:-1:1)], [upperEnv' lowerEnv(end:-1:1)'], [1 0.8 0.9], 'Linestyle', 'none')
    title(titleC, 'interpreter', 'none', 'fontsize', 12)
    xlim([latency(1) latency(end)]); ylim(plotMinMax);
    line(xlim, [0 0], 'Color', [0.75 0.75 0.75])
    line([0 0], ylim, 'Color', [0.75 0.75 0.75])
    fillHandle = fill([fillTimes fillTimes(end:-1:1)], [repmat(plotMinMax(1), length(fillTimes)) repmat(plotMinMax(2), length(fillTimes))], statsWindowColor);
    set(fillHandle, 'EdgeColor', statsWindowColor);
    uistack(fillHandle, 'bottom');
    uistack(linePlotHandle, 'top');
    hold off
    
    % Diff
    axes(handles.differenceAxes);
    linePlotHandle = plot(latency, diffMean, 'Color', [0 0.75 0.25], 'LineWidth', 2);
    hold on
    upperEnv = diffMean+diffErr;
    lowerEnv = diffMean-diffErr;
    fill([latency latency(end:-1:1)], [upperEnv' lowerEnv(end:-1:1)'], [0.8 1 0.9], 'Linestyle', 'none')
    title([titleA ' - ' titleC], 'interpreter', 'none', 'fontsize', 12)
    xlim([latency(1) latency(end)]); ylim(plotMinMax);
    line(xlim, [0 0], 'Color', [0.75 0.75 0.75])
    line([0 0], ylim, 'Color', [0.75 0.75 0.75])
    fillHandle = fill([fillTimes fillTimes(end:-1:1)], [repmat(plotMinMax(1), length(fillTimes)) repmat(plotMinMax(2), length(fillTimes))], statsWindowColor);
    set(fillHandle, 'EdgeColor', statsWindowColor);
    uistack(fillHandle, 'bottom');
    uistack(linePlotHandle, 'top');
    hold off
    
    % Prepare 1x2 data
    data1 = data{1,1};
    data1Window = data1(statsWindowIdx,:);
    data3 = data{3,1};
    data3Window = data3(statsWindowIdx,:);
    if     get(handles.meanOrPeak, 'value') == 1
        data1AmpForStats = mean(data1Window);
        data3AmpForStats = mean(data3Window);
    else
        
        if     get(handles.meanOrPeak, 'value') == 2
            [data1AmpForStats, latency1Idx] = max(data1Window);
            [data3AmpForStats, latency3Idx] = max(data3Window);
            
            % Find negative peak
        elseif get(handles.meanOrPeak, 'value') == 3
            [data1AmpForStats, latency1Idx] = min(data1Window);
            [data3AmpForStats, latency3Idx] = min(data3Window);
        end
        data1LatencyForStats = latency(statsWindowIdx(latency1Idx));
        data3LatencyForStats = latency(statsWindowIdx(latency3Idx));
    end

    % Perform t-test
    if nnz([pairingStatus(5) pairingStatus(6)])==2 % paired test
            [F_amp, df_amp, pval_amp, surro_amp] = statcond({data1AmpForStats data3AmpForStats},         'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'on');
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data1LatencyForStats data3LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'on');
        end
    else % unpaired test
            [F_amp, df_amp, pval_amp, surog_amp] = statcond({data1AmpForStats data3AmpForStats},         'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'off');
        if get(handles.meanOrPeak, 'value') == 2 | get(handles.meanOrPeak, 'value') == 3
            [F_lat, df_lat, pval_lat, surro_lat] = statcond({data1LatencyForStats data3LatencyForStats}, 'method', surroMethod, 'naccu', 100000, 'tail', 'both', 'paired', 'off');
        end    
    end
       
%     % Generate a report.
%     if get(handles.meanOrPeak, 'value') == 1
%         str21 = sprintf('     Simple Effect C-D: C = %.3f (%.3f), D = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data3AmpForStats), std(data3AmpForStats), mean(data4AmpForStats), std(data4AmpForStats), df_amp, F_amp, pval_amp);
%         %str21 = sprintf('     Simple Effect C-D: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n',                  mean(data3AmpForStats)-mean(data4AmpForStats), df_amp, F_amp, pval_amp);
%     else
%         str21 = sprintf('     Simple Effect C-D: C = %.3f (%.3f), D = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data3AmpForStats), std(data3AmpForStats), mean(data4AmpForStats), std(data4AmpForStats), df_amp, F_amp, pval_amp);
%         str22 = sprintf('     Simple Effect C-D: C = %.3f (%.3f), D = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data3LatencyForStats), std(data3LatencyForStats), mean(data4LatencyForStats), std(data4LatencyForStats), df_lat, F_lat, pval_lat);
%         %str21 = sprintf('     Simple Effect C-D: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n', mean(data3AmpForStats)-mean(data4AmpForStats), df_amp, F_amp, pval_amp);
%         %str22 = sprintf('     Simple Effect C-D: M_diff = %.2f, t(%.1f) = %.2f, p = %.4f\n', mean(data3LatencyForStats)-mean(data4LatencyForStats), df_lat, F_lat, pval_lat);
%     end
    
    % Generate a report.
    if get(handles.meanOrPeak, 'value') == 1
        str1 = sprintf('Pertumation Test on Mean Ampitude\n');
        str2 = sprintf('     Main Effect A-C: A = %.3f (%.3f), C = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data1AmpForStats), std(data1AmpForStats), mean(data3AmpForStats), std(data3AmpForStats), df_amp, F_amp, pval_amp);
       %str2 = sprintf('     Main Effect A-C: M (SD) = %.2f (%.2f) vs. %.2f (%.2f), t(%.1f) = %.2f, p = %.4f\n', mean(data1AmpForStats), std(data1AmpForStats), mean(data3AmpForStats), std(data3AmpForStats), df_amp, F_amp, pval_amp);
        str = [str1 str2];
    else
        str1 = sprintf('Pertumation Test on Peak Ampitude\n');
        str2 = sprintf('     Main Effect A-C: A = %.3f (%.3f), C = %.3f (%.3f), t(%.1f) = %.2f, p = %.4f\n', mean(data1AmpForStats), std(data1AmpForStats), mean(data3AmpForStats), std(data3AmpForStats), df_amp, F_amp, pval_amp);
        str3 = sprintf('Pertumation Test on Peak latency\n');
        str4 = sprintf('     Main Effect A-C: A = %.0f (%.0f), C = %.0f (%.0f), t(%.1f) = %.2f, p = %.4f\n', mean(data1LatencyForStats), std(data1LatencyForStats), mean(data3LatencyForStats), std(data3LatencyForStats), df_lat, F_lat, pval_lat);
       %str4 = sprintf('     Main Effect A-C: M (SD) = %.0f (%.2f) vs. %.2f (%.2f), t(%.1f) = %.2f, p = %.4f\n', mean(data1LatencyForStats), std(data1LatencyForStats), mean(data3LatencyForStats), std(data3LatencyForStats), df_lat, F_lat, pval_lat);
        str = [str1 str2 str3 str4];
    end
    if nnz([pairingStatus(5) pairingStatus(6)])~=2
        str = sprintf([str '     Degrees of Freedom was adjusted for inhomogenous variance.']);
    end
end

% Display the report.
set(handles.resultReportText, 'String', str)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract subject info %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A
if any(design11Idx) & any(design12Idx)
    setinds = STUDY.cluster(currentCluster).setinds{design11Idx, design12Idx};
    datasetIdx = cell2mat({STUDY.design(STUDY.currentdesign).cell(setinds).dataset})';
    subjectNames = {ALLEEG(1,datasetIdx).subject}';
    if get(handles.meanOrPeak, 'value') == 1
        printIndividualSubjectA = [subjectNames num2cell(data1AmpForStats')]
    else
        printIndividualSubjectA = [subjectNames num2cell(data1AmpForStats') num2cell(data1LatencyForStats')]
    end
else
    printIndividualSubjectA = {};
end
% B
if any(design21Idx) & any(design22Idx)
    setinds = STUDY.cluster(currentCluster).setinds{design21Idx, design22Idx};
    datasetIdx = cell2mat({STUDY.design(STUDY.currentdesign).cell(setinds).dataset})';
    subjectNames = {ALLEEG(1,datasetIdx).subject}';
    if get(handles.meanOrPeak, 'value') == 1
        printIndividualSubjectB = [subjectNames num2cell(data2AmpForStats')]
    else
        printIndividualSubjectB = [subjectNames num2cell(data2AmpForStats') num2cell(data2LatencyForStats')]
    end
else
    printIndividualSubjectB = {};
end
% C
if any(design31Idx) & any(design32Idx)
    setinds = STUDY.cluster(currentCluster).setinds{design31Idx, design32Idx};
    datasetIdx = cell2mat({STUDY.design(STUDY.currentdesign).cell(setinds).dataset})';
    subjectNames = {ALLEEG(1,datasetIdx).subject}';
    if get(handles.meanOrPeak, 'value') == 1
        printIndividualSubjectC = [subjectNames num2cell(data3AmpForStats')]
    else
        printIndividualSubjectC = [subjectNames num2cell(data3AmpForStats') num2cell(data3LatencyForStats')]
    end
else
    printIndividualSubjectC = {};
end
% D
if any(design41Idx) & any(design42Idx)
    setinds = STUDY.cluster(currentCluster).setinds{design41Idx, design42Idx};
    datasetIdx = cell2mat({STUDY.design(STUDY.currentdesign).cell(setinds).dataset})';
    subjectNames = {ALLEEG(1,datasetIdx).subject}';
    if get(handles.meanOrPeak, 'value') == 1
        printIndividualSubjectD = [subjectNames num2cell(data4AmpForStats')]
    else
        printIndividualSubjectD = [subjectNames num2cell(data4AmpForStats') num2cell(data4LatencyForStats')]
    end
else
    printIndividualSubjectD = {};
end

% Store individual subjects' results under STUDY.etc.erpStudio
STUDY.etc.erpStudio.individualSubjectsA = printIndividualSubjectA;
STUDY.etc.erpStudio.individualSubjectsB = printIndividualSubjectB;
STUDY.etc.erpStudio.individualSubjectsC = printIndividualSubjectC;
STUDY.etc.erpStudio.individualSubjectsD = printIndividualSubjectD;
assignin('base', 'STUDY', STUDY);
disp('Individual subjects'' data are stored under STUDY.etc.erpStudio.')

% Show how many unique subjects and number of ICs par group
numUniqueSubj1 = length(unique(STUDY.cluster(currentCluster).setinds{design11Idx, design12Idx}));
numUniqueSubj3 = length(unique(STUDY.cluster(currentCluster).setinds{design31Idx, design32Idx}));
numICs1 = length(data1AmpForStats);
numICs3 = length(data3AmpForStats);
str1 = sprintf('A: %.0f unique subjects, %.0f ICs.\n', numUniqueSubj1, numICs1);
str3 = sprintf('C: %.0f unique subjects, %.0f ICs.\n', numUniqueSubj3, numICs3);

if any(data2Mean) | any(data4Mean)
    numUniqueSubj2 = length(unique(STUDY.cluster(currentCluster).setinds{design21Idx, design22Idx}));
    numUniqueSubj4 = length(unique(STUDY.cluster(currentCluster).setinds{design41Idx, design42Idx}));
    numICs2 = length(data2AmpForStats);
    numICs4 = length(data4AmpForStats);
    str2 = sprintf('B: %.0f unique subjects, %.0f ICs.\n', numUniqueSubj2, numICs2);
    str4 = sprintf('D: %.0f unique subjects, %.0f ICs.', numUniqueSubj4, numICs4);
    str = [str1 str2 str3 str4];
else
    str = [str1 str3 ];
end
set(handles.uniqueSubjectReport, 'String', str);

% Change figure title
set(gcf, 'Name', 'std_erpStudio()')

disp('Done.')

% Choose default command line output for std_erpStudio
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Executes on selection change in design22.
function design22_Callback(hObject, eventdata, handles)
% hObject    handle to design22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design22 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design22


% --- Executes during object creation, after setting all properties.
function design22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in design21.
function design21_Callback(hObject, eventdata, handles)
% hObject    handle to design21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design21 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design21


% --- Executes during object creation, after setting all properties.
function design21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in design12.
function design12_Callback(hObject, eventdata, handles)
% hObject    handle to design12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design12


% --- Executes during object creation, after setting all properties.
function design12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in design11.
function design11_Callback(hObject, eventdata, handles)
% hObject    handle to design11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design11


% --- Executes during object creation, after setting all properties.
function design11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baselineEdit_Callback(hObject, eventdata, handles)
% hObject    handle to baselineEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baselineEdit as text
%        str2double(get(hObject,'String')) returns contents of baselineEdit as a double


% --- Executes during object creation, after setting all properties.
function baselineEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baselineEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function clusterList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

clusterList = evalin('base', 'STUDY.cluster');
allClusterName = {clusterList.name}';
allClusterName = allClusterName(2:end);
set(hObject, 'String', allClusterName);

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in selectStatisticsMethod.
function selectStatisticsMethod_Callback(hObject, eventdata, handles)
% hObject    handle to selectStatisticsMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectStatisticsMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectStatisticsMethod


% --- Executes during object creation, after setting all properties.
function selectStatisticsMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectStatisticsMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makeBackgroundWhitecheckbox.
function makeBackgroundWhitecheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to makeBackgroundWhitecheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeBackgroundWhitecheckbox
if get(hObject,'Value')==1
    set(gcf, 'color', [1 1 1]);
    set([handles.text23 ...
         handles.text24 ...
         handles.text25 ...
         handles.text26 ...
         handles.text27 ...
         handles.text28 ...
         handles.text18 ...
         handles.text19 ...
         handles.resultReportText ...
         handles.makeBackgroundWhitecheckbox ...
         handles.uniqueSubjectReport], 'BackgroundColor', [1 1 1])
else
    set(gcf, 'color', [0.66 0.76 1]);
    set([handles.text23 ...
         handles.text24 ...
         handles.text25 ...
         handles.text26 ...
         handles.text27 ...
         handles.text28 ...
         handles.text18 ...
         handles.text19 ...
         handles.resultReportText ...
         handles.makeBackgroundWhitecheckbox ...
         handles.uniqueSubjectReport], 'BackgroundColor', [0.66 0.76 1])
end



function statsWindowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to statsWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of statsWindowEdit as text
%        str2double(get(hObject,'String')) returns contents of statsWindowEdit as a double


% --- Executes during object creation, after setting all properties.
function statsWindowEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statsWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in design42.
function design42_Callback(hObject, eventdata, handles)
% hObject    handle to design42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design42 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design42


% --- Executes during object creation, after setting all properties.
function design42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in design41.
function design41_Callback(hObject, eventdata, handles)
% hObject    handle to design41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design41 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design41


% --- Executes during object creation, after setting all properties.
function design41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in design32.
function design32_Callback(hObject, eventdata, handles)
% hObject    handle to design32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design32 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design32


% --- Executes during object creation, after setting all properties.
function design32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in design31.
function design31_Callback(hObject, eventdata, handles)
% hObject    handle to design31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns design31 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design31


% --- Executes during object creation, after setting all properties.
function design31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowPassFilterEdit_Callback(hObject, eventdata, handles)
% hObject    handle to lowPassFilterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowPassFilterEdit as text
%        str2double(get(hObject,'String')) returns contents of lowPassFilterEdit as a double


% --- Executes during object creation, after setting all properties.
function lowPassFilterEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowPassFilterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in meanOrPeak.
function meanOrPeak_Callback(hObject, eventdata, handles)
% hObject    handle to meanOrPeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns meanOrPeak contents as cell array
%        contents{get(hObject,'Value')} returns selected item from meanOrPeak


% --- Executes during object creation, after setting all properties.
function meanOrPeak_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanOrPeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in permOrBootstrapPopupmenu.
function permOrBootstrapPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to permOrBootstrapPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns permOrBootstrapPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from permOrBootstrapPopupmenu


% --- Executes during object creation, after setting all properties.
function permOrBootstrapPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to permOrBootstrapPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
