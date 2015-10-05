function varargout = GUI_Analyse_MD(varargin)
% GUI_ANALYSE_MD M-file for GUI_Analyse_MD.fig
%      GUI_ANALYSE_MD, by itself, creates a new GUI_ANALYSE_MD or raises the existing
%      singleton*.
%
%      H = GUI_ANALYSE_MD returns the handle to a new GUI_ANALYSE_MD or the handle to
%      the existing singleton*.
%
%      GUI_ANALYSE_MD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ANALYSE_MD.M with the given input arguments.
%
%      GUI_ANALYSE_MD('Property','Value',...) creates a new GUI_ANALYSE_MD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Analyse_MD_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Analyse_MD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Analyse_MD

% Last Modified by GUIDE v2.5 21-Jan-2010 10:16:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Analyse_MD_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Analyse_MD_OutputFcn, ...
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


% --- Executes just before GUI_Analyse_MD is made visible.
function GUI_Analyse_MD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Analyse_MD (see VARARGIN)

GUI_InitializeData( hObject );
handles = guidata(hObject);

GUI_InitializeControls( hObject );
handles = guidata(hObject);

% UIWAIT makes GUI_Analyse_MD wait for user response (see UIRESUME)
% uiwait(handles.figure1);

DrawCoordAxes( hObject );

return;

function GUI_InitializeData( hObject )

global Data X;

handles = guidata(hObject);

% Choose default command line output for GUI_Analyse_MD
handles.output = hObject;
handles.Data = Data;
handles.X = X;
handles.NumberOfDataPoints = size(X,1);
handles.RestrictIdxs = [];
handles.CoordAxes_Pts = [];
handles.CoordAxes_Pts2= [];
if size(handles.Data.Coords,2)>=3,
    handles.DisplayCoordsIdxs = [1,2,3];
else
    handles.DisplayCoordsIdxs = [1,2];
end;
handles.DisplayCoords = handles.Data.Coords(:,handles.DisplayCoordsIdxs); %handles.Data.G.EigenVecs(:,handles.DisplayCoordsIdxs);
handles.CoordAxes_SelPtIdxs = [];
handles.CoordAxes_SelPtHandles = [];
handles.CoordAxes_NNInfo = [];
handles.CoordAxes_PtColor = 1;
handles.CoordAxesColorMode = 0;
handles.CoordAxes2D = true;

if ~isfield(handles.Data,'Stats'),      handles.Data.Stats = [];    end;
if ~isfield(handles.Data,'Nets'),       handles.Data.Nets = [];     end;
if ~isfield(handles.Data.Stats,'S'),    handles.Data.Stats.S = [];  end;

% Update handles structure
guidata(hObject, handles);

return;


function GUI_InitializeControls( hObject )

handles = guidata(hObject);

% Color bar for CoordAxes
axes(handles.CoordAxes);


% Strings for Color eigenfunctions  popup menu
lStrings = [];
for k= 1:size(handles.Data.G.EigenVecs,2),
    lStrings{k} = sprintf('Eigfcn %d',k);
end;

set(handles.pmnuColorEigenfunctions,'String',lStrings);
set(handles.pmnuColorEigenfunctions,'Value',4);

guidata( hObject, handles);

pmnuColorEigenfunctions_Callback(handles.pmnuColorEigenfunctions);

handles = guidata( hObject );

% Strings for Color SVD, NSVD, JFr, Res popup menu
idx = 1;
lStrings = [];
lStringsN = [];
lStringsJ = [];
lStringsL = [];
lStringsResVar = [];
lStringsNResVar = [];
lUserData = [];
for k = 1:size(handles.Data.Stats.S,3),
    for j = 1:size(handles.Data.Stats.S,1),
        lStrings{idx} = sprintf('SVD(%.3f,%d)',handles.Data.NetsOpts.Delta(j),k);
        lStringsResVar{idx} = sprintf('Res(%.3f,%d)',handles.Data.NetsOpts.Delta(j),k);
        lStringsNResVar{idx} = sprintf('NRes(%.3f,%d)',handles.Data.NetsOpts.Delta(j),k);
        lUserData(idx,:) = [j,k];
        idx = idx+1;
    end;
end;
idx = 1;
for k = 2:size(handles.Data.Stats.S,3),
    for j = 1:size(handles.Data.Stats.S,1),
        lStringsN{idx} = sprintf('MSVD(%.3f,%d)',handles.Data.NetsOpts.Delta(j),k);
        lUserData(idx,:) = [j,k];
        idx = idx+1;
    end;
end;
%for i = 1:size(handles.Data.Stats.JFr,2),
%    lStringsJ{i} = sprintf('JFr(%d)',i);
%end;
if isfield(handles.Data,'Labels'),
    for i = 1:size(handles.Data.Labels,2),
        lStringsL{i} = sprintf('Lbl %d',i);
    end;
end;
set(handles.pmnuColorSVD,'String',lStrings);
set(handles.pmnuColorSVD,'UserData',lUserData);
set(handles.pmnuColorNSVD,'String',lStringsN);
set(handles.pmnuColorNSVD,'UserData',lUserData);
set(handles.pmnuColorResVar,'String',lStringsResVar);
set(handles.pmnuColorResVar,'UserData',lUserData);
set(handles.pmnuColorNResVar,'String',lStringsNResVar);
set(handles.pmnuColorNResVar,'UserData',lUserData);
set(handles.pmnuLblCoordAxes,'UserData',lUserData);
%set(handles.pmnuColorJFr,'String',lStringsJ);
if ~isempty(lStringsL),
    set(handles.pmnuLblCoordAxes,'String',lStringsL);
else
    set(handles.pmnuLblCoordAxes,'String','Empty');
end;

lStrings = [];
for k = 1:length(handles.Data.Nets),
    lStrings{k} = sprintf('Net(%d).InvIdxs',k);
end;
set(handles.pmnuColorNetInvIdxs,'String',lStrings);

lStrings = [];
for k = 1:length(handles.Data.Nets),
    lStrings{k} = sprintf('Volume (%d)',k);
end;
set(handles.pmnuVolume,'String',lStrings);

lStrings = [];
lStrings{1} = 'All';
for k = 1000:1000:handles.NumberOfDataPoints,
    lStrings{k/1000+1} = sprintf('%d',k);
end;
set(handles.pmnuRestrictToIndices,'String',lStrings);


%
%
% Initialize data for second and third axes
%
%
handles.Axes{1}.handle = handles.MultiscaleSVDAxes;
handles.Axes{1}.Data{1}.X = handles.Data.NetsOpts.Delta(1:handles.Data.NetsOpts.NumberOfScales);
handles.Axes{1}.Data{1}.Xfixdim = [];
handles.Axes{1}.Data{1}.Xfixparam = [];
handles.Axes{1}.Data{1}.Y = handles.Data.Stats.S;
handles.Axes{1}.Data{1}.Yfixdim = 2;
handles.Axes{1}.Data{1}.Yfixparam = 1;
handles.Axes{1}.Data{1}.PlotStyle = 'o:';
try
handles.Axes{1}.GoodScales = handles.Data.EstDimStats.GoodScales;
catch
end;
handles.Axes{1}.Data{2}.X = handles.Data.NetsOpts.Delta(1:handles.Data.NetsOpts.NumberOfScales);
handles.Axes{1}.Data{2}.Xfixdim = [];
handles.Axes{1}.Data{2}.Xfixparam = [];
try
    handles.Axes{1}.Data{2}.Y = handles.Data.Stats.Tychonov.S;
catch
    handles.Axes{1}.Data{2}.Y = handles.Data.Stats.S;
end;
handles.Axes{1}.Data{2}.Yfixdim = 2;
handles.Axes{1}.Data{2}.Yfixparam = 1;
handles.Axes{1}.Data{2}.PlotStyle = 'x-';
handles.Axes{1}.Xlbl = 'Scale';
handles.Axes{1}.Ylbl = 'Singular values';
handles.Axes{1}.title = 'Singular values and their regularization';

handles.Axes{2}.handle = handles.DataAxes;
handles.Axes{2}.Data{1}.X = handles.Data.NetsOpts.Delta(1:handles.Data.NetsOpts.NumberOfScales);
handles.Axes{2}.Data{1}.Xfixdim = [];
handles.Axes{2}.Data{1}.Xfixparam = [];
handles.Axes{2}.Data{1}.Y = handles.Data.Stats.S;
handles.Axes{2}.Data{1}.Yfixdim = 2;
handles.Axes{2}.Data{1}.Yfixparam = 1;
handles.Axes{2}.Data{1}.PlotStyle = 'o:';
handles.Axes{2}.Data{2}.X = handles.Data.NetsOpts.Delta(1:handles.Data.NetsOpts.NumberOfScales);
handles.Axes{2}.Data{2}.Xfixdim = [];
handles.Axes{2}.Data{2}.Xfixparam = [];
try
    handles.Axes{2}.Data{2}.Y = handles.Data.Stats.Tychonov.S;
catch
    handles.Axes{2}.Data{2}.Y = handles.Data.Stats.S;
end;
handles.Axes{2}.Data{2}.Yfixdim = 2;
handles.Axes{2}.Data{2}.Yfixparam = 1;
handles.Axes{2}.Data{2}.PlotStyle = 'x-';
handles.Axes{2}.Xlbl = 'Scale';
handles.Axes{2}.Ylbl = 'Singular values';
handles.Axes{2}.title = 'Singular values and their regularization';


guidata( hObject, handles );

return;


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Analyse_MD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over axes background.
function CoordAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to CoordAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%axes(handles.CoordAxes);
lCurPt = get(handles.CoordAxes,'CurrentPoint');
if handles.CoordAxes2D,
    lCurPt = lCurPt([1,3]);%lXRange = get(gca,'XLim');lYRange=get(gca,'YLim');%lCurPt = lCurPt.*[lXRange(2)-lXRange(1),lYRange(2)-lYRange(1)] + [lXRange(1),lYRange(1)],
    lCurCoordIdxs = 1:2;
    % Find the nearest point to the clicked point
    [count,lClosestPtIdx,dists,handles.CoordAxes_NNInfo] = ...
        nrsearch(handles.DisplayCoords(:,lCurCoordIdxs), lCurPt, 1, 0, [], ...
        struct('ReturnAsArrays',true,'NNInfo',handles.CoordAxes_NNInfo,'FastNNSearcher',[],'Tolerance',0,'XIsTransposed',true));
else
%    lCurPt = mean(lCurPt,1);       % Mean point on the line from front to back
%    lCurPt = lCurPt([1,3,5]);      % Front point on the line from front to back
    lLinePts = [linspace(lCurPt(1,1),lCurPt(2,1),30)',linspace(lCurPt(1,2),lCurPt(2,2),30)',linspace(lCurPt(1,3),lCurPt(2,3),30)'];
    lCurCoordIdxs = 1:3;
    % Find the nearest point to the clicked point
    [count,lClosestPtIdx,dists,handles.CoordAxes_NNInfo] = ...
        nrsearch(handles.DisplayCoords(:,lCurCoordIdxs), lLinePts, 1, 0, [], ...
        struct('ReturnAsArrays',true,'NNInfo',handles.CoordAxes_NNInfo,'FastNNSearcher',[],'Tolerance',0,'XIsTransposed',true));
    [lTmpMinDist,lTmpMinDistIdx] = min(dists);
    lClosestPtIdx = lClosestPtIdx(lTmpMinDistIdx);
end;

% Add closest point to list of clicked points
handles.CoordAxes_SelPtIdxs = [handles.CoordAxes_SelPtIdxs,lClosestPtIdx];

% Update handles structure
guidata( hObject, handles);

%DrawCoordAxes( hObject, eventdata, handles );
DrawCoordAxesSelPt( hObject );

% Draw multiscale svd statistics
%DrawMultiscaleSVDAxes( hObject );

DrawGUIDataAxes( hObject, 1 );
DrawGUIDataAxes( hObject, 2 );

% Draw molecular configuration
DrawMoleculeAxes( hObject );

return;


%
%
% Updates the visualization the selected point in CoordAxes
%
%
function DrawCoordAxesSelPt( hObject )

handles = guidata( hObject );

%axes(handles.CoordAxes);hold on;

% Delete selection of previously selected point,
% display newly selected point.
if ~isempty(handles.CoordAxes_SelPtIdxs),
    % Delete previous point
    try
    delete(handles.CoordAxes_SelPtHandles);
    catch
    end;
    % Plot current point    
    lCoords = handles.DisplayCoords(handles.CoordAxes_SelPtIdxs(length(handles.CoordAxes_SelPtIdxs)),1:min([3,size(handles.DisplayCoords,2)]));
    handles.CoordAxes_SelPtHandles = scatter3(lCoords(1),lCoords(2),lCoords(3),200,0,'filled','Marker','d');        
end;

% Update handles structure
guidata( hObject, handles);

return;


%
%
% Updates the visualization in MoleculeAxes
%
%
function DrawMoleculeAxes( hObject )

try         % MM : not elegant it will do for now.
    handles = guidata( hObject );

    axes(handles.MoleculeAxes);cla;

    idx = handles.CoordAxes_SelPtIdxs(length(handles.CoordAxes_SelPtIdxs));

    if length(handles.CoordAxes_SelPtIdxs)>1,
        previdx = handles.CoordAxes_SelPtIdxs(length(handles.CoordAxes_SelPtIdxs)-1);
        [lTmp,lAlignedConfiguration] = MolecularMetric( handles.Data.X_norm(:,previdx),handles.Data.X_norm(:,idx) );
        % Display the corresponding molecular configuration
        feval( @DisplayMolecule, [], handles.Data.X_norm(:,previdx), handles.MoleculeAxes, struct('Color',0) );hold on;
    else
        lAlignedConfiguration = handles.Data.X_norm(:,idx);
    end;

    feval( @DisplayMolecule, [], lAlignedConfiguration, handles.MoleculeAxes, struct('Color',1) );
catch
end;

try
idx = handles.CoordAxes_SelPtIdxs(length(handles.CoordAxes_SelPtIdxs));
axes(handles.MoleculeAxes);
if isfield(handles.Data,'XCoords') && ~isempty(handles.Data.XCoords),
    lImageSize = sqrt(size(handles.Data.XCoords,1));
    subimage(double(reshape(handles.X(idx,:)*handles.Data.XCoords',[lImageSize,lImageSize])),gray)
    lV = handles.Data.XCoords*squeeze(handles.Data.Nets(5).V(handles.Data.Nets(5).InvIdxs(idx),:,:));
else
    lImageSize = sqrt(size(handles.X,2));
    subimage(double(reshape(handles.X(idx,:),[lImageSize,lImageSize])),gray);
    %lV = handles.X(idx,:)*squeeze(handles.Data.Nets(5).V(handles.Data.Nets(5).InvIdxs(idx),:,:));
    lV = squeeze(handles.Data.Nets(5).V(handles.Data.Nets(5).InvIdxs(idx),:,:));
end;

% Also display the svd directions

%axes(handles.axTan1); subimage(double(reshape(lV(:,1)*(squeeze(handles.Data.Nets(5).V(handles.Data.Nets(5).InvIdxs(idx),:,1)))',[lImageSize,lImageSize])),jet);
axes(handles.axTan1); imagesc(double(reshape(lV(:,1),[lImageSize,lImageSize])));%colormap(gray);
axes(handles.axTan2); imagesc(double(reshape(lV(:,2),[lImageSize,lImageSize])));%colormap(gray);
axes(handles.axTan3); imagesc(double(reshape(lV(:,3),[lImageSize,lImageSize])));%colormap(gray);
axes(handles.axTan4); imagesc(double(reshape(lV(:,4),[lImageSize,lImageSize])));%colormap(gray);
catch
end;

return;


%
%
% Updates the visualization in the MultiscaleSVDAxes
%
%
function DrawMultiscaleSVDAxes( hObject )

handles = guidata( hObject );

if length(handles.CoordAxes_SelPtIdxs)>0,
    idx = handles.CoordAxes_SelPtIdxs(length(handles.CoordAxes_SelPtIdxs));

    for k =  1:length(handles.Axes{1}.Data),
        if isempty(handles.Axes{1}.Data{k}),
            continue;
        end;
        handles.Axes{1}.Data{k}.Yfixparam = idx;
    end;

    guidata( hObject, handles);

    DrawDataAxes( handles.Axes{1} );
end;

return;


%
%
% Updates the visualization in the MultiscaleSVDAxes
%
%
function DrawGUIDataAxes( hObject, cAxesIdx )

handles = guidata( hObject );

if length(handles.CoordAxes_SelPtIdxs)>0,
    idx = handles.CoordAxes_SelPtIdxs(length(handles.CoordAxes_SelPtIdxs));

    for k =  1:length(handles.Axes{cAxesIdx}.Data),
        if isempty(handles.Axes{cAxesIdx}.Data{k}),
            continue;
        end;
        handles.Axes{cAxesIdx}.Data{k}.Yfixparam = idx;
        if isfield(handles,'XLogScale'),
            handles.Axes{cAxesIdx}.Data{k}.XLogScale = handles.XLogScale;
        end;
        if isfield(handles,'YLogScale'),
            handles.Axes{cAxesIdx}.Data{k}.YLogScale = handles.YLogScale;
        end;
    end;

    guidata( hObject, handles);

    DrawDataAxes( handles.Axes{cAxesIdx} );
end;

return;



%
% DrawDataAxes
%
% Rather general function to display plots of data
%
% IN:
%   cData: structure with the following fields:
%           Data : cell array of structures with the following fields:
%                   X   : data for the X axis
%                   Y   : data for the Y axis
%                   Xfixdim,
%                   Xfixparam : will squeeze out the Xfixdim-th dimension of X by setting its value to Xfixparam. 
%                               Right now works only if scalar or empty.
%                   Yfixdim,
%                   Yfixparam : will squeeze out the Yfixdim-th dimension of Y by setting its value to Yfixparam. 
%                               Right now works only if scalar or empty.
%                   PlotStyle : plot style passed to plot command
%                   XLogScale : plot x axis in log scale
%                   YLogScale : plot y axis in log scale
%           Xlbl : label of X axis
%           Ylbl : label of Y axis
%           title: title of axis
%
% USES:
%   EliminateFixDim
%
% Comments: only works if X, after squeezing the Xfixdim-th dimension, is one-dimensional. No error checking.
%           only works if Y, after squeezing the Xfixdim-th dimension, is one- or two-dimensional. No error checking.
%
%

function DrawDataAxes( cData )

axes(cData.handle);cla;
hold on;

% Draw the various data plots
for k = 1:length(cData.Data),
    lX  = EliminateFixDim( cData.Data{k}.X, cData.Data{k}.Xfixdim, cData.Data{k}.Xfixparam );
    lY  = EliminateFixDim( cData.Data{k}.Y, cData.Data{k}.Yfixdim, cData.Data{k}.Yfixparam );
    if isempty(lY), continue; end;
    if isfield(cData.Data{k},'XLogScale'),
        if cData.Data{k}.XLogScale,
            lX = log(lX);
        end;
    end;
    if isfield(cData.Data{k},'YLogScale'),
        if cData.Data{k}.YLogScale,
            lY = log(double(lY));
            lY(find(isinf(lY)))=NaN;
        end;
    end;
    plot( lX, lY, cData.Data{k}.PlotStyle );
    axis tight;
    try
    if isfield(cData,'GoodScales'),
        YLims = get(gca,'YLim');
        if iscell(cData.GoodScales),
            lGoodScales = cData.GoodScales{cData.Data{k}.Yfixparam};
        else
            lGoodScales = cData.GoodScales;
        end;
        p=patch([min(lGoodScales),min(lGoodScales),max(lGoodScales),max(lGoodScales)],[min(YLims),max(YLims),max(YLims),min(YLims)],'g');
        set(p,'FaceAlpha',0.1);
    end;
    catch
    end;
end;
grid on;
% Axes labels and title
xlabel(cData.Xlbl);ylabel(cData.Ylbl);
title(cData.title);
axis tight;

return;




function vX = EliminateFixDim( cX, cFixDim, cFixParam )

vX = [];

if isempty(cFixDim),
    vX = cX;
else
    lNDims = length(size(cX));
    lString = 'cX(';
    for k = 1:lNDims,
        if k~=cFixDim,
            lString = [lString,':'];
        else
            lString = [lString,'cFixParam'];
        end;
        if k<lNDims,
            lString = [lString, ','];
        end;
    end;
    lString = ['vX=squeeze(',lString,'));'];
    try
    eval(lString);
    catch
    end;
end;

return;

%
%
% Updates the visualization in the CoordAxes
%
%
function DrawCoordAxes( hObject )

handles = guidata(hObject);

try
    handles.CoordAxes_NNInfo = [];

    %handles.DisplayCoords = handles.Data.G.EigenVecs(:,[handles.DisplayCoordsIdxs(1),handles.DisplayCoordsIdxs(2),handles.DisplayCoordsIdxs(3)]);

    axes(handles.CoordAxes);hold on;
    lButtonDownFcn=get(gca,'ButtonDownFcn');
    try
        delete(handles.CoordAxes_Pts);
        delete(handles.CoordAxes_Pts2);
    catch
    end;
    if isempty(handles.RestrictIdxs),    
        x = handles.DisplayCoords(:,1);
        y = handles.DisplayCoords(:,2);
        z = handles.DisplayCoords(:,3);
        color = handles.CoordAxes_PtColor;
    else
        x = handles.DisplayCoords(handles.RestrictIdxs,1);
        y = handles.DisplayCoords(handles.RestrictIdxs,2);
        z = handles.DisplayCoords(handles.RestrictIdxs,3);
        color = handles.CoordAxes_PtColor(handles.RestrictIdxs);
    end;
    handles.CoordAxes_Pts = scatter3(x,y,z,20,color,'filled');
    %[XI,YI] = meshgrid([min(x):(max(x)-min(x))/40:max(x)],[min(y):(max(y)-min(y))/40:max(y)]); %ZI = griddata(x,y,z,XI,YI,'cubic'); %handles.CoordAxes_Pts = surf(XI,YI,ZI); %hold %shading interp %view(3);
    axis tight; grid on;
    h=colorbar('SouthOutside');    
    set(gca,'ButtonDownFcn',lButtonDownFcn);
    set(handles.CoordAxes_Pts,'ButtonDownFcn',lButtonDownFcn)
    handles.CoordAxes_Pts2 = scatter3(x,y,min(z)*ones(length(x),1),2,zeros(size(x)));
    set(handles.CoordAxes_Pts2,'ButtonDownFcn',lButtonDownFcn)

    if handles.CoordAxes2D,
        title(sprintf('Diffusion map (\\Phi_{%d},\\Phi_{%d})',handles.DisplayCoordsIdxs(1),handles.DisplayCoordsIdxs(2)));
    else
        title(sprintf('Diffusion map (\\Phi_{%d},\\Phi_{%d},\\Phi_{%d})',handles.DisplayCoordsIdxs(1),handles.DisplayCoordsIdxs(2),handles.DisplayCoordsIdxs(3)));
    end;
%    colormap('hot');
catch
    % Update handles structure
    guidata( gcf, handles);
end;

% Update handles structure
guidata( gcf, handles);


return;


% --- Executes on button press in btnCoordAxesYDown.
function btnCoordAxesYDown_Callback(hObject, eventdata, handles)
% hObject    handle to btnCoordAxesYDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(gcbo);

handles.DisplayCoordsIdxs(2) = max([handles.DisplayCoordsIdxs(2)-1,1]);
handles.DisplayCoords(:,2) = handles.Data.G.EigenVecs(:,handles.DisplayCoordsIdxs(2));

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- Executes on button press in btnCoordAxesYup.
function btnCoordAxesYup_Callback(hObject, eventdata, handles)
% hObject    handle to btnCoordAxesYup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(gcbo);

handles.DisplayCoordsIdxs(2) = min([handles.DisplayCoordsIdxs(2)+1,size(handles.Data.G.EigenVecs,2)]);
handles.DisplayCoords(:,2) = handles.Data.G.EigenVecs(:,handles.DisplayCoordsIdxs(2));

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- Executes on button press in btnCoordAxesXup.
function btnCoordAxesXup_Callback(hObject, eventdata, handles)
% hObject    handle to btnCoordAxesXup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(gcbo);

handles.DisplayCoordsIdxs(1) = min([handles.DisplayCoordsIdxs(1)+1,size(handles.Data.G.EigenVecs,2)]);
handles.DisplayCoords(:,1) = handles.Data.G.EigenVecs(:,handles.DisplayCoordsIdxs(1));

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- Executes on button press in btnCoordAxesXdown.
function btnCoordAxesXdown_Callback(hObject, eventdata, handles)
% hObject    handle to btnCoordAxesXdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(gcbo);

handles.DisplayCoordsIdxs(1) = max([handles.DisplayCoordsIdxs(1)-1,1]);
handles.DisplayCoords(:,1) = handles.Data.G.EigenVecs(:,handles.DisplayCoordsIdxs(1));

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;



% --- Executes on selection change in pmnuColorEigenfunctions.
function pmnuColorEigenfunctions_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuColorEigenfunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');

if handles.CoordAxes2D
    handles.CoordAxes_PtColor = double(handles.Data.G.EigenVecs(:,lSelIdx));
else
    handles.DisplayCoords(:,3)= double(handles.Data.G.EigenVecs(:,lSelIdx));
    handles.DisplayCoordsIdxs(3) = lSelIdx;
end;


guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- Executes during object creation, after setting all properties.
function pmnuColorEigenfunctions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuColorEigenfunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pmnuColorSVD.
function pmnuColorSVD_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuColorSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorSVD;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
lUserData = get(hObject,'UserData');

if get(handles.chkDenoisedSVD,'Value')==0,
    handles.CoordAxes_PtColor =  ...
        double(squeeze(handles.Data.Stats.S(lUserData(lSelIdx,1),:,lUserData(lSelIdx,2))));
else
    handles.CoordAxes_PtColor =  ...
        double(squeeze(handles.Data.Stats.Tychonov.S(lUserData(lSelIdx,1),:,lUserData(lSelIdx,2))));
end;

handles.CoordAxesColorMode = 1;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- Executes during object creation, after setting all properties.
function pmnuColorSVD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuColorSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pmnuColorNSVD.
function pmnuColorNSVD_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuColorNSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorNSVD;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
lUserData = get(hObject,'UserData');

warning off;
if get(handles.chkDenoisedSVD,'Value')==0,
    handles.CoordAxes_PtColor =  double( ...
        squeeze(handles.Data.Stats.S(lUserData(lSelIdx,1),:,lUserData(lSelIdx,2)))./ ...
        squeeze(handles.Data.Stats.S(lUserData(lSelIdx,1),:,1)));
else
    handles.CoordAxes_PtColor =  double( ...
        squeeze(handles.Data.Stats.Tychonov.S(lUserData(lSelIdx,1),:,lUserData(lSelIdx,2)))./ ...
        squeeze(handles.Data.Stats.Tychonov.S(lUserData(lSelIdx,1),:,1)));
end;
warning on;
handles.CoordAxes_PtColor(find(isnan(handles.CoordAxes_PtColor)))=0;

handles.CoordAxesColorMode = 2;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- Executes during object creation, after setting all properties.
function pmnuColorNSVD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuColorNSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btnRotate.
function btnRotate_Callback(hObject, eventdata, handles)
% hObject    handle to btnRotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

rotate3d(handles.CoordAxes);

return;


% --- Executes on button press in btn23D.
function btn23D_Callback(hObject, eventdata, handles)
% hObject    handle to btn23D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

axes(handles.CoordAxes);

lCurview = view;

if lCurview(3,3)==-1,
    view(3);
    handles.CoordAxes2D = false;
    set(handles.txtZaxis,'Visible','on');
else
    view(2);
    handles.CoordAxes2D = true;
    set(handles.txtZaxis,'Visible','off');
end;

guidata( hObject, handles );

DrawCoordAxes( hObject );

return;



% --- Executes on button press in chkDenoisedSVD.
function chkDenoisedSVD_Callback(hObject, eventdata, handles)
% hObject    handle to chkDenoisedSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( gcbo );

if handles.CoordAxesColorMode == 1,
    pmnuColorSVD_Callback;
elseif handles.CoordAxesColorMode == 2,
    pmnuColorNSVD_Callback;
end;

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;



% --- Executes on slider movement.
function sldDenoising_Callback(hObject, eventdata, handles)
% hObject    handle to sldDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

handles.Data.Stats.Opts.Lambda = get(hObject,'Value');
handles.Data.Stats = GetStatisticsMultiscaleSVD( [], handles.Data.Nets, handles.Data.Stats.Opts );

guidata( hObject, handles );

DrawMultiscaleSVDAxes( hObject );

return;

% --- Executes during object creation, after setting all properties.
function sldDenoising_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sldDenoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on selection change in pmnuColorJFr.
function pmnuColorJFr_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuColorJFr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorJFr;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');

handles.CoordAxes_PtColor =  double(squeeze(handles.Data.Stats.JFr(:,lSelIdx)));

handles.CoordAxesColorMode = 2;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;


% --- Executes during object creation, after setting all properties.
function pmnuColorJFr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuColorJFr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btnFreeMemory.
function btnFreeMemory_Callback(hObject, eventdata, handles)
% hObject    handle to btnFreeMemory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

if isfield(handles.Data.G,'T'),
    handles.Data.G.T = [];
end;
if isfield(handles.Data.G,'W'),
    handles.Data.G.W = [];
end;
if isfield(handles.Data.G,'P'),
    handles.Data.G.P = [];
end;
if isfield(handles.Data.Stats,'J'),
    handles.Data.Stats.J = [];
end;
if isfield(handles.Data.Stats,'Volume'),
    handles.Data.Stats.Volume = uint32(handles.Data.Stats.Volume);
end;

guidata( hObject, handles );

return;


% --- Executes on selection change in pmnuColorResVar.
function pmnuColorResVar_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuColorResVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorResVar;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
lUserData = get(hObject,'UserData');
if ~isfield(handles.Data.Stats,'ResVar'),
    handles.Data.Stats = GetStatisticsMultiscaleSVD( [], handles.Data.Nets, struct('NumberOfPoints',size(handles.Data.Coords,1),'NetOpts',handles.Data.NetsOpts,'Lambda',40,'Iter',20) );
end;
handles.CoordAxes_PtColor =  double(squeeze(handles.Data.Stats.ResVar(lUserData(lSelIdx,1),:,lUserData(lSelIdx,2))));

handles.CoordAxesColorMode = 3;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;


% --- Executes during object creation, after setting all properties.
function pmnuColorResVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuColorResVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pmnuColorNResVar.
function pmnuColorNResVar_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuColorNResVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorNResVar;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
lUserData = get(hObject,'UserData');
if ~isfield(handles.Data.Stats,'ResVar'),
    handles.Data.Stats = GetStatisticsMultiscaleSVD( [], handles.Data.Nets, struct('NumberOfPoints',size(handles.Data.Coords,1),'NetOpts',handles.Data.NetsOpts,'Lambda',40,'Iter',20) );
end;
handles.CoordAxes_PtColor =  double(squeeze(handles.Data.Stats.ResVar(lUserData(lSelIdx,1),:,lUserData(lSelIdx,2))));
warning off;
handles.CoordAxes_PtColor = handles.CoordAxes_PtColor./double(squeeze(handles.Data.Stats.ResVar(lUserData(lSelIdx,1),:,1)));
warning on;

handles.CoordAxes_PtColor(find(isnan(handles.CoordAxes_PtColor)))=0;

handles.CoordAxesColorMode = 3;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;



% --- Executes during object creation, after setting all properties.
function pmnuColorNResVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuColorNResVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function mnuFile_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuFileExit_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFileExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf);


% --- Executes on selection change in pmnuSecondAxis.
function pmnuSecondAxis_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuSecondAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
lUserData = get(hObject,'UserData'); lUserData = lUserData(lSelIdx,:);
lStrings = get(hObject,'String');

lVarName = lStrings{lSelIdx};

PrepareAxesData( hObject, 1, lVarName, lUserData );

DrawMultiscaleSVDAxes( hObject );

return;



function PrepareAxesData( hObject, cAxisIdx, cVarName, cUserData )

handles = guidata( hObject );

handles.Axes{cAxisIdx}.Data{1}.X = handles.Data.NetsOpts.Delta(1:handles.Data.NetsOpts.NumberOfScales);
handles.Axes{cAxisIdx}.Data{1}.Xfixdim = [];
handles.Axes{cAxisIdx}.Data{1}.Xfixparam = [];
eval(['handles.Axes{cAxisIdx}.Data{1}.Y = handles.',cVarName,';']);
handles.Axes{cAxisIdx}.Data{1}.Yfixdim = cUserData(1);
handles.Axes{cAxisIdx}.Data{1}.Yfixparam = cUserData(2);
handles.Axes{cAxisIdx}.Data{1}.PlotStyle = 'x-';

cVarName2 = strrep(cVarName,'Data.Stats.','Data.Stats.Tychonov.');

try
    handles.Axes{cAxisIdx}.Data{2}.X = handles.Data.NetsOpts.Delta(1:handles.Data.NetsOpts.NumberOfScales);
    handles.Axes{cAxisIdx}.Data{2}.Xfixdim = [];
    handles.Axes{cAxisIdx}.Data{2}.Xfixparam = [];
    eval(['handles.Axes{cAxisIdx}.Data{2}.Y = handles.',cVarName2,';']);
    handles.Axes{cAxisIdx}.Data{2}.Yfixdim = cUserData(1);
    handles.Axes{cAxisIdx}.Data{2}.Yfixparam = cUserData(2);
    handles.Axes{cAxisIdx}.Data{2}.PlotStyle = 'x-';
    handles.Axes{cAxisIdx}.Data{1}.PlotStyle = 'o:';
catch
    handles.Axes{cAxisIdx}.Data = cell({handles.Axes{cAxisIdx}.Data{1}});
end;

handles.Axes{cAxisIdx}.Xlbl = 'Scale';
handles.Axes{cAxisIdx}.Ylbl = cVarName;
handles.Axes{cAxisIdx}.title = cVarName;

guidata( hObject, handles );

return;

% --- Executes during object creation, after setting all properties.
function pmnuSecondAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuSecondAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

return;


% --- Executes on selection change in pmnuColorNetInvIdxs.
function pmnuColorNetInvIdxs_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuColorNetInvIdxs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorNResVar;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
handles.CoordAxes_PtColor =  double(handles.Data.Nets(lSelIdx).InvIdxs);

handles.CoordAxesColorMode = 3;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;





% --- Executes during object creation, after setting all properties.
function pmnuColorNetInvIdxs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuColorNetInvIdxs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pmnuThirdAxes.
function pmnuThirdAxes_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuThirdAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
lUserData = get(hObject,'UserData'); 
lUserData = lUserData(lSelIdx,:);
lStrings = get(hObject,'String');

lVarName = lStrings{lSelIdx};

PrepareAxesData( hObject, 2, lVarName, lUserData );

DrawGUIDataAxes( hObject, 2 );

return;

% --- Executes during object creation, after setting all properties.
function pmnuThirdAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuThirdAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pmnuVolume.
function pmnuVolume_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorNResVar;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
handles.CoordAxes_PtColor =  double(handles.Data.Stats.Volume(lSelIdx,:));

handles.CoordAxesColorMode = 3;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;





% --- Executes during object creation, after setting all properties.
function pmnuVolume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in pmnuRestrictToIndices.
function pmnuRestrictToIndices_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuRestrictToIndices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

lRandPerm = randperm(handles.NumberOfDataPoints);
lString = get(hObject,'String');
lValue = get(hObject,'Value');

if lValue>1,
    handles.RestrictIdxs = lRandPerm(1:str2num(lString{lValue}));
else
    handles.RestrictIdxs = [];
end;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- Executes during object creation, after setting all properties.
function pmnuRestrictToIndices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuRestrictToIndices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btnLogScale.
function btnLogScale_Callback(hObject, eventdata, handles)
% hObject    handle to btnLogScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'XLogScale'),
    handles.XLogScale = false;
end;

handles.XLogScale = ~handles.XLogScale;

guidata( hObject, handles );

DrawGUIDataAxes( hObject, 1 );
DrawGUIDataAxes( hObject, 2 );

return;



% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'YLogScale'),
    handles.YLogScale = false;
end;

handles.YLogScale = ~handles.YLogScale;

guidata( hObject, handles );

DrawGUIDataAxes( hObject, 1 );
DrawGUIDataAxes( hObject, 2 );

return;





% --- Executes on selection change in pmnuLblCoordAxes.
function pmnuLblCoordAxes_Callback(hObject, eventdata, handles)
% hObject    handle to pmnuLblCoordAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if nargin<1,
    hObject = gcf;
    handles = guidata( gcf );
    hObject = handles.pmnuColorNResVar;
end;

handles = guidata( hObject );

lSelIdx = get(hObject,'Value');
lUserData = get(hObject,'UserData');
handles.CoordAxes_PtColor =  double(squeeze(handles.Data.Labels(:,lUserData(lSelIdx,1))));
handles.CoordAxes_PtColor(find(isnan(handles.CoordAxes_PtColor)))=0;

handles.CoordAxesColorMode = 3;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;




% --- Executes during object creation, after setting all properties.
function pmnuLblCoordAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnuLblCoordAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata( hObject );

lStatus = get(hObject,'Value');

if lStatus,
    handles.DisplayCoords = handles.Data.XCoords;
else
    handles.DisplayCoords = handles.Data.Coords;
end;

guidata( hObject, handles );

DrawCoordAxes( hObject );
DrawCoordAxesSelPt( hObject );

return;

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over checkbox2.
function checkbox2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


