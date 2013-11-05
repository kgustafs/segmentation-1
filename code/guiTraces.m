function varargout = guiTraces(varargin)
% GUITRACES MATLAB code for guiTraces.fig
%      GUITRACES, by itself, creates a new GUITRACES or raises the existing
%      singleton*.
%
%      H = GUITRACES returns the handle to a new GUITRACES or the handle to
%      the existing singleton*.
%
%      GUITRACES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUITRACES.M with the given input arguments.
%
%      GUITRACES('Property','Value',...) creates a new GUITRACES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiTraces_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiTraces_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiTraces

% Last Modified by GUIDE v2.5 28-Feb-2013 16:22:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiTraces_OpeningFcn, ...
                   'gui_OutputFcn',  @guiTraces_OutputFcn, ...
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


% --- Executes just before guiTraces is made visible.
function guiTraces_OpeningFcn(hObject, eventdata, handles, varargin)

    % Choose default command line output for guiTraces
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    
    %set(gcf,'doublebuffer','off');
    set(gcf, 'WindowButtonMotionFcn', @mouseMove);
    set(gcf, 'WindowKeyPressFcn', @KeyPress)

    global k t;
    k=1;
    t=1;

    global moviePath Ntraces Nframes; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% MAIN PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    moviePath = [pwd() '/'] %'/Users/bieler/Desktop/matlab/13sept13_longdaysin1to5uM/movie35/';
 
    
    
    %moviePath =  get(handles.edit1,'string');    
        
    set(handles.edit1,'string',moviePath);    
         
    global axesData axesGif;
    axesData = handles.axes1;
    axesGif = handles.axes2;

    global tracks signal traj ind longTraces peakMatrix divMatrix divisions modK speed;
    
    modK=0;
    speed = 1;
      
    % load stuff
    load tracks.mat tracks
    load signal.mat signal
    load traj.mat traj
    load ind.mat ind
    
    Nframes = size(signal,2);
    
    load divisions.mat divisions

    if(exist('longTraces.mat','file'))
       load longTraces.mat 
    end
    
    if(exist('hadPeaksChanged.mat','file'))
       load hadPeaksChanged.mat       
    else
       load peakMatrix.mat
    end
    
    if(exist('divMatrixFinal.mat','file'))
       load divMatrixFinal.mat       
    else
       load divMatrix.mat
    end

    Ntraces = length(longTraces);     
    
    global gif gmap;    
    %[gif gmap] = imread([moviePath 'snapShots/' num2str(longTraces(k)) '.gif'],'frames','all');
    
    set(handles.text3,'String',num2str(longTraces(k)));
    updatePlot();
    updatePlotGif();

function  updatePlotGif()
    
    global gif t axesGif gmap longTraces k;      

    %expe.t = linspace(0,72,168);
    
    axes(axesGif);
    %imagesc(medfilt2( double( gif(:,:,:,t) ), 3*[1 1]))
    
    if(exist(['snapShots/' num2str(longTraces(k)) '_' num2str(t) '.png'],'file'))
        im = imread(['snapShots/' num2str(longTraces(k)) '_' num2str(t) '.png']);
        imagesc(im);
        caxis([min(im(:)) 50])
    end
    
    %imagesc(gif(:,:,1,t))
    %caxis([20, 230])
    colormap hot
    
function KeyPress(Source, EventData)
    
    handles = guidata(gcf);
    global t speed gmap axesData moviePath k tracks signal traj ind longTraces peakMatrix divMatrix Nframes divisions Ntraces gif;
    
    hasKChanged = 0;
    hasTChanged = 0;
    
    hadPeaksChanged = 0;
    hadDivsChanged = 0;


    switch EventData.Key
        
        case 'shift'
            
            if(speed == 1)
                speed = 5;
            else
                speed = 1;
            end
        
        case 'leftarrow'
            
            t = max(t-speed,1);
            updatePlotGif();
            
        case 'rightarrow'
            
            t = min(t+speed,Nframes);
            updatePlotGif();
        
    end
    
    switch EventData.Character
        
        %remove everything
        case 'k'
            divMatrix(longTraces(k),:) = 0;
            hadDivsChanged=1;
            peakMatrix(longTraces(k),:) = 0;
            
            hadPeaksChanged=1;
        
        %remove division        
        case 'a'
            
            sel = max(t-5,1):min(t+5,Nframes);            
            divMatrix(longTraces(k),sel) = 0;
            
            hadDivsChanged=1;
            
        %add division        
        case 'd'
            
            sel = max(t-5,1):min(t+5,Nframes);            
            divMatrix(longTraces(k),sel) = 0;
            divMatrix(longTraces(k),t) = 1;
            
            hadDivsChanged=1;
            
        %remove peak        
        case 'q'
            
            sel = max(t-5,1):min(t+5,Nframes);            
            peakMatrix(longTraces(k),sel) = 0;
            
            hadPeaksChanged=1;
            
        %add peak        
        case 'e'
            
            sel = max(t-5,1):min(t+5,Nframes);            
            peakMatrix(longTraces(k),sel) = 0;
            peakMatrix(longTraces(k),t) = 1;
            
            hadPeaksChanged=1;
                    
        case 't'
            
            C = get(axesData, 'CurrentPoint');
            title(axesData, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);

            mouseX = C(1,1);
            mouseY = C(1,2);
        
            frame = round( Nframes*mouseX/72 );
            frame = min(frame, Nframes);
            frame = max(frame, 1);
            
            t = frame;
            hasTChanged = 1;

        case '1'

            k = max(k-1,1);
            hasKChanged=1;
            
            
        case '2'

            k=min(k+1,Ntraces);
            hasKChanged=1;
                        
    end
    
    updatePlot();
    
    if(hadDivsChanged)
        save divMatrixFinal.mat divMatrix
    end
    
    if(hadPeaksChanged)
        save hadPeaksChanged.mat peakMatrix
    end
    
    if(hasKChanged)
                
        set(handles.text3,'String',[num2str(longTraces(k)) ' (' num2str(k) '/' num2str(length(longTraces)) ')'] );
        
       % [gif gmap] = imread([moviePath 'snapShots/' num2str(longTraces(k)) '.gif'],'frames','all');
        updatePlotGif();
    end
    if(hasTChanged)
        updatePlotGif();
    end
    
    
function mouseMove(hObject,eventdata)
   
    global axesData t Nframes modK;
    
    modK = modK+1;
    
    %if( mod(modK,3)==0 )
    if(1) 
    
        modK = 0;

        C = get(axesData, 'CurrentPoint');
        %title(axesData, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);

        mouseX = C(1,1);
        mouseY = C(1,2);

        frame = round( Nframes*mouseX/72 );
        frame = min(frame, Nframes);
        frame = max(frame, 1);

        t = frame;

        updatePlotGif();

    end
    
    
    
    
%     %load data
%     handles = guidata(gcf);
%     
%     seg = handles.seg;
%     axesSeg = handles.axesSeg;
%     dKey = handles.dKey;
%     previousMouseX = handles.previousMouseX;
%     previousMouseY = handles.previousMouseY;
%     
%     %%change cursor
%     if(dKey == 1)
%         set(gcf,'Pointer','cross');
%     elseif(aKey==1)
%         set(gcf,'Pointer','crosshair');
%     else
%         set(gcf,'Pointer','arrow');
%     end
        
%     C = get(axesData, 'CurrentPoint');
%     %title(axesSeg, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
% 
%     mouseX = C(1,1);
%     mouseY = C(1,2);
%     

    
function updatePlot()
    
    global signal t Ntraces longTraces k axesData peakMatrix divMatrix ind divisions Nframes;     

    expe.t = linspace(0,72,Nframes);
    
    axes(axesData);
    cla;
    hold on;
            
    idx = longTraces(k);
    tmp = signal(idx,:);
         
    maxp = find( peakMatrix(idx,:))';
         
    nonZeroIdx = find(ind(idx,:) > 0);
    nonZeroIdx = nonZeroIdx(1);
    
    divInd = find( idx == [divisions.sisterInd] );
    
    if(length(divInd) == 1)
       
        plot(expe.t(1:divisions(divInd).motherFrame), signal(divisions(divInd).motherInd,1:divisions(divInd).motherFrame),'color','k','lineWidth',1)
                
    end
    
    for j=1:size(maxp,1)
        plot([expe.t(maxp(j,1)) expe.t(maxp(j,1))],[0 300],'color',0.6*[1 1 1])
    end
    per=[];
    for j=1:size(maxp,1)-1
        per(j) = expe.t(maxp(j+1,1)) - expe.t(maxp(j,1));

        offset = 5+5*rand;
        plot([expe.t(maxp(j,1)) expe.t(maxp(j,1)) + per(j)],tmp(maxp(j,1))*[1 1],'k')
        text(expe.t(maxp(j,1)) + per(j)/2,tmp(maxp(j,1))*1-offset-0.05,num2str(per(j)))
    end
            
    plot(expe.t(1:length(tmp)), tmp,'color','k','lineWidth',2)
    plot(expe.t(t), tmp(t),'ks')
    
    % plot division
    dd = find(divMatrix(idx,:) > 0 );
    for j=1:length(dd)
       
        offset = 5*randn;           
           
        ff = dd(j);
        plot([expe.t(ff) expe.t(ff)],[0 300],'color',[1 0.6 0.6])
               
        if( ~isempty(maxp) )
        
            ppeaks = maxp(:,1);
            ppeaks = ppeaks(ppeaks < ff);            
            divTiming = find(abs( ff - ppeaks) == min(abs(ff - ppeaks)) );
            %divTiming = divTiming(1);
            
            if(isempty(divTiming))
                ppeaks = maxp(:,1);
                ppeaks = ppeaks(ppeaks < ff);            
                divTiming = find(abs( ff - ppeaks) == min(abs(ff - ppeaks)) );                
            end
                        
            if(~isempty(divTiming))
                divTiming = divTiming(1);

                plot([expe.t(ff) expe.t(maxp(divTiming,1))],10+max(tmp)*[1 1]+offset,'r')

                divTimingHours = expe.t(ff) - expe.t(maxp(divTiming,1));

                text(mean([expe.t(ff) expe.t(maxp(divTiming,1))]),[max(tmp)+15+offset],num2str(divTimingHours))
            end
                        
        end
    end
        
    xlabel('time [h]')
    ylabel('mean pixel intensity')
    
%    set(gca,'XTick',  0:2:expe.t(NToTrack) )
    
    axis([0 expe.t(length(tmp)) 0.9*min(tmp) 1.2*(max(tmp))])
    

% --- Outputs from this function are returned to the command line.
function varargout = guiTraces_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
