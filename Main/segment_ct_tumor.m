function varargout = segment_ct_tumor(varargin)
% SEGMENT_CT_TUMOR MATLAB code for segment_ct_tumor.fig
%      SEGMENT_CT_TUMOR, by itself, creates a new SEGMENT_CT_TUMOR or raises the existing
%      singleton*.
%
%      H = SEGMENT_CT_TUMOR returns the handle to a new SEGMENT_CT_TUMOR or the handle to
%      the existing singleton*.
%
%      SEGMENT_CT_TUMOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENT_CT_TUMOR.M with the given input arguments.
%
%      SEGMENT_CT_TUMOR('Property','Value',...) creates a new SEGMENT_CT_TUMOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segment_ct_tumor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segment_ct_tumor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segment_ct_tumor

% Last Modified by GUIDE v2.5 27-Jun-2013 15:15:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segment_ct_tumor_OpeningFcn, ...
                   'gui_OutputFcn',  @segment_ct_tumor_OutputFcn, ...
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


% --- Executes just before segment_ct_tumor is made visible.
function segment_ct_tumor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segment_ct_tumor (see VARARGIN)

    % Choose default command line output for segment_ct_tumor
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    setappdata(0  , 'hMainGui'    , gcf);
    clc
    
    % UIWAIT makes segment_ct_tumor wait for user response (see UIRESUME)
    % uiwait(handles.tumorQ);
    
    % enforce use of native open dialog boxes (nicer than Matlab defaults)
    setappdata(0,'UseNativeSystemDialogs',0);  % http://www.mathworks.com/matlabcentral/newsreader/view_thread/73094
    
    % load a blank default image
    app.blank = imread('setuvo.jpg');
    axes(handles.axes1);
    imagesc(app.blank);axis off;axis equal;    
 
    clearGlobal;
    
    
    
    
% --- Outputs from this function are returned to the command line.
function varargout = segment_ct_tumor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Clear everything from global variable 'app'
function clearGlobal

    global app;
    
    % set up common global variable 'app' fields
    app.roiPath = [];
    app.imgPath = [];
    app.imgPrefix = [];
    app.alpha = 0;
    app.beta = 0;
    app.eta = 0;
    app.T = 0;
    app.maxits = 0;
    app.outputPath = [];
    app.tumor.startZ = 0;
    app.tumor.endZ = 0;
    app.tumor.rect = [];
    app.tumor.initX = 0;
    app.tumor.initY = 0;
    app.tumor.initZ = 0;
    app.vol.image = [];
    app.vol.mouse = [];
    app.vol.tumor = [];
    app.vol.mini = [];
    app.vol.mseg = [];
    app.vol.seg = [];
    app.vol.phi = [];
    app.vol.tmap = [];
    app.vol.init = [];
    app.vol.lp = [];
    app.currentZ = 0;
    app.stage = 0;      % 0 = whole image, 1 = tumor image 2 = seg done
    app.dx = 0;
    app.dy = 0;
    app.dz = 0;
    app.ls_vols = [];
    app.fulldx = 0;
    app.fulldy = 0;
    app.fulldz = 0;
    app.minidx = 0;
    app.minidy = 0;
    app.minidz = 0;
    app.imgType = 0;    % 0 = no image; 1 = IMG; 2 = DCM (from RT_image)
    app.downsize = 0; % 0 = run at native resolution, 1 = downsample x2
    
    
% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
% hObject    handle to loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global app
    
    set(handles.status,'String','Loading CT image data...');
    pause(0.1);
    
    app.imgPath = get(handles.imgPath,'String');   
        
    if app.imgPath(end-2:end) == 'img'
        app.fulldx  = str2num(get(handles.dx,'String'));
        app.fulldy  = str2num(get(handles.dy,'String'));
        app.fulldz  = str2num(get(handles.dz,'String'));
        tic
        fh = fopen(app.imgPath,'r','l');
        bin = fread(fh,'int16');
        fclose(fh);
        vol = reshape(bin,[app.fulldx app.fulldy app.fulldz]);
        clear bin
        tLoad = toc;
        app.imgType = 1;
    elseif app.imgPath(end-2:end) == 'dcm'
        tic;
        vol = dicomread(app.imgPath);
        vol = double(squeeze(vol));
        vol = flipdim(vol,2); % do a LR flip
        [app.fulldx app.fulldy app.fulldz] = size(vol);
        tLoad = toc;
        app.imgType = 2;
    end
    
    set(handles.dx,'String',num2str(app.fulldx));
    set(handles.dy,'String',num2str(app.fulldy));
    set(handles.dz,'String',num2str(app.fulldz));    
    
    app.vol.image = vol;
    [app.dx app.dy app.dz] = size(app.vol.image);
    app.currentZ = round(app.fulldz/2);
    app.stage = 0;
    updateView;
        
    set(handles.status,'String',['Loaded CT dataset in ' num2str(tLoad) ' seconds']);
    set(handles.setZ,'Value',0.5);
    
    set(handles.rectROI,'Enable','ON');
    set(handles.startZ,'Enable','ON');
    set(handles.endZ,'Enable','ON');
    set(handles.initPoint,'Enable','ON');
    

function updateView

    global app
    
    hMainGui = getappdata(0, 'hMainGui');
    handles  = guidata(hMainGui); 

    set(handles.axes1,'HandleVisibility','ON');
    axes(handles.axes1);
    if app.stage == 0           % show whole image
        if app.imgType == 1
            imagesc(rotateCT(app.vol.image(:,:,app.currentZ)));
        elseif app.imgType == 2
            imagesc(app.vol.image(:,:,app.currentZ));
        end
    elseif app.stage == 1       % show tumor seg volume
        if app.imgType == 1
            imagesc(rotateCT(app.vol.tumor(:,:,app.currentZ)));
        elseif app.imgType == 2
            imagesc(app.vol.tumor(:,:,app.currentZ));
        end
    elseif app.stage == 2       % show LS seg vs initialisation
        if app.imgType == 1
            fullimg = superimpose_binary_map(rotateCT(app.vol.tumor(:,:,app.currentZ)), ...
                rotateCT(app.vol.seg(:,:,app.currentZ)),rotateCT(app.vol.init(:,:,app.currentZ)));
        elseif app.imgType == 2
            fullimg = superimpose_binary_map(app.vol.tumor(:,:,app.currentZ), ...
                app.vol.seg(:,:,app.currentZ),app.vol.init(:,:,app.currentZ));
        end
        imagesc(fullimg);
    elseif app.stage == 3       % show LS seg vs manual seg
        if app.imgType == 1
            fullimg = superimpose_binary_map(rotateCT(app.vol.tumor(:,:,app.currentZ)), ...
                rotateCT(app.vol.seg(:,:,app.currentZ)),rotateCT(app.vol.mseg(:,:,app.currentZ)));
        elseif app.imgType == 2
            fullimg = superimpose_binary_map(app.vol.tumor(:,:,app.currentZ), ...
                app.vol.seg(:,:,app.currentZ),app.vol.mseg(:,:,app.currentZ));            
        end
        imagesc(fullimg);
    end
    % set(handles.axes1,'CLim',[cmin cmax]);
    colormap('gray');
    axis image, axis off; axis tight;  
    set(handles.axes1,'HandleVisibility','OFF');   


function preprocess_data
        
    global app
    
    hMainGui = getappdata(0, 'hMainGui');
    handles  = guidata(hMainGui); 
    
    disp('Calculating local phase...');  
    
    [A le lp] = monogenic_3D(app.vol.tumor,1,3.25,0.25);
    lp = norm_volume(lp);
    disp('Done with calculating local phase...');
            
    disp('Setting up exclusion zones');
    vol_exclude = zeros(size(lp));
    vol_bones   = zeros(size(lp));
    vol_bones(app.vol.mini > 0.5) = 1;
    vol_bones = dilateBinaryVolume(vol_bones);
    vol_bones = dilateBinaryVolume(vol_bones);

    vol_exclude(vol_bones == 1) = 1;
    vol_exclude(app.vol.tumor < 0.1) = 1;

    clear vol_bones
    
    % remove the bed
    % vol_exclude(app.vol.tumor < 0.25) = 1;
    
    % test out effect of gradient magnitude isolation of tissue-air voxels
    disp('Calculating gradient magnitude...');   
    [gx gy gz] = gradient(app.vol.tumor);
    gm = sqrt(gx.^2 + gy.^2 + gz.^2);
    
    
    gm1 = gm;
    gm1(gm < 0.04) = 0;
    lp = lp + (app.beta*gm1);
    
%     h = fspecial('gaussian',10,10);
%     gm1 = imfilter(gm,h);
%     gm1(gm < 0.02) = 0;
%     lp = lp + (12*gm1);
    
    app.vol.gm = gm;
    
    % isolate all air / bone regions from speed term
    
    I = lp;
    I(vol_exclude == 1) = 0;
    m = zeros(size(lp));
    
    s = 5;  
    mx = app.tumor.initX;
    my = app.tumor.initY;
    mz = app.tumor.initZ;
    
    m(my-s:my+s,mx-s:mx+s,mz-s:mz+s) = 1;
       
    app.vol.lp = I;
    app.vol.init = m;
    
    
function postprocess_data

    global app
    
    hMainGui = getappdata(0, 'hMainGui');
    handles  = guidata(hMainGui); 
    
    strResultsLabel = app.outputPath;
    disp(['Final volume at original resolution is ' num2str(sum(sum(sum(app.vol.seg)))) ' voxels']);
    
    set(handles.status,'String',['Generating plots and binaries...']);
    % do plot of seg volume vs iterations
    figure;
    plot(app.ls_vols);
    xlabel('Iterations');ylabel('Volume');title('Size of Zero Level Set over Time');
    saveas(gcf,['results/' strResultsLabel '_convergence_plot.png'],'png');

    % create movie flythrough of 3D segmentation
    disp('Creating movie of results');
    figure;
    vidObj = VideoWriter(['results/' strResultsLabel '_3Dflythru.avi']);
    open(vidObj);
    for i = 1 : app.dz
        if app.imgType == 1
            fullimg = superimpose_binary_map(rotateCT(app.vol.tumor(:,:,i)),rotateCT(app.vol.seg(:,:,i)),rotateCT(app.vol.init(:,:,i)));
        elseif app.imgType == 2
            fullimg = superimpose_binary_map(app.vol.tumor(:,:,i),app.vol.seg(:,:,i),app.vol.init(:,:,i));
        end
        imagesc(fullimg);axis image
        M(i) = getframe;
        writeVideo(vidObj,M(i));  
    end
    close(vidObj);    
       
    % save segmentation result as binary file
    disp('Writing results to hard drive as RAW file');    
    seg_final = uint8(zeros(app.fulldx,app.fulldy,app.fulldz));
    seg_final = rect_volume_insert(seg_final, app.vol.seg, app.tumor.rect, app.tumor.startZ, app.tumor.endZ);
    fh = fopen(['results/' strResultsLabel '_fullsize_seg.bin'],'w');
    fwrite(fh,seg_final,'uint8');
    fclose(fh);
    disp('Finished writing results to hard drive as RAW file');
    clear seg_final 

    % write out the tmap (output seg map where each value represents LS iteration number)
    tmap_final = uint8(zeros(app.fulldx,app.fulldy,app.fulldz));
    tmap_final = rect_volume_insert(tmap_final, app.vol.tmap, app.tumor.rect, app.tumor.startZ, app.tumor.endZ);  
    fh = fopen(['results/' strResultsLabel '_tmap.bin'],'w');
    fwrite(fh,tmap_final,'uint8');
    fclose(fh);
    disp('Finished writing tmap results to hard drive as RAW file');
    clear tmap_final
    
    set(handles.status,'String',['Results and plots saved']);
            
    
% --- Executes on slider movement.
function setZ_Callback(hObject, eventdata, handles)
% hObject    handle to setZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    global app
    z = round(get(hObject,'Value')*app.dz);
    if z > app.dz
        z = app.dz;
    end
    if z < 1
        z = 1;
    end
    app.currentZ = z;
    set(handles.status,'String',['Current z = ' num2str(z)]);    
    updateView;
    

% --- Executes during object creation, after setting all properties.
function setZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in startZ.
function startZ_Callback(hObject, eventdata, handles)
% hObject    handle to startZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    global app
    if app.stage == 0
        app.tumor.startZ = app.currentZ;
        set(handles.status,'String',['Tumor Start Z set to ' num2str(app.tumor.startZ)]);
    end

    
% --- Executes on button press in endZ.
function endZ_Callback(hObject, eventdata, handles)
% hObject    handle to endZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global app
    if app.stage == 0
        app.tumor.endZ = app.currentZ;
        set(handles.status,'String',['Tumor End Z set to ' num2str(app.tumor.endZ)]);
    end

    
% --- Executes on button press in rectROI.
function rectROI_Callback(hObject, eventdata, handles)
% hObject    handle to rectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global app
    if app.stage == 0
        set(handles.status,'String','Select tumor with rectangle ROI');
        
        set(handles.axes1,'HandleVisibility','ON');hold on
        axes(handles.axes1);
        rect_tumor = getrect;
        hold off;        
        set(handles.axes1,'HandleVisibility','OFF');
        
        if app.imgType == 1
            app.tumor.rect = convert_CT_rect(floor(rect_tumor),app.dx);
        elseif app.imgType == 2
            app.tumor.rect = floor(rect_tumor);
        end
        
        % check if we need to flip around the startZ and endZ values
        if app.tumor.startZ > app.tumor.endZ
            startZ = app.tumor.endZ;
            endZ   = app.tumor.startZ;
            app.tumor.startZ = startZ;
            app.tumor.endZ   = endZ;
        end
            
        % crop CT image to show tumor ROI
        app.vol.tumor = rect_volume_select(app.vol.image,app.tumor.rect,app.tumor.startZ,app.tumor.endZ);
        app.vol.image = [];
        app.vol.tumor = norm_volume(app.vol.tumor);
        [app.dx app.dy app.dz] = size(app.vol.tumor);
        [app.minidx app.minidy app.minidz] = size(app.vol.tumor);        
        app.currentZ = round(app.dz/2);        
        app.stage = 1;
        
        set(handles.rectROI,'Enable','OFF');
        set(handles.startZ,'Enable','OFF');
        set(handles.endZ,'Enable','OFF');
        set(handles.initPoint,'Enable','ON');
        set(handles.dx,'String',num2str(app.minidx));
        set(handles.dy,'String',num2str(app.minidy));
        set(handles.dz,'String',num2str(app.minidz));   
    end
    
    updateView;
    set(handles.status,'String','Image cropped.');
    set(handles.setZ,'Value',0.5);


% --- Executes on button press in initPoint.
function initPoint_Callback(hObject, eventdata, handles)
% hObject    handle to initPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    global app
    set(handles.status,'String','Select initialisation point inside tumor');
    set(handles.axes1,'HandleVisibility','ON');hold on
    axes(handles.axes1);
    [xpts ypts] = getpts;
    hold off;        
    set(handles.axes1,'HandleVisibility','OFF');
        
    if app.imgType == 1
        app.tumor.initX = floor(ypts(1));
        app.tumor.initY = app.dx-floor(xpts(1));    
    elseif app.imgType == 2
        app.tumor.initX = floor(xpts(1));
        app.tumor.initY = floor(ypts(1));
    end
    app.tumor.initZ = app.currentZ;
        
    set(handles.initPoint,'Enable','OFF');
    set(handles.segment,'Enable','ON');
    set(handles.status,'String','Ready to begin segmentation');
    
    

% --- Executes on button press in segment.
function segment_Callback(hObject, eventdata, handles)
% hObject    handle to segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global app
    
    app.imgPrefix = get(handles.savedPath,'String');    
    
    maxits = str2num(get(handles.its,'String'));
    alpha  = str2num(get(handles.alpha,'String'));
    beta   = str2num(get(handles.beta,'String'));
    eta    = str2num(get(handles.eta,'String'));
    T      = str2num(get(handles.T,'String'));
    label  = get(handles.outputPath,'String');
    app.maxits = maxits;
    app.alpha  = alpha;
    app.beta   = beta;
    app.eta    = eta;
    app.T      = T;
    app.outputPath  = label;           
            
    %if app.stage == 1
        set(handles.status,'String','Pre-processing for segmentation...');
        preprocess_data;
    %elseif app.stage >= 2
    %        app.vol.seg = [];
    %        app.vol.phi = [];
    %        app.vol.tmap = [];
    %        app.ls_vols = [];
    %end
        
    seg  = uint8(zeros(size(app.vol.lp))); % saves seg results in viewable orientation    
    set(handles.status,'String',['Running level set for ' num2str(app.maxits) ' iterations']);
    tic;
    %[seg,phi,ls_vols,tmap] = levelset3DC(app.vol.lp, app.vol.init, app.maxits, app.eta, app.T, app.alpha);
    [seg,phi,ls_vols,tmap] = levelset3DC(app.vol.lp, app.vol.init, app.maxits, app.eta, app.T, app.alpha, 10);
    tseg = toc;
    disp(['Level set completed in ' num2str(tseg) ' secs']);
    disp('=======================================');
    disp('=======================================');
    app.vol.seg  = seg;
    app.vol.phi  = phi;
    app.vol.tmap = tmap;        
    
    app.ls_vols = ls_vols;
    app.stage = 2;
    
    postprocess_data;
    
    set(handles.status,'String',['Level set output saved in results folder']);
    save(['results/' app.outputPath '.mat'],'app');
    set(handles.savedPath,'String',app.outputPath);
    
    updateView;
    


% --- Executes on button press in validate.
function validate_Callback(hObject, eventdata, handles)
% hObject    handle to validate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    global app
    
    set(handles.status,'String',['Loading results data for analysis... may take a few minutes...']);
    pause(0.1);
    
    strResultsLabel = get(handles.savedPath,'String');
    strManualSegLabel = get(handles.roiPath,'String');
    
    % load the saved segmentation results
    load(['results/' strResultsLabel '.mat']);

    % load the manual segmentation ROI file
    mseg = load_roi_file(strManualSegLabel,app.tumor.rect,app.tumor.startZ,app.tumor.endZ,app.fulldx,app.fulldy,app.fulldz);         
    
%     tumor_idx = mseg(app.tumor.initX,app.tumor.initY,app.tumor.initZ);
%     mseg(mseg ~= tumor_idx) = 0;
%     mseg(mseg == tumor_idx) = 1;
    
    % clear out all regions from mask except for largest one
    r = bwconncomp(mseg);
    [biggest,idx] = max(cellfun(@numel,r.PixelIdxList));
    for i = 1 : length(r.PixelIdxList)
        if i ~= idx
            mseg(r.PixelIdxList{i}) = 0;
        end
    end
        
    app.vol.mseg = mseg;
    save(['results/' strResultsLabel '.mat'],'app');
    
    app.stage = 3;
    
    set(handles.status,'String',['Generating plots and binaries...']);
        
    % compute how DICE overlap varies as level set evolves
    maxval = max(max(max(app.vol.tmap)));
    D = zeros(maxval,1);
    for i = 1 : maxval
       tempseg = uint8(app.vol.tmap);
       tempseg1 = uint8(zeros(size(tempseg)));
       % tmap is noisy, have to mask values < 7
       tempseg1(tempseg > 6) = 1;
       tempseg1(tempseg > i) = 0;
       tempseg1(tempseg <= 6) = 0;
       % tempseg1(tempseg == 1) = 1;
       D(i) = dice(tempseg1, mseg);
       clear tempseg tempseg1      
    end
    
    maxD = max(D);
    its = (1:maxval)*10;
    optIt = its(find(D == maxD));
    disp(['Maximum DICE coefficient = ' num2str(maxD) ' at ' num2str(optIt)]);
    disp(['Final DICE coefficient = ' num2str(D(end)) ' at ' num2str(its(end))]);
    figure;
    plot([1:1:maxval]*10,D);
    xlabel('Iterations');
    ylabel('DICE coefficient');
    title('Variation in segmentation agreement with manual boundary vs level set iterations');
    saveas(gcf,['results/' strResultsLabel '_dice_vs_t.png'],'png');    
    
    % create movie flythrough of 3D segmentation
    disp('Creating movie of results');
    set(handles.status,'String','Select initialisation point inside tumor');
    set(handles.axes1,'HandleVisibility','ON');hold on
    axes(handles.axes1);
    vidObj = VideoWriter(['results/' strResultsLabel '_3Dflythru_vs_manual.avi']);
    open(vidObj);
    for i = 1 : app.dz
        if app.imgType == 1
            fullimg = superimpose_binary_map(rotateCT(app.vol.tumor(:,:,i)),rotateCT(app.vol.seg(:,:,i)),rotateCT(mseg(:,:,i)));
        elseif app.imgType == 2
            fullimg = superimpose_binary_map(app.vol.tumor(:,:,i),app.vol.seg(:,:,i),mseg(:,:,i));
        end
        imagesc(fullimg);axis image
        M(i) = getframe;
        writeVideo(vidObj,M(i));  
    end
    close(vidObj);
    hold off;        
    set(handles.axes1,'HandleVisibility','OFF');
       
    set(handles.status,'String',['Validation vs manual segmentation finished']);
    disp('=======================================');


% function edit3_Callback(hObject, eventdata, handles)
% % hObject    handle to edit3 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit3 as text
% %        str2double(get(hObject,'String')) returns contents of edit3 as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit3_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit3 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end




function imgPath_Callback(hObject, eventdata, handles)
% hObject    handle to imgPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgPath as text
%        str2double(get(hObject,'String')) returns contents of imgPath as a double


% --- Executes during object creation, after setting all properties.
function imgPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roiPath_Callback(hObject, eventdata, handles)
% hObject    handle to roiPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roiPath as text
%        str2double(get(hObject,'String')) returns contents of roiPath as a double


% --- Executes during object creation, after setting all properties.
function roiPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function savedPath_Callback(hObject, eventdata, handles)
% hObject    handle to savedPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savedPath as text
%        str2double(get(hObject,'String')) returns contents of savedPath as a double


% --- Executes during object creation, after setting all properties.
function savedPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savedPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eta_Callback(hObject, eventdata, handles)
% hObject    handle to eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eta as text
%        str2double(get(hObject,'String')) returns contents of eta as a double


% --- Executes during object creation, after setting all properties.
function eta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_Callback(hObject, eventdata, handles)
% hObject    handle to T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T as text
%        str2double(get(hObject,'String')) returns contents of T as a double


% --- Executes during object creation, after setting all properties.
function T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function its_Callback(hObject, eventdata, handles)
% hObject    handle to its (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of its as text
%        str2double(get(hObject,'String')) returns contents of its as a double


% --- Executes during object creation, after setting all properties.
function its_CreateFcn(hObject, eventdata, handles)
% hObject    handle to its (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outputPath_Callback(hObject, eventdata, handles)
% hObject    handle to outputPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputPath as text
%        str2double(get(hObject,'String')) returns contents of outputPath as a double


% --- Executes during object creation, after setting all properties.
function outputPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dz_Callback(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dz as text
%        str2double(get(hObject,'String')) returns contents of dz as a double


% --- Executes during object creation, after setting all properties.
function dz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double


% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_Callback(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta as text
%        str2double(get(hObject,'String')) returns contents of beta as a double


% --- Executes during object creation, after setting all properties.
function beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
