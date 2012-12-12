function image2mesh_gui_tab()

data = struct();
gui = struct();
createGUI();

    function createGUI()
        f = figure();

        % Tab layout and first panel
        tabLayout = uiextras.TabPanel('Parent',f);
        % Main layout for tab 1
        basicLayout = uiextras.VBox('Parent',tabLayout);
        tabLayout.TabNames = {'Basic'};
        tabLayout.SelectedChild = 1;
        tabLayout.Padding = 5;

        % Add Sub layouts
        % Input files panel
        inputPanel = uiextras.Panel('Title','Input info','Parent',basicLayout);
        inputFilenamesLayout = uiextras.HBox('Parent', inputPanel);
        set(inputFilenamesLayout, 'Padding', 4);

        imageFilenameLayout = uiextras.VBox('Parent', inputFilenamesLayout);
        gui.ImageFilenameText = uicontrol('Parent',imageFilenameLayout,...
            'Style','text','HorizontalAlignment','left','String','Image/Segmentation File name:');
        imageFilenameSubLayout = uiextras.HBox('Parent',imageFilenameLayout); 
        gui.ImageFilenameEdit = uicontrol('Parent',imageFilenameSubLayout,'Style','edit');
        gui.BrowseInputImageButton = uicontrol('parent', imageFilenameSubLayout, 'Style', 'pushbutton','String','Browse...');
        set(gui.BrowseInputImageButton, 'Callback', @onBrowseInputImageClicked);
        set(imageFilenameSubLayout, 'Sizes', [-1 70]);
        set(imageFilenameLayout, 'Sizes', [20 40]);

        sdFilenameLayout = uiextras.VBox('Parent', inputFilenamesLayout);
        gui.SDFilenameText = uicontrol('Parent',sdFilenameLayout,'Style','text',...
            'HorizontalAlignment','left','String','Source/Detector File name:');
        sdFilenameSubLayout = uiextras.HBox('Parent',sdFilenameLayout);
        gui.SDFilenameEdit = uicontrol('Parent',sdFilenameSubLayout,'Style','edit');
        gui.BrowseSDButton = uicontrol('parent', sdFilenameSubLayout, 'Style', 'pushbutton','String','Browse...');
        set(gui.BrowseSDButton, 'Callback', @onBrowseSDButtonClicked);
        set(sdFilenameSubLayout, 'Sizes', [-1 70]);
        set(sdFilenameLayout, 'Sizes', [20 40]);
        %% Image size panels
        secondRowLayout = uiextras.HBox('Parent', basicLayout);
        set(secondRowLayout, 'Padding', 5);

        % Image info panel
        mainImageInfoPanel = uiextras.Panel('Title', 'Image Info', 'Parent', secondRowLayout);
        mainImageInfoLayout = uiextras.VBox('Parent',mainImageInfoPanel);
        set(mainImageInfoLayout,'Padding',5);

        panelTop = uiextras.Panel('Title','Pixel Size (mm)','Parent',mainImageInfoLayout);
        imageSizesLayout = uiextras.VBox('Parent',panelTop);
        pixelInfoLayout = uiextras.Grid('Parent', imageSizesLayout, 'Spacing', 5);
        gui.XText = uicontrol('Parent', pixelInfoLayout, 'Style', 'text','HorizontalAlignment','center','String','X');
        gui.XEdit = uicontrol('Parent', pixelInfoLayout, 'Style', 'edit');
        gui.YText = uicontrol('Parent', pixelInfoLayout, 'Style', 'text','HorizontalAlignment','center','String','Y');
        gui.YEdit = uicontrol('Parent', pixelInfoLayout, 'Style', 'edit');
        gui.ZText = uicontrol('Parent', pixelInfoLayout, 'Style', 'text','HorizontalAlignment','center','String','Z');
        gui.ZEdit = uicontrol('Parent', pixelInfoLayout, 'Style', 'edit');
        set(pixelInfoLayout, 'ColumnSizes', [-1 -1 -1], 'RowSizes', [20 -1]);

        panelButtom = uiextras.Panel('Title','Voxel Info','Parent',mainImageInfoLayout);
        voxelSizesLayout = uiextras.VBox('Parent',panelButtom);
        voxelInfoLayout = uiextras.Grid('Parent', voxelSizesLayout, 'Spacing', 5);
        gui.VXText = uicontrol('Parent', voxelInfoLayout, 'Style', 'text','HorizontalAlignment','center','String','Rows');
        gui.VXEdit = uicontrol('Parent', voxelInfoLayout, 'Style', 'edit');
        gui.VYText = uicontrol('Parent', voxelInfoLayout, 'Style', 'text','HorizontalAlignment','center','String','Columns');
        gui.VYEdit = uicontrol('Parent', voxelInfoLayout, 'Style', 'edit');
        gui.VZText = uicontrol('Parent', voxelInfoLayout, 'Style', 'text','HorizontalAlignment','center','String','Slices');
        gui.VZEdit = uicontrol('Parent', voxelInfoLayout, 'Style', 'edit');
        set(voxelInfoLayout, 'ColumnSizes', [-1 -1 -1], 'RowSizes', [20 -1]);

        % Status text panel
        statusPanel = uiextras.Panel('Parent',secondRowLayout,'Title','Status');
        statusTextLayout = uiextras.VBox('Parent',statusPanel);
        set(statusTextLayout, 'Padding', 5);
        gui.StatusText = uicontrol('Parent',statusTextLayout,'Style','text');

        %% Mesh Sizing
        meshSizePanel = uiextras.Panel('Title', 'Mesh Sizing', 'Parent', basicLayout);
        meshSizeLayout = uiextras.Grid('Parent',meshSizePanel, 'Spacing', 5);
        gui.CellSizeText = uicontrol('Parent',meshSizeLayout,'Style','text','HorizontalAlignment','center','String','Global Tetrahedron Size');
        gui.CellSizeEdit = uicontrol('Parent',meshSizeLayout,'Style','edit');
        gui.RegionIDsText = uicontrol('Parent',meshSizeLayout,'Style','text','HorizontalAlignment','center','String','Region IDs to refine');
        gui.RegionIDsEdit = uicontrol('Parent',meshSizeLayout,'Style','edit');
        gui.RegionSizesText = uicontrol('Parent',meshSizeLayout,'Style','text','HorizontalAlignment','center','String','Region Sizes');
        gui.RegionSizesEdit = uicontrol('Parent',meshSizeLayout,'Style','edit');
        set(meshSizeLayout, 'ColumnSizes', [-1 -1 -1], 'RowSizes', [20 35]);

        %%
        outputPanel = uiextras.Panel('Title','Output Settings','Parent', basicLayout,'Padding',5);
        outputLayout = uiextras.Grid('Parent',outputPanel, 'Spacing',5);
        gui.OutputFilenameText = uicontrol('Parent',outputLayout,'Style','text','HorizontalAlignment','center','String','Output file name:');
        gui.OutputFilenameEdit = uicontrol('Parent',outputLayout,'Style','edit');
        uiextras.Empty('Parent',outputLayout);
        gui.OutputFilenameButton = uicontrol('Parent',outputLayout,'Style','pushbutton','String','Browse...');
        gui.NirfastMeshTypeText = uicontrol('Parent',outputLayout,'Style','text','HorizontalAlignment','center','String','Nirfast mesh type');
        gui.NirfastMeshTypePopup = uicontrol('Parent',outputLayout,'Style','popupmenu',...
            'String',{'Standard','Fluorescence','Spectral','SPN'});
        uiextras.Empty('Parent',outputLayout);
        gui.RunMeshGenerator = uicontrol('Parent',outputLayout,'Style','pushbutton','String','Run Mesh Generator...');
        set(outputLayout, 'ColumnSizes', [200 70 120 -1], 'RowSizes', [20 40]);

        set(f,'Position', [440   122   649   460]);
        set(basicLayout, 'Sizes', [80 -1 80 80]);
    end
%%
    function onBrowseInputImageClicked(~, ~)
        [fn, pathname] = uigetfile( ...
            {'*.bmp;*.jpg;*.tif;*.gif;*.mha;*.mhd','Image Files (*.bmp,*.jpg,*.tif,*.gif,*.mha,*.mhd)';'*.*','All Files (*.*)'}, ...
           'Pick a file');
        if isequal(fn,0)
            error('You need to select an image file!');
        end
        UpdateInputFileInfo(fullfile(pathname,fn));
    end
    function UpdateInputFileInfo(inputfn)
        flag=false;
        if isfield(data,'infilename') && ~strcmp(inputfn,data.infilename) || ...
                ~isempty(get(gui.ImageFilenameEdit, 'String'))
            flag=true;
        end
        set(gui.ImageFilenameEdit, 'String', inputfn);
        data.infilename = inputfn;
        
        set(gui.StatusText, 'String', {'Status:';'';'Loading Image';'Please wait...'});
        drawnow('update');
        UpdateImageInformation();
        UpdateOutputFn();
        if flag
            [f1 f2] = fileparts(data.infilename);
            if exist(fullfile(f1,[f2 '.txt']),'file')
                data.sdfilename = fullfile(f1,[f2 '.txt']);
                UpdateSDFileInfo();
            elseif exist(fullfile(f1,[f2 '.csv']),'file')
                data.sdfilename = fullfile(f1,[f2 '.csv']);
                UpdateSDFileInfo();
            else
                set(gui.SDFilenameEdit, 'String', 'Can not find S/D coordinate info!');
                set(gui.SDFilenameEdit, 'ForegroundColor', [1 0 0]);
                data.sdfilename = '';
            end
        end
    end
    function UpdateImageInformation()
        param.medfilter = 1;
        [mask info] = GetImageStack(data.infilename, param);
        regions = unique(mask(:));
        s = {''};
        if ~isempty(info)
            s{1} = 'Info from file header:';
            s{2} = sprintf('Rows: %d, Cols: %d, Slices: %d', info.Dimensions(1), ...
                info.Dimensions(2), info.Dimensions(3));
            if isfield(info,'PixelDimensions')
                set(gui.XEdit, 'String', num2str(info.PixelDimensions(1)));
                set(gui.YEdit, 'String', num2str(info.PixelDimensions(2)));
                set(gui.ZEdit, 'String', num2str(info.PixelDimensions(3)));
                s{end+1} = sprintf('Pixel size: [%.3f %.3f %.3f]',info.PixelDimensions);
                fook = min(info.PixelDimensions(1),info.PixelDimensions(2));
                UpdateTetAndFacetSize(min(info.PixelDimensions(3),fook));
            else
                set(gui.XEdit, 'String', '');
                set(gui.XEdit, 'String', '');
                set(gui.XEdit, 'String', '');
            end
            if isfield(info,'Dimensions')
                set(gui.VXEdit, 'String', num2str(info.Dimensions(1)));
                set(gui.VYEdit, 'String', num2str(info.Dimensions(2)));
                set(gui.VZEdit, 'String', num2str(info.Dimensions(3)));
            end
            if isfield(info,'Offset')
                s{end+1} = sprintf('Offset: [%.2f %.2f %.2f]', info.Offset);
            end
            set(gui.StatusText, 'String', s)
        end
        s{end+1} = '';
        s{end+1} = sprintf('Region IDs (total: %d):', length(regions));
        fmt = ['%d,' repmat(' %d,', 1, length(regions)-1)];
        s{end+1} = sprintf(fmt, regions);
        set(gui.StatusText, 'String', s);
        [tf idx] = ismember(0, regions);
        if tf, regions(idx) = []; end
        
        foo = get(gui.RegionIDsEdit, 'String');
        if isempty(foo) || (numel(foo)==1 && str2double(foo) == 0) || ...
                (isfield(data,'regions') && length(regions) ~= length(data.regions_torefine))
            set(gui.RegionIDsEdit, 'String', '0');
            data.regions_torefine = 0;
        end
        
        if isempty(get(gui.RegionSizesEdit, 'String')) ||...
                (isfield(data,'regions') && length(regions) ~= length(data.regions_torefine))
            set(gui.RegionSizesEdit, 'String', '0');
            data.region_sizes = 0;
        end
        data.regions_torefine = regions;
        data.mask = uint8(mask);
        data.maskinfo = info;
    end
    function UpdateOutputFn()
        s = remove_extension(data.infilename);
        s = add_extension(s, '.ele');
        data.outputfn = s;
        set(gui.OutputFilenameEdit, 'String', data.outputfn);
    end
    function UpdateTetAndFacetSize(minsize)
        % set(gui.facet_size , 'String', num2str(2*minsize));
        set(gui.CellSizeEdit , 'String', num2str(2*minsize));
        data.cell_size = minsize;
    end
    function onBrowseSDButtonClicked(~, ~)
        [fn, pathname] = uigetfile( ...
            {'*.txt;*.csv','Text Files (*.txt,*.csv)';'*.*','All Files (*.*)'}, ...
           'Pick a file');
        if isequal(fn,0)
            warning('image2mesh_gui:NoInputFile','You need to select an image file!');
        end
        data.sdfilename = fullfile(pathname,fn);
        UpdateSDFileInfo();
    end
    function UpdateSDFileInfo(hObject,eventdata,handles)
        set(gui.SDFilenameEdit, 'ForegroundColor', [0 0 0]);
        fid=fopen(data.sdfilename,'rt');
        s=textscan(fid,'%f,%f,%f,%f');
        data.sdcoords = [s{2} s{3} s{4}];
        s = get(gui.StatusText, 'String');
        s{end+1}='';
        s{end+1} = sprintf('Number of source/detectors: %d',size(data.sdcoords,1));
        set(gui.StatusText ,'String',s);
    end
end