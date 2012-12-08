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
set(imageFilenameSubLayout, 'Sizes', [-1 70]);
set(imageFilenameLayout, 'Sizes', [20 40]);

sdFilenameLayout = uiextras.VBox('Parent', inputFilenamesLayout);
gui.SDFilenameText = uicontrol('Parent',sdFilenameLayout,'Style','text',...
    'HorizontalAlignment','left','String','Source/Detector File name:');
sdFilenameSubLayout = uiextras.HBox('Parent',sdFilenameLayout);
gui.SDFilenameEdit = uicontrol('Parent',sdFilenameSubLayout,'Style','edit');
gui.BrowseSDButton = uicontrol('parent', sdFilenameSubLayout, 'Style', 'pushbutton','String','Browse...');
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
% set(basicLayout, 'Sizes', [80 180 80 80]);

% genInfoLayout = uiextras.HBox('Parent', inputImageInfoLayout);
% gui.GenInfoText = uicontrol('Parent', genInfoLayout, 'Style', 'text');
