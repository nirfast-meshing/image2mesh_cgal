function [e p] = image2mesh_cgal(fn,param,outfn)
% Tries to read stack of 2D images (with file name 'fn') and create a
% tetrahedral mesh using CGAL library. It returns the mesh in 'e' and 'p',
% tetrahedral elements and node coordinates.
% 
% param.pad : add padding to 4 sides of each 2D image
% param.medfilter : apply a median filter to iron out speckle noises of images
% param.xpixelsize = x pixel size
% param.ypixelsize = y pixel size
% param.zpixelsize = z pixel size
% 
% subdomain label (defined above)
% param.facet_angle = minimum angle of triangles on surface of mesh (25)
% param.facet_size = controls size of triangles on surface of mesh
% param.facet_distance = controls how accurate the surface mesh mimics the
%                        actual model
% param.cell_radius_edge = controls the quality of tetrahedrons
% param.cell_size = controls the general size of tetrahedrons
% 
% param.special_subdomain_label = label ID (grayscale value) of subdomain
%                    that user wants to have a different tetrahedron size
% param.special_subdomain_size = size of tetrahedron for the special
% 
% outfn: (optional) specifies the prefix for .ele/.node that the resulting
%        mesh will be written into

if nargin==0
    [fn, pathname] = uigetfile( ...
    {'*.bmp;*.jpg;*.tif;*.gif;*.mha','Image Files';
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file');
    if isequal(fn,0)
        error(sprintf('\nYou need to select an image file!'));
    end
    fn = fullfile(pathname,fn);
    param.xpixelsize = 0.703;
    param.ypixelsize = 0.703;
    param.zpixelsize = 1.5;
    param.pad = 0;
    param.medfilter = 1;
end

[mask info] = GetImageStack(fn,param);

[nrow ncol nslice] = size(mask);
mask = uint8(mask);

medfilter = 0;
if isfield(param,'medfilter') && param.medfilter==1
    medfilter = 1;
end
for i=1:nslice
    foo = mask(:,:,i);
%     foo(foo==198)=50;
%     foo(foo==199)=200;
    if medfilter == 1
        foo = medfilt2(foo,[5 5]);
    end
    mask(:,:,i) = foo;
end

if isfield(info,'PixelDimensions')
    stackInfo.PixelSpacing(1) = info.PixelDimensions(1);
    stackInfo.PixelSpacing(2) = info.PixelDimensions(2);
    stackInfo.SliceThickness  = info.PixelDimensions(3);
elseif isfield(param,'xpixelsize')
    stackInfo.PixelSpacing(1) = param.xpixelsize;
    if isfield(param,'ypixelsize'), stackInfo.PixelSpacing(2) = param.ypixelsize; end
    if isfield(param,'zpixelsize'), stackInfo.SliceThickness  = param.zpixelsize; end
else
    error('No pixel dimension information is provided');
end

savefn = add_extension(fn,'.inr');
saveinr(mask,savefn,stackInfo);

% Set up the necessary parameters for meshing
facet_angle = 25; facet_size = 3; facet_distance = 2;
cell_radius_edge = 3; cell_size = 3; % general tet size of all regions
special_subdomain_label = 0; % label of region to be refined
special_size = 0; % tet size of the region 'special_subdomain_label'
if isfield(param,'facet_angle'), facet_angle = param.facet_angle; end
if isfield(param,'facet_size'),  facet_size  = param.facet_size; end
if isfield(param,'facet_distance'), facet_distance = param.facet_distance; end
if isfield(param,'cell_radius_edge'), cell_radius_edge = param.cell_radius_edge; end
if isfield(param,'cell_size'), cell_size = param.cell_size; end
if isfield(param,'special_subdomain_label'), special_subdomain_label = param.special_subdomain_label; end
if isfield(param,'special_size'), special_size = param.special_size; end

% Write up the parameter files
cgalparam_fn = [pwd filesep 'criteria.txt'];
fid = fopen(cgalparam_fn,'wt');
fprintf(fid,'%f\n',facet_angle);
fprintf(fid,'%f\n',facet_size);
fprintf(fid,'%f\n',facet_distance);
fprintf(fid,'%f\n',cell_radius_edge);
fprintf(fid,'%f\n',cell_size);
fprintf(fid,'%d\n',special_subdomain_label);
fprintf(fid,'%f\n',special_size);
fclose(fid);
tmpmeshfn='out.mesh';
% Run the executable
syscommand = GetSystemCommand('image2mesh_cgal');
if ~ispc
    eval(['! chmod u+x "' syscommand '"']);
end
makemeshcommand = ['! "' syscommand '" ' savefn ' ' cgalparam_fn ' ' tmpmeshfn];
eval(makemeshcommand);

% Read the resulting mesh
[e p] = readMEDIT(tmpmeshfn);
if isfield(info,'Offset')
    p = p + repmat(info.Offset,size(p,1),1);
end
if nargin<3
    outfn=[fn '-tetmesh'];
end
outfn = add_extension(outfn,'.ele');
writenodelm_nod_elm(outfn,e,p,[],1);
warning('off','MATLAB:DELETE:FileNotFound');
delete(cgalparam_fn,tmpmeshfn);
warning('on','MATLAB:DELETE:FileNotFound');

