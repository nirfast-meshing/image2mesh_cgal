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
% outfn: (optional) specifies the prefix for .ele/.node that the resulting
%        mesh will be written into

if nargin==0
    [fn, pathname] = uigetfile( ...
    {'*.bmp;*.jpg;*.tif;*.gif','Image Files';
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

mask = GetImageStack(fn,param);

[nrow ncol nslice] = size(mask);
mask = uint8(mask);

medfilter = 0;
if isfield(param,'medfilter') && param.medfilter==1
    medfilter = 1;
end
for i=1:nslice
    foo = mask(:,:,i);
    if medfilter == 1
        foo = medfilt2(foo,[5 5]);
    end
    mask(:,:,i) = foo;
end
stackInfo.PixelSpacing(1) = param.xpixelsize;
stackInfo.PixelSpacing(2) = param.ypixelsize;
stackInfo.SliceThickness = param.zpixelsize;

fn = remove_extension(fn);
savefn = add_extension(fn,'.inr');
saveinr(mask,savefn,stackInfo);

% Write the meshing parameters file
facet_angle = 25; facet_size = 3; facet_distance = 2;
cell_radius_edge = 3; cell_size = 3; % general tet size of all regions
special_subdomain_label = 0; % label of region to be refined
special_size = 0; % tet size of the region 'special_subdomain_label'
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


[e p] = readMEDIT(tmpmeshfn);
if nargin<3
    outfn=[fn '-tetmesh'];
end
outfn = add_extension(outfn,'.ele');
writenodelm_nod_elm(outfn,e,p,[],1);
warning('off','MATLAB:DELETE:FileNotFound');
delete(cgalparam_fn,tmpmeshfn);
warning('on','MATLAB:DELETE:FileNotFound');
