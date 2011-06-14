function [e p] = RunCGALMeshGenerator(mask,param)
% Runs CGAL mesh generator on 3D matrix 'mask'
% For info on 'param' and 'info' refer to 'image2mesh_cgal.m'
% 
% Written by:
%   Hamid Ghadyani May 2011

tmpmeshfn = '._out.mesh';
tmpinrfn  = '._cgal_mesher.inr';
savefn = add_extension(tmpinrfn,'.inr');
saveinr(mask,savefn,param);

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

% Run the executable
syscommand = GetSystemCommand('image2mesh_cgal');
if ~ispc
    eval(['! chmod u+x "' syscommand '"']);
end
makemeshcommand = ['! "' syscommand '" ' savefn ' ' cgalparam_fn ' ' tmpmeshfn];
eval(makemeshcommand);

% Read the resulting mesh
[e p] = readMEDIT(tmpmeshfn);
if isfield(param,'Offset')
    p = p + repmat(param.Offset,size(p,1),1);
end
warning('off','MATLAB:DELETE:FileNotFound');
delete(cgalparam_fn,tmpmeshfn,tmpinrfn);
warning('on','MATLAB:DELETE:FileNotFound');

