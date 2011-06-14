function LaunchNirfastMesher(inrfn)

inrfn=inrfn(2:end-1);
h=image2mesh_gui;
data=guidata(h);
data.inputfn=inrfn;
guidata(h,data);
image2mesh_gui('UpdateInputFileInfo',h,[],data);
data=guidata(h);
image2mesh_gui('callimage2mesh_cgal_Callback',h,[],data);
