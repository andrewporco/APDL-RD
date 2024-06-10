clc
clear
startx = 0; %mm
starty = 0; %mm
startz = 0; %mm
endx = 30; %mm
endy = 60; %mm
endz = 30; %mm
gridX = (startx:0.2:startx+endx); %22
gridY = (starty:0.2:starty+endy); %96
gridZ = (startz:0.2:startz+endz); %63

[gridOUTPUT, gridCOx, gridCOy, gridCOz] = VOXELISE(gridX, gridY, gridZ, 'VideoPart.stl');
% Create fxn to save grid

[vol_handle] = VoxelPlotter(gridOUTPUT,1);

