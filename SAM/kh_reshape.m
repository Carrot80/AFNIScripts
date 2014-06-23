
% um aus zwei dimensionen 4d zu machen
% aus Yuval's VS2Brik.m entnommen
% nochmal ausprobieren, ob ich auf gleiches Ergebnis komme, wenn ich
% WriteBrik mit rehape bzw. VS2Brik ohne Rehape verwende 
cfg.boxSize = [-120 120 -90 90 -20 150];
cfg.step    = 5;

xyzMin=cfg.boxSize([1 3 5]);
xyzMax=cfg.boxSize([2 4 6]);
xsize=length(xyzMin(1):cfg.step:xyzMax(1));
ysize=length(xyzMin(2):cfg.step:xyzMax(2));
zsize=length(xyzMin(3):cfg.step:xyzMax(3));

tsize=size(vs,2);
vsRs=reshape(vs,[zsize,ysize,xsize,tsize]); % 35*37*49

% 
[V_ERF, Info_ERF]=BrikLoad('ERF_1000ms_Trial_1+orig'); % 49*37*35

% unklar, ob ich origin setzen muss (siehe VS2Brik)