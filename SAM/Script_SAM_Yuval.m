


% Create Parameter File for all Subjects:

% !!  Automatisierung für Patienten siehe kh_Sam_Beamformer

createPARAM('M400','ERF','VG',[0.35 0.45],'VG',[-0.1 0],[1 50],[-0.1 0.7]);
% if not done yet, do:
!meshnorm ortho_brainhull.ply > hull.shape

% then compute covariance:
!SAMcov64 -r  zzz_ms -d hb_lf_c,rfhp0.1Hz -m M400 -v

% then compute weights:
!SAMwts64 -r  zzz_ms -d hb_lf_c,rfhp0.1Hz -m M400 -c VGa -v

%% did not work, maybe leave it out: ERPs
% !SAMerf64 -r  zzz_ht -d hb_tr_lf_c,rfhp0.1Hz -m M400 -v
% !SAMerf64 -r  zzz_ht -d hb_tr_lf_c,rfhp0.1Hz -m M400 -z 1 -v
% !afni &
%%
% change matlab dir to Path of weights (SAMdir):
% multiply weights with avgBL:
[SAMHeader, ActIndex, ActWgts]=readWeights('M400,1-50Hz,VGa.wts');
load('avgBL.mat')

vs=ActWgts*avgBL.avg;
ns=mean(abs(ActWgts),2);
vs=vs./repmat(ns,1,size(vs,2));

cfg=[];
cfg.step=5;
cfg.boxSize=[-120 120 -90 90 -20 150];
cfg.prefix='ERF'; % change prefix
cfg.torig=-500;   %  comment if you want to sum up activity of specific time intervall
cfg.TR=1/1.01725; % comment if you want to sum up activity of specific time intervall
% VS2Brik(cfg,vs);
% max(max(vs))
VS2Brik(cfg,1e+13*abs(vs)); % =>creates ERF+orig.Brik+Head 
save Workspace_SAM

!afni &