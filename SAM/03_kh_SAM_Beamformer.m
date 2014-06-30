function ForAll()

ControlsFolder = '/home/kh/data/controls_SAM';

createPARAM('M400','ERF','VG',[0.35 0.45],'VG',[-0.1 0],[1 50],[-0.1 0.7]);

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %# returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
   kh_SAM( strcat(ControlsFolder, filesep, nameFolds(i,1)), nameFolds(i))

end

end




function kh_SAM(SubjectPath, SubjectName)

Path_vs = (strcat(SubjectPath, filesep, SAM, filesep, 'Workspace_SAM'));

if exist(Path_vs, 'file')
    return
end


[FileName]=findFileName(PatientPath, PatientName)

PathAll = '/home/kh/data/controls_SAM/';

cd (PathAll)

% compute covariance
!SAMcov64 -r  SubjectName -d FileName -m M400 -v

% then compute weights:
!SAMwts64 -r  SubjectName -d FileName.1Hz -m M400 -c VGa -v


% change matlab dir to Path of weights (SAMdir):
PathSAM = strcat(SubjectPath, filesep, SAM, filesep);

% to multiply weights with avgBL, read weights :
[SAMHeader, ActIndex, ActWgts]=readWeights('M400,1-50Hz,VGa.wts');

% load avg:
PathAVG = strcat(SubjectPath, filesep, 'avgBL_', PatientName);
load(PathAVG)

vs=ActWgts*avgBL.avg;
% dies nächsten beiden Zeilen auskommentieren und
% schauen, ob Ergebnisse besser sind, da noisenormalisierung evtl.
% Kopfposition mit in Sourcespace übernimmt; Bias to Center ist dann zwar
% vorhanden, UTest sollte aber dennoch funktionieren bzw. Lateralisierung
% ns=mean(abs(ActWgts),2); 
% vs=vs./repmat(ns,1,size(vs,2));

cfg=[];
cfg.step=5;
cfg.boxSize=[-120 120 -90 90 -20 150];
cfg.prefix='ERF'; % change prefix
cfg.torig=-500;   %  comment if you want to sum up activity of specific time intervall
cfg.TR=1/1.01725; % comment if you want to sum up activity of specific time intervall
% VS2Brik(cfg,vs);
% max(max(vs))
VS2Brik(cfg,1e+13*abs(vs)); % =>creates ERF+orig.Brik+Head 


save (strcat(SubjectPath, filesep, SAM, filesep, 'Workspace_SAM')); % oder nur einzelne files speichern

end

function [FileName]=findFileName(SubjectPath, SubjectName)

   
     PathfileNameHBcor=strcat(SubjectPath, filesep,'hb_tr_lf_c,rfhp0.1Hz') ;
     PathfileNameHBcorNoTr=strcat(SubjectPath, filesep,'hb_lf_c,rfhp0.1Hz') ;
     PathfileNameNoHBcor=strcat(SubjectPath, filesep,'tr_lf_c,rfhp0.1Hz') ;
     PathfileNameLfcor=strcat(SubjectPath, filesep,'lf_c,rfhp0.1Hz') ;

     if exist (PathfileNameHBcor, 'file')
         fileName=PathfileNameHBcor;
     elseif exist (PathfileNameHBcorNoTr, 'file')
         fileName=PathfileNameHBcorNoTr;
     elseif exist(PathfileNameNoHBcor, 'file')
         fileName=PathfileNameNoHBcor;
     elseif exist(PathfileNameLfcor, 'file')
         fileName=PathfileNameLfcor;
         
     else 
         
         fileName = [];
     end

end
