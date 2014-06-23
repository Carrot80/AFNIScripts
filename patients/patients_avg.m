
function for_all_pat()

Patientpath = '/home/kh/ShareWindows/data/patients/patients_SAM';

Dir = dir(Patientpath);
isub = [Dir(:).isdir]; 
nameFolds = {Dir(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i=6:length(nameFolds)

kh_SAMpatients (Patientpath, nameFolds{i} )

end


end


function kh_SAMpatients (PatientPath, PatientName)
% 
% if 0 == strcmp (PatientName, 'Pat_19_13055eg')
%     
%     return 
%     
% end

SAMPath = strcat(PatientPath, filesep, PatientName, filesep, 'SAM');
cd (SAMPath)

% to multiply weights with avgBL, read weights :
[SAMHeader, ActIndex, ActWgts]=readWeights('M400,1-50Hz,VGa.wts');

% load avg:
PathAVG = strcat(PatientPath, filesep, PatientName, filesep, 'avgBL');
load(PathAVG)

vs=ActWgts*avgBL.avg;
% dies nächsten beiden Zeilen auskommentieren und
% schauen, ob Ergebnisse besser sind, da noisenormalisierung evtl.
% Kopfposition mit in Sourcespace übernimmt; Bias to Center ist dann zwar
% vorhanden, UTest sollte aber dennoch funktionieren bzw. Lateralisierung
ns=mean(abs(ActWgts),2); 
vs=vs./repmat(ns,1,size(vs,2));



cd(strcat(PatientPath, filesep, PatientName))

cfg=[];
cfg.step=5;
cfg.boxSize=[-120 120 -90 90 -20 150];
cfg.prefix='ERF_noise'; % change prefix
cfg.torig=-500;   %  comment if you want to sum up activity of specific time intervall
cfg.TR=1/1.01725; % comment if you want to sum up activity of specific time intervall
% VS2Brik(cfg,vs);
% max(max(vs))
VS2Brik(cfg,1e+13*abs(vs)); % =>c
 
end


