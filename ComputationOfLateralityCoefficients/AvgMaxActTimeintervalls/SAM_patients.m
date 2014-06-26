
function for_all_pat()

Patientpath = '/home/kh/ShareWindows/data/patients/patients_SAM';

Dir = dir(Patientpath);
isub = [Dir(:).isdir]; 
nameFolds = {Dir(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

% for i=21:length(nameFolds)
%  for i=[2:5, 8, 9]
for i=19

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

[fileName]=findFileName(PatientPath, PatientName)


% if isempty (fileName)
%     return
% end
% 
% if exist (strcat(PatientPath, filesep, PatientName, filesep, 'SAM', filesep, 'M400,1-50Hz,VGa.wts'))
%    return
% end
% 
% if exist (strcat(PatientPath, filesep, PatientName, filesep, 'SAM'), 'dir')
%     % compute weights:
%     eval(['!SAMwts64 -r ', PatientName, ' -d ', fileName, ' -m M400 -c VGa -v'])
%     return
% end

% covariance
eval (['!SAMcov64 -r ',  PatientName, ' -d ', fileName, ' -m M400 -v'])

% compute weights: 
eval(['!SAMwts64 -r ', PatientName, ' -d ', fileName, ' -m M400 -c VGa -v'])



 
end


function [fileName]=findFileName(PatientPath, PatientName)

   
     PathfileNameHBcor=strcat(PatientPath, filesep,PatientName, filesep, 'hb_tr_lf_c,rfhp0.1Hz') ;
     PathfileNameHBcorNoTr=strcat(PatientPath, filesep, PatientName, filesep, 'hb_lf_c,rfhp0.1Hz') ;
     PathfileNameNoHBcor=strcat(PatientPath, filesep, PatientName,filesep,  'tr_lf_c,rfhp0.1Hz') ;
     PathfileNameLfcor=strcat(PatientPath, filesep, PatientName, filesep, 'lf_c,rfhp0.1Hz') ;

     if exist (PathfileNameHBcor, 'file')
         fileName='hb_tr_lf_c,rfhp0.1Hz';
     elseif exist (PathfileNameHBcorNoTr, 'file')
         fileName='hb_lf_c,rfhp0.1Hz';
     elseif exist(PathfileNameNoHBcor, 'file')
         fileName='tr_lf_c,rfhp0.1Hz';
     elseif exist(PathfileNameLfcor, 'file')
         fileName='lf_c,rfhp0.1Hz';
         
     else 
         
         fileName = [];
     end

end
