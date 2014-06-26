
function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
   kh_ERF_transform_MNI( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .47)
%    kh_ERF_transform_MNI( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .31) 
%    kh_ERF_transform_MNI( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .6) 
%    
   
end


end


function kh_ERF_transform_MNI (SubjectPath, SubjectName, TimeBeg, TimeEnd)

% if 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ca') || 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ht')
%     return
% end

SAMPath = strcat(SubjectPath, filesep, 'SAM');
cd (SAMPath)

FileNameOldBRIK = strcat('ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc.BRIK');
if exist (FileNameOldBRIK, 'file')
    delete (FileNameOldBRIK)
end

FileNameOldHEAD = strcat('ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc.HEAD');
if exist (FileNameOldHEAD, 'file')
    delete (FileNameOldHEAD)
end



FileName = strcat('ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+orig');
%  disp(['!@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ', NewFileName,' -dxyz 5']) % 
eval(['!@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ', FileName,' -dxyz 5']) % 

end