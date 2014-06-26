

function for_all ()

% Script written by Yuval:

ControlsFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i=1:length(nameFolds)
MRI2Afni (nameFolds{i}, ControlsFolder)

end
end

function MRI2Afni (PatientName, PatientPath)

PathName = strcat(PatientPath, filesep, PatientName)
cd (PathName)

if ~exist (strcat('s_', PatientName, '.nii'), 'file')
    return
end

eval(['!3dcopy ', 's_', PatientName, '.nii ', 'anat+orig'])
% !3dcopy s_Pat_03_13014bg.nii anat+orig


end