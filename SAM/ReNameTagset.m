
function for_all ()

% Script written by Yuval:

ControlsFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i=1:length(nameFolds)
RenameTagfile (nameFolds{i}, ControlsFolder)

end
end

function RenameTagfile (PatientName, PatientPath)

PathName = strcat(PatientPath, filesep, PatientName)
cd (PathName)

if ~exist (strcat('MNI305.tag'), 'file')
    return
end

OldName = 'MNI305.tag';
NewName = (strcat( PatientName, '.tag'))

movefile(OldName,NewName)
%  eval(['!rename ' OldName ' ' NewName])



end