function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
TimeInt = [.32, .6];

for i= 1:size(nameFolds)
    
   kh_convert2nifti( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, TimeInt)
    
end


end

function kh_convert2nifti(SubjectPath, SubjectName, TimeInt)


PathName = strcat(SubjectPath, filesep, 'SAM');
cd (PathName)

PathERF = strcat('br_z_transf_brain01ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, '+tlrc');

% eval(['!@auto_tlrc -apar ', SubjectPath, 'orthoMNI_avg152T+tlrc',' -input ', PathERF, ' -dxyz 5']);

eval (['!3dcopy ', PathERF, ' ',strcat('br_z_transf_brain01ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, 'MNI.nii')])
% eval (['!3dcopy ', strcat('ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, '+tlrc'), ' ',strcat('ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, '_MNI.nii')])


PathKeptTrials=strcat(SubjectPath, filesep, 'keptTrials')
cd(PathKeptTrials)

UTest= strcat('brain01_Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 'sMNI+tlrc')
eval (['!3dcopy ', UTest, ' ', strcat('brain01_Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 'sMNI.nii')])
% eval (['!3dcopy ', strcat('brain01_Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's+orig'), ' ', strcat('brain01_Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's.nii')])

eval (['!3dcopy ', 'orthoMNI_avg152T+tlrc ', 'orthoMNI_avg152T.nii'])

cd (PathName)
eval (['!3dcopy ', 'ortho+orig ', 'ortho.nii'])
end