function patients ()
   
    PatientFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';
    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

    forAll (PatientFolder, 'Patients')
    forAll (ControlsFolder, 'Controls')

end


function forAll(Folder, group)
    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    
    TimeInt=[.2, .6];
    for i= 1:size(nameFolds)
       kh_deleteFiles(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, TimeInt)
    end    
end

function kh_deleteFiles(SubjectPath, SubjectName, TimeInt)
if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg')
    return
end
cd(strcat(SubjectPath, filesep, 'RelativeAct'))

if exist(strcat('br01RelAct_',num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.BRIK'), 'file')
    delete(strcat('br01RelAct_',num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.HEAD'))
    delete(strcat('br01RelAct_',num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.BRIK'))
end
end