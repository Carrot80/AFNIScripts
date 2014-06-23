
function for_all ()

% function created textfiles for Maxvalue in ROI, zum Ausrechnen des LI's


    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];


    for i= 1:size(nameFolds)
   

    TimeInt = [.32, .6];
        
     kh_Afni2Nifti(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, TimeInt)
    
  

    end

 

end

function kh_Afni2Nifti (SubjectPath, SubjectName, TimeInt)

Path2oldROI = strcat( SubjectPath, filesep, 'SAM')
cd (Path2oldROI)



PathERF = strcat('br_z_transf_brain01ERF_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_', SubjectName, '+tlrc');


eval (['!3dcopy ', PathERF, ' ',strcat('br_z_transf_brain01ERF', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_', SubjectName, 'MNI.nii')])
eval (['!3dcopy ', 'orthoMNI_avg152T+tlrc', ' ',strcat('orthoMNI_avg152T.nii')])



end