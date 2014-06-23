
function for_all ()

% function created textfiles for Maxvalue in ROI


    ControlsFolder = '/home/kh/data/controls_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    for i= 1:size(nameFolds)
        
    TimeInt = [.32, .6];
        
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, '/home/kh/data/ROIsForAfni/AAL/L_Broca_3voxdil3Drecalc+tlrc', '/home/kh/data/ROIsForAfni/AAL/R_Broca_3voxdil3Drecalc+tlrc', 'Broca_left_dil', 'Broca_right_dil', TimeInt)
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, '/home/kh/data/ROIsForAfni/AAL/L_AngSuprTempSup_3voxdil3Drecalc+tlrc', '/home/kh/data/ROIsForAfni/AAL/R_AngSuprTempSup_3voxdil3Drecalc+tlrc', 'Wernicke_left_dil', 'Wernicke_right_dil', TimeInt)
%         kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'TempInfplusPols_left', 'TempInfplusPols_right_recalc', 'TempInfplusPols', TimeInt)

    end

end

function kh_extractActROI (SubjectPath, SubjectName, ROI_leftOld, ROI_rightOld, ROI_leftNEW, ROI_rightNEW, TimeInt)

% SubjectToAnalyse = '/home/kh/data/controls_SAM/zzz_ms';
% if 1 == strcmp (SubjectPath, SubjectToAnalyse)
%     return
% end

PathUTest = strcat (SubjectPath, filesep, 'UTest', filesep);
cd (PathUTest)

Path2UValues = strcat (SubjectPath, filesep, 'UTest', filesep, 'Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 'sMNI+tlrc' );
eval(['!3dresample -master ', Path2UValues, ' -prefix ', ROI_leftNEW, ' -inset ', ROI_leftOld ])  
eval(['!3dresample -master ', Path2UValues, ' -prefix ', ROI_rightNEW, ' -inset ', ROI_rightOld ]) 


eval(['!3dresample -master ', 'z_transf_brain01ERF_noise_0.32-0.6s_Pat_01_13021km+tlrc', ' -prefix ', 'Left_Brainmask_MEG_1mm', ' -inset ', 'Left_Brainmask_fMRI+tlrc' ]) 





end
