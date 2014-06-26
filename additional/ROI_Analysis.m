function for_all ()

% function created textfiles for Maxvalue in ROI


ControlsFolder = '/home/kh/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];


for i= 1:size(nameFolds)
    
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .4, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Broca+tlrc', 'Broca_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .4, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Broca+tlrc', 'Broca_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .4, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .4, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .4, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Pols_Inf_Temp+tlrc', 'TempInfplusPols_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .4, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Pols_Inf_Temp+tlrc', 'TempInfplusPols_right_recalc' )
   
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .4, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Broca+tlrc', 'Broca_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .4, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Broca+tlrc', 'Broca_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .4, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .4, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .4, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Pols_Inf_Temp+tlrc', 'TempInfplusPols_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .4, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Pols_Inf_Temp+tlrc', 'TempInfplusPols_right_recalc' )
   
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .3, .5, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Broca+tlrc', 'Broca_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .3, .5, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Broca+tlrc', 'Broca_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .3, .5, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .3, .5, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .3, .5, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Pols_Inf_Temp+tlrc', 'TempInfplusPols_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .3, .5, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Pols_Inf_Temp+tlrc', 'TempInfplusPols_right_recalc' )
   
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .5, .7, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Broca+tlrc', 'Broca_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .5, .7, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Broca+tlrc', 'Broca_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .5, .7, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .5, .7, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .5, .7, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Pols_Inf_Temp+tlrc', 'TempInfplusPols_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .5, .7, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Pols_Inf_Temp+tlrc', 'TempInfplusPols_right_recalc' )
   
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Broca+tlrc', 'Broca_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Broca+tlrc', 'Broca_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Pols_Inf_Temp+tlrc', 'TempInfplusPols_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .6, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Pols_Inf_Temp+tlrc', 'TempInfplusPols_right_recalc' )
 
   
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .55, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Broca+tlrc', 'Broca_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .55, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Broca+tlrc', 'Broca_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .55, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .55, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .55, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Pols_Inf_Temp+tlrc', 'TempInfplusPols_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .55, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Pols_Inf_Temp+tlrc', 'TempInfplusPols_right_recalc' )
   
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .65, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Broca+tlrc', 'Broca_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .65, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Broca+tlrc', 'Broca_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .65, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .65, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Ang_Supr_Mid_Sup_Temp+tlrc', 'Wernicke_right_recalc' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .65, '/home/kh/data/ROIsForAfni/AAL/AAL_left_Pols_Inf_Temp+tlrc', 'TempInfplusPols_left' )
   GetExtremValue( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .65, '/home/kh/data/ROIsForAfni/AAL/AAL_rightfliped_Pols_Inf_Temp+tlrc', 'TempInfplusPols_right_recalc' )
      
end


end

function GetExtremValue( SubjectPath, SubjectName, TimeBeg, TimeEnd, ROI_Mask, NewMaskName)

% take extrema from z-scores:
PathSAM = strcat (SubjectPath, filesep, 'SAM');
cd (PathSAM)

FileName= strcat ('br_z_transf_brain01ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's_', SubjectName, '+tlrc');   


% eval(['!3dresample -master ', FileName, ' -prefix Broca_left -inset AAL_left_Broca+tlrc' ])
eval(['!3dresample -master ', FileName, ' -prefix ', NewMaskName, ' -inset ', ROI_Mask ])  
% use ROIMaskrecalc.m here, as right ROIS don't consist of ones and zeros:
%  !3dcalc -a Broca_right+tlrc -exp 'ispositive(a)' -prefix Broca_right_recalc 

Textfile = strcat (NewMaskName, '_', num2str(TimeBeg), '_', num2str(TimeEnd), 's_Max.txt');

% closure: also takes extrema at the border of mask:

eval(['!3dExtrema -volume -sep_dist 30 -closure -mask_file ', NewMaskName, '+tlrc ', FileName, ' > ', Textfile ]) 






end
    