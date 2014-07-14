
function for_all ()

% for AVG Controls

    ControlsFolder = '/home/kh/ShareWindows/data/controls/AVGcontrols';
    TimeInt = [.32, .47];
        
        kh_extractActROI(ControlsFolder, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt)
        kh_extractActROI(ControlsFolder, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt)
end

function kh_extractActROI (ControlsFolder, ROI_left, ROI_right, ROIName, TimeInt)

cd (ControlsFolder)

Path2UValues = strcat (ControlsFolder, filesep, 'Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's+tlrc' );

[V_UValues, Info_UValues] = BrikLoad (Path2UValues);


PathMask_left = strcat ('/home/kh/ShareWindows/data/patients', filesep, ROI_left, '+tlrc');
[Mask_left, Info_MASK_left] = BrikLoad (PathMask_left);


PathMask_right = strcat ('/home/kh/ShareWindows/data/patients', filesep, ROI_right, '+tlrc');
[Mask_right, Info_MASK_right] = BrikLoad (PathMask_right);


Left_Voxels = find(Mask_left==1);

leftAct=V_UValues(Left_Voxels);

[signLeftVoxels]=find(leftAct>=.95);


Right_Voxels = find(Mask_right==1);

rightAct=V_UValues(Right_Voxels);

[signRightVoxels]=find(rightAct>=.95);


% calculate LI:

 LI.signRightVoxels  = length(signRightVoxels);
 LI.signLeftVoxels = length(signLeftVoxels);
 LI.SizeROI = size(Left_Voxels, 1);
 LI.relActLeft = length(signLeftVoxels)./size(Left_Voxels, 1);
 LI.relActRight = length(signRightVoxels)./size(Right_Voxels, 1);
 LI.Max_percchange = (length(signLeftVoxels)-length(signRightVoxels))./length(signRightVoxels);
 LI.Classic = (length(signLeftVoxels)-length(signRightVoxels))./(length(signLeftVoxels)+length(signRightVoxels));
 if LI.Classic >= .2
      LI.Lateralization = 'left';
 elseif LI.Classic <= -.2 
    LI.Lateralization = 'right';
 elseif LI.Classic > -.2 && LI.Classic < +.2
 LI.Lateralization = 'bilateral';
 end
 

    
 PathFile = strcat (ControlsFolder, filesep, 'LI_', ROIName, 'dil_SumOfSignVoxels_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)),  '_s.mat' );
 save (PathFile, 'LI')

end
