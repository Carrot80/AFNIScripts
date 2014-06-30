
function for_all ()

% column output 1 to 4: LI_Voxvaluep05Broca LI_Voxcountp05Broca LeftVoxelsp05Broca RightVoxelsp05Broca
% column output 5 to 8: LI_Voxvaluep05Wernicke LI_Voxcountp05Wernicke LeftVoxelsp05Wernicke RightVoxelsp05Wernicke
% function created textfiles for Maxvalue in ROI, zum Ausrechnen des LI's

    ControlsFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';
    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    LI_ALL_voxelvalues=zeros(24, 24)

    for i= 1:size(nameFolds)
        TimeInt = [.32, .6];
        [LI_ALL_voxelvalues]= kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt, i, LI_ALL_voxelvalues, 1:4, 9:12, 17:20)
        [LI_ALL_voxelvalues]= kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt, i , LI_ALL_voxelvalues, 5:8, 13:16, 21:24)
    end
    Path_LI_All_noise_abs=strcat('/home/kh/ShareWindows/data/', filesep, 'patients', 'LI_All_voxelvalues_noise_abs')
    save (Path_LI_All_noise_abs, 'LI_ALL_voxelvalues')

end

function [LI_ALL_voxelvalues]=kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt, i, LI_ALL_voxelvalues, A, B, C)

% time dimension: to plot left and right activity over time
if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg') 

    return
end
Path2oldROI = strcat( SubjectPath, filesep, 'TimeIntervalls')
cd (Path2oldROI)

PathERF = strcat('br_z_transf_brain01ERF_noise_abs_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_', SubjectName, '+tlrc');

[V_ERF, Info_ERF] = BrikLoad (PathERF);

PathMask_left = strcat ('/home/kh/ShareWindows/data/patients/', ROI_left, '+tlrc');
[Mask_left, Info_MASK_left] = BrikLoad (PathMask_left); 

PathMask_right = strcat ('/home/kh/ShareWindows/data/patients/', ROI_right, '+tlrc');
[Mask_right, Info_MASK_right] = BrikLoad (PathMask_right);

    LeftAct = Mask_left.*V_ERF;
    RightAct = Mask_right.*V_ERF;
%% p>.01
    
[ind_LeftAct]=find(LeftAct>2.33)
zscores_LeftAct=LeftAct(ind_LeftAct)

[ind_RightAct]=find(RightAct>2.33)
zscores_RightAct=RightAct(ind_RightAct)

LI.LI_Voxelvalue_p01=(sum(zscores_LeftAct)-sum(zscores_RightAct))./(sum(zscores_LeftAct)+sum(zscores_RightAct))
LI.LI_Voxelcount_p01=(length(zscores_LeftAct)-length(zscores_RightAct))./(length(zscores_LeftAct)+length(zscores_RightAct))
LI.Voxelcount_p01_LeftVox=length(zscores_LeftAct);
LI.Voxelcount_p01_RightVox=length(zscores_RightAct);

%% p>.001

[ind_LeftAct_p001]=find(LeftAct>3.090232)
zscores_LeftAct_p001=LeftAct(ind_LeftAct_p001)

[ind_RightAct_p001]=find(RightAct>3.090232)
zscores_RightAct_p001=RightAct(ind_RightAct_p001)

LI.LI_Voxelvalue_p001=(sum(zscores_LeftAct_p001)-sum(zscores_RightAct_p001))./(sum(zscores_LeftAct_p001)+sum(zscores_RightAct_p001))
LI.LI_Voxelcount_p001=(length(zscores_LeftAct_p001)-length(zscores_RightAct_p001))./(length(zscores_LeftAct_p001)+length(zscores_RightAct_p001))
LI.Voxelcount_p001_LeftVox=length(zscores_LeftAct_p001);
LI.Voxelcount_p001_RightVox=length(zscores_RightAct_p001);

LI_Path=strcat(SubjectPath, filesep, 'LI_', ROI, '_Voxelvalue_noise_abs_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_Path, 'LI') 

%% p>.05

[ind_LeftAct_p05]=find(LeftAct>1.644853)
zscores_LeftAct_p05=LeftAct(ind_LeftAct_p05)

[ind_RightAct_p05]=find(RightAct>1.644853)
zscores_RightAct_p05=RightAct(ind_RightAct_p05)

LI.LI_Voxelvalue_p05=(sum(zscores_LeftAct_p05)-sum(zscores_RightAct_p05))./(sum(zscores_LeftAct_p05)+sum(zscores_RightAct_p05));
LI.LI_Voxelcount_p05=(length(zscores_LeftAct_p05)-length(zscores_RightAct_p05))./(length(zscores_LeftAct_p05)+length(zscores_RightAct_p05));
LI.Voxelcount_p05_LeftVox=length(zscores_LeftAct_p05);
LI.Voxelcount_p05_RightVox=length(zscores_RightAct_p05);

LI_Path=strcat(SubjectPath, filesep, 'LI_', ROI, '_Voxelvalue_noise_abs_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_Path, 'LI') 

LI_ALL_voxelvalues(i,A)= [LI.LI_Voxelvalue_p05 LI.LI_Voxelcount_p05 LI.Voxelcount_p05_LeftVox LI.Voxelcount_p05_RightVox ]
LI_ALL_voxelvalues(i,B)= [LI.LI_Voxelvalue_p01 LI.LI_Voxelcount_p01 LI.Voxelcount_p01_LeftVox LI.Voxelcount_p01_RightVox ]
LI_ALL_voxelvalues(i,C)= [LI.LI_Voxelvalue_p001 LI.LI_Voxelcount_p001 LI.Voxelcount_p001_LeftVox LI.Voxelcount_p001_RightVox ]
end
