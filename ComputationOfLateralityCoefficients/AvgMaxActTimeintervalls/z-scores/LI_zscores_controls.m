
function for_all ()

% function created textfiles for Maxvalue in ROI, zum Ausrechnen des LI's


    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    LI_ALL_voxelvalues=zeros(10, 24)

    for i= 1:size(nameFolds)
   

    TimeInt = [.32, .6];
        
       [LI_ALL_voxelvalues]= kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt, i, LI_ALL_voxelvalues, 1:4, 9:12, 17:20)
       [LI_ALL_voxelvalues]= kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt, i , LI_ALL_voxelvalues, 5:8, 13:16, 21:24)

    end

    
    Path_LI_All_noise_abs=strcat('/home/kh/ShareWindows/data/', filesep, 'LI_controls_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_voxelvalues_noise_abs.mat')
    save (Path_LI_All_noise_abs, 'LI_ALL_voxelvalues')

end

function [LI_ALL_voxelvalues]=kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt, i, LI_ALL_voxelvalues, A, B, C)

% time dimension: to plot left and right activity over time

Path2oldROI = strcat( SubjectPath, filesep, 'SAM')
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
    
[ind_LeftAct]=find(LeftAct>2.326348)
zscores_LeftAct_p01=LeftAct(ind_LeftAct)

[ind_RightAct]=find(RightAct>2.326348)
zscores_RightAct_p01=RightAct(ind_RightAct)

LI.LI_Voxelvalue_p01=[sum(zscores_LeftAct_p01)-sum(zscores_RightAct_p01)]./[sum(zscores_LeftAct_p01)+sum(zscores_RightAct_p01)]



LI.LI_Voxelcount_p01=[length(zscores_LeftAct_p01)-length(zscores_RightAct_p01)]./[length(zscores_LeftAct_p01)+length(zscores_RightAct_p01)];
LI.Voxelcount_p01_LeftVox=length(zscores_LeftAct_p01);
LI.Voxelcount_p01_RightVox=length(zscores_RightAct_p01);



%% p>.001

[ind_LeftAct_p001]=find(LeftAct>3.090232)
zscores_LeftAct_p001=LeftAct(ind_LeftAct_p001)

[ind_RightAct_p001]=find(RightAct>3.090232)
zscores_RightAct_p001=RightAct(ind_RightAct_p001)

LI.LI_Voxelvalue_p001=(sum(zscores_LeftAct_p001)-sum(zscores_RightAct_p001))./(sum(zscores_LeftAct_p001)+sum(zscores_RightAct_p001));
LI.LI_Voxelcount_p001=(length(zscores_LeftAct_p001)-length(zscores_RightAct_p001))./(length(zscores_LeftAct_p001)+length(zscores_RightAct_p001))
LI.Voxelcount_p001_LeftVox=length(zscores_LeftAct_p001);
LI.Voxelcount_p001_RightVox=length(zscores_RightAct_p001);

LI_Path=strcat(SubjectPath, filesep, 'LI_', ROI, '_Voxelvalue_noise_abs_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_Path, 'LI') 


%%

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


function findextrema(SubjectPath, ROI, SubjectName)


PathMask_left = strcat ('/home/kh/ShareWindows/data/patients/', ROI, '+tlrc');
% eval(['!3dcalc -a ERF_noise_0.32-0.6s_Pat_01_13021km+tlrc -b ', PathMask_left, ' -exp a*b -prefix ',ROI, '_zscores'])
eval(['!3dExtrema -closure -volume -output ''Location_extrema', ROI,''' -mask_file ', PathMask_left, ' ERF_noise_0.32-0.6s_', SubjectName,'+tlrc ']) 

end

function [LI_All_noise]=collect_LI (SubjectPath, i, LI_All_noise )



load(strcat(SubjectPath, filesep,  'LI_Max_squared_noise_abs_Broca_0.32_0.6ms.mat'))
LI_All(1,1)=LI_Max_sqared;
clear LI_Max_sqared

load(strcat(SubjectPath, filesep, 'LI_squared_noise_abs_Broca_SumMax_0.32_0.6ms.mat'))
LI_All(1,2)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep, 'LI_squared_noise_abs__Broca_SumMax_0.32_0.47ms.mat'))
LI_All(1,3)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep,  'LI_squared_noise_abs__Broca_SumMax_0.4_0.6ms.mat'))
LI_All(1,4)=LI_SumMax;
clear LI_SumMax


load(strcat(SubjectPath, filesep,   'LI_Max_squared_noise_abs__Wernicke_0.32_0.6ms.mat'))
LI_All(1,5)=LI_Max_sqared;
clear LI_Max

load(strcat(SubjectPath, filesep, 'LI_squared_noise_abs__Wernicke_SumMax_0.32_0.6ms.mat'))
LI_All(1,6)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep, 'LI_squared_noise_abs__Wernicke_SumMax_0.32_0.47ms.mat'))
LI_All(1,7)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep,  'LI_squared_noise_abs__Wernicke_SumMax_0.4_0.6ms.mat'))
LI_All(1,8)=LI_SumMax;
clear LI_SumMax

Path_LI_All = strcat(SubjectPath, filesep, 'LI_All_Maxima_noise_abs_')
save(Path_LI_All, 'LI_All')

LI_All_noise_abs_(i,:)=LI_All

end



