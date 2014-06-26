function patients ()
   
PatientFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';
ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

for_all (PatientFolder, 'Patients')
for_all (ControlsFolder, 'Controls')

end


function for_all (Folder, group)

% function created textfiles for Maxvalue in ROI, zum Ausrechnen des LI's

    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    LI_All_noise=zeros(size(nameFolds), 8);
    for i= 1:size(nameFolds)
   

    TimeInt = [.32, .6];
        
        kh_extractActROI(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt, group)
        kh_extractActROI(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt, group)



      [LI_All_noise]= collect_LI (strcat(Folder, filesep, nameFolds{i,1}), i, LI_All_noise, nameFolds{i,1})


    end

    Path_LI_All_noise=strcat('/home/kh/ShareWindows/data/', group, filesep, 'LI_All_noise_abs')
    save (Path_LI_All_noise, 'LI_All_noise')

end



function kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt, group)

if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg') 
    return
end
% time dimension: to plot left and right activity over time
if 1==strcmp(group, 'Controls')
    Path = strcat( SubjectPath, filesep, 'SAM')
    cd (Path)
    FileName = strcat('ERF+tlrc');
else
    cd (SubjectPath)
    FileName = strcat('ERF_noise+tlrc');
end

% Path2MNI = strcat( SubjectPath, filesep, 'keptTrials', filesep, 'orthoMNI_avg152T+tlrc')
% NewFileName = 'ERF_noise+orig';
% eval(['!@auto_tlrc -apar ', strcat(SubjectPath, filesep, 'keptTrials', filesep, 'orthoMNI_avg152T+tlrc'), ' -input ', NewFileName,' -dxyz 5']) % 

[V_ERF, Info_ERF] = BrikLoad (FileName);

PathMask_left = strcat ('/home/kh/ShareWindows/data/patients/', ROI_left, '+tlrc');
[Mask_left, Info_MASK_left] = BrikLoad (PathMask_left); 

PathMask_right = strcat ('/home/kh/ShareWindows/data/patients/', ROI_right, '+tlrc');
[Mask_right, Info_MASK_right] = BrikLoad (PathMask_right);

for i= 1:length(V_ERF)
    LeftAct(:,:,:,i) = Mask_left.*V_ERF(:,:,:,i);
end

for i= 1:length(V_ERF)
    RightAct(:,:,:,i) = Mask_right.*V_ERF(:,:,:,i);
end

%%
for i=1:length(RightAct)
    Max_RightAct(i) = max(max(max(RightAct(:,:,:,i))));
end

for i=1:length(LeftAct)
    Max_LeftAct(i) = max(max(max(LeftAct(:,:,:,i))));
end

fs = 1017.25;
Time = 1:length(LeftAct);

TimeSec =Time./fs;
%  SampleOffset = nearest(TimeSec, 0.500);
TimeSec_Offset = TimeSec-.5;

SampleBeg=nearest(TimeSec_Offset,TimeInt(1,1))
SampleEnd=nearest(TimeSec_Offset,TimeInt(1,2))

Max_RightAct_TimeInt = Max_RightAct(SampleBeg:SampleEnd);

Max_LeftAct_TimeInt = Max_LeftAct(SampleBeg:SampleEnd);

Max_right=max(Max_RightAct_TimeInt);
Max_left=max(Max_LeftAct_TimeInt);

% LI=(Max_left-Max_right)/(Max_left+Max_right)
LI_squared=(Max_left^2-Max_right^2)/(Max_left^2+Max_right^2);
LI_sqrt=sqrt(abs(LI_squared));
if LI_squared <0
    LI_Max_sqared=LI_sqrt*(-1);
else LI_Max_sqared=LI_sqrt;
end

LI_Max_squared_Path=strcat(SubjectPath, filesep, 'LI_Max_squared_noise_', ROI, '_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat');
save (LI_Max_squared_Path, 'LI_Max_sqared') 

LI_Max=(Max_left-Max_right)/(Max_left+Max_right);

LI_Max_Path=strcat(SubjectPath, filesep, 'LI_Max_noise_', ROI, '_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat');
save (LI_Max_Path, 'LI_Max') 


%summe funktioniert nciht, da nur bilaterale LIs rauskommen:
Sum_Max_LeftAct_TimeInt     = sum(Max_LeftAct_TimeInt);
Sum_Max_RightAct_TimeInt    = sum(Max_RightAct_TimeInt);
LI_SumMax_squared           = (Sum_Max_LeftAct_TimeInt^2-Sum_Max_RightAct_TimeInt^2)/(Sum_Max_LeftAct_TimeInt^2+Sum_Max_RightAct_TimeInt^2);

LI_SumMax_sqrt=sqrt(abs(LI_SumMax_squared))
if LI_SumMax_squared <0
    LI_SumMax=LI_SumMax_sqrt*(-1)
else LI_SumMax=LI_SumMax_sqrt
end

LI_SumMax_Path=strcat(SubjectPath, filesep, 'LI_squared_noise_', ROI, '_SumMax_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_SumMax_Path, 'LI_SumMax') 

LI_SumMax_noise           = (Sum_Max_LeftAct_TimeInt-Sum_Max_RightAct_TimeInt)./(Sum_Max_LeftAct_TimeInt+Sum_Max_RightAct_TimeInt);
LI_SumMax_Path=strcat(SubjectPath, filesep, 'LI_noise_', ROI, '_SumMax_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_SumMax_Path, 'LI_SumMax_noise') 

end


function [LI_All_noise]=collect_LI (SubjectPath, i, LI_All_noise, SubjectName )

if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg') 

    return
end

load(strcat(SubjectPath, filesep,  'LI_Max_squared_noise_Broca_0.32_0.6ms.mat'))
LI_All(1,1)=LI_Max_sqared;
clear LI_Max_sqared

load(strcat(SubjectPath, filesep, 'LI_squared_noise_Broca_SumMax_0.32_0.6ms.mat'))
LI_All(1,2)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep, 'LI_squared_noise_Broca_SumMax_0.32_0.47ms.mat'))
LI_All(1,3)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep,  'LI_squared_noise_Broca_SumMax_0.4_0.6ms.mat'))
LI_All(1,4)=LI_SumMax;
clear LI_SumMax


load(strcat(SubjectPath, filesep,   'LI_Max_squared_noise_Wernicke_0.32_0.6ms.mat'))
LI_All(1,5)=LI_Max_sqared;
clear LI_Max

load(strcat(SubjectPath, filesep, 'LI_squared_noise_Wernicke_SumMax_0.32_0.6ms.mat'))
LI_All(1,6)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep, 'LI_squared_noise_Wernicke_SumMax_0.32_0.47ms.mat'))
LI_All(1,7)=LI_SumMax;
clear LI_SumMax

load(strcat(SubjectPath, filesep,  'LI_squared_noise_Wernicke_SumMax_0.4_0.6ms.mat'))
LI_All(1,8)=LI_SumMax;
clear LI_SumMax

Path_LI_All = strcat(SubjectPath, filesep, 'LI_All_Maxima_noise')
save(Path_LI_All, 'LI_All')

LI_All_noise(i,:)=LI_All

end



