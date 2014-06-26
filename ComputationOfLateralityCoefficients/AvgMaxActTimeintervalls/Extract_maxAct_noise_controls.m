
function for_all ()

% function created textfiles for Maxvalue in ROI, zum Ausrechnen des LI's


    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    LI_All_noise=zeros(24, 8)
    for i= 3:4 %:size(nameFolds)
   

    TimeInt = [.32, .6; .32, .47; .4, .6];
        
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt(1,:))
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt(1,:))
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt(2,:))
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt(2,:))
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt(3,:))
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt(3,:))
  
%       [LI_All_noise]= collect_LI (strcat(ControlsFolder, filesep, nameFolds{i,1}), i, LI_All_noise)


    end

%     Path_LI_All_noise=strcat('/home/kh/ShareWindows/data/patients', filesep, 'LI_All_noise')
%     save (Path_LI_All_noise, 'LI_All_noise')

end

function kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt)

% time dimension: to plot left and right activity over time

Path2oldROI = strcat( SubjectPath, filesep, 'SAM')
cd (Path2oldROI)

Path2MNI = strcat( SubjectPath, filesep, 'SAM', filesep, 'orthoMNI_avg152T+tlrc')

NewFileName = 'ERF_noise+orig';
eval(['!@auto_tlrc -apar ', strcat(SubjectPath, filesep, 'keptTrials', filesep, 'orthoMNI_avg152T+tlrc'), ' -input ', NewFileName,' -dxyz 5']) % 

FileName = strcat('ERF_noise+tlrc');
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


SampleBeg = (TimeInt(1,1)*fs)+512;
SampleEnd = (TimeInt(1,2)*fs)+512;

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

LI_SumMax_noise           = (Sum_Max_LeftAct_TimeInt-Sum_Max_RightAct_TimeInt)/(Sum_Max_LeftAct_TimeInt+Sum_Max_RightAct_TimeInt);
LI_SumMax_Path=strcat(SubjectPath, filesep, 'LI_noise_', ROI, '_SumMax_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_SumMax_Path, 'LI_SumMax_noise') 




end


function [LI_All_noise]=collect_LI (SubjectPath, i, LI_All_noise )



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



