
function for_all ()

% function created textfiles for Maxvalue in ROI


    ControlsFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    for i= 8:size(nameFolds)
        
    TimeInt = [.32, .6];
        
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt)
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt)
%     

    end

end

function kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt)

% time dimension: to plot left and right activity over time

Path2oldROI = strcat( SubjectPath, filesep)
cd (Path2oldROI)

Path2MNI = strcat( SubjectPath, filesep, 'keptTrials', filesep, 'orthoMNI_avg152T+tlrc')

eval(['!@auto_tlrc -apar ',  Path2MNI, ' -input ERF+orig -dxyz 5'])

FileName = strcat('ERF+tlrc');
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

figure
plot(TimeSec_Offset, Max_LeftAct, 'r')
hold on
plot(TimeSec_Offset, Max_RightAct)
legend('left', 'right')
title(strcat(ROI, '_', SubjectName));
axis tight
set(gca, 'FontSize', 14)
Path_fig = strcat (SubjectPath, filesep, ROI, '_extractedActivity');
saveas (gcf, Path_fig, 'fig')
print ('-dpng', Path_fig)
close all


% SampleBeg = (TimeInt(1,1)*fs)+512;
% SampleEnd = (TimeInt(1,2)*fs)+512;
% 
% Sum_Max_RightAct = Max_RightAct(1, SampleBeg:SampleEnd);
% 
% Sum_Max_LeftAct = Max_LeftAct(1, SampleBeg:SampleEnd);








end
