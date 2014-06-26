
function for_all ()

% function created textfiles for Maxvalue in ROI


    ControlsFolder = '/home/kh/data/controls_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    for i= 1:size(nameFolds)
        
    TimeInt = [.32, .55];
        
%         kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left', 'Broca_right_recalc', 'Broca', TimeInt)
%         kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left', 'Wernicke_right_recalc', 'Wernicke', TimeInt)
        kh_extractActROI(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'TempInfplusPols_left', 'TempInfplusPols_right_recalc', 'TempInfplusPols', TimeInt)

    end

end

function kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt)

% time dimension: to plot left and right activity over time

Path2oldROI = strcat( SubjectPath, filesep, 'SAM', filesep)
cd (Path2oldROI)

FileName = strcat('ERF+tlrc');
[V_ERF, Info_ERF] = BrikLoad (FileName);

PathMask_left = strcat (Path2oldROI, ROI_left, '+tlrc');
[Mask_left, Info_MASK_left] = BrikLoad (PathMask_left);


PathMask_right = strcat (Path2oldROI, ROI_right, '+tlrc');
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

SampleBeg = (TimeInt(1,1)*fs)+512;
SampleEnd = (TimeInt(1,2)*fs)+512;

Sum_Max_RightAct = Max_RightAct(1, SampleBeg:SampleEnd);

Sum_Max_LeftAct = Max_LeftAct(1, SampleBeg:SampleEnd);


Path_fig = strcat (SubjectPath, filesep, 'LI', filesep, ROI, '_extractedActivity');


saveas (gcf, Path_fig, 'fig')
print ('-dpng', Path_fig)



end
