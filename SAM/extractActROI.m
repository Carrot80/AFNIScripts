function for_all ()

% function created textfiles for Maxvalue in ROI


ControlsFolder = '/home/kh/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
TimeIntvls = [.2, .4; .3, .5; .4, .6; .5, .7; .32, .6; .25, .55; .25, .65 ];

for i= 1:size(nameFolds)
    
    kh_recalc(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, TimeIntvls, 'Broca_left', 'Broca_right_recalc')
    kh_recalc(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_right')
    kh_recalc(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, 'TempInfplusPols_right')
    
end

end

function kh_recalc(SubjectPath, SubjectName, TimeIntvls, ROI_left, ROI_right)

Path2oldROI = strcat( SubjectPath, filesep, 'SAM', filesep)
cd (Path2oldROI)

FileName = strcat('br_z_transf_brain01ERF_', num2str(TimeIntvls(5,1)), '-', num2str(TimeIntvls(5,2)), 's', '_', SubjectName, '+tlrc');

[V, Info] = BrikLoad (FileName);

PathMask_left = strcat (Path2oldROI, ROI_left, '+tlrc');

[Mask_left, Info_MASK_left] = BrikLoad (PathMask_left);

Left_Voxels = find(Mask_left==1);

leftAct=V(Left_Voxels);


PathMask_right = strcat (Path2oldROI, ROI_right, '+tlrc');

[Mask_right, Info_MASK_right] = BrikLoad (PathMask_right);

Right_Voxels = find(Mask_right==1);

rightAct=V(Right_Voxels);


median_rightAct = median(rightAct);
median_leftAct = median(leftAct);

LeftAboveThr=find(leftAct>median_leftAct) % 
RightAboveThr=find(rightAct>median_rightAct)

sum(rightAct(RightAboveThr))
sum(leftAct(LeftAboveThr))

mean(leftAct)
mean(rightAct)

boxplot(rightAct)
title('rightAct')
figure
boxplot(leftAct)  % evtl. Maske vergrößern? Boxplots zum Vergleich abspeichern
title('leftAct')

%% time dimension: to plot left and right activity over time


[V_ERF, Info_ERF] = BrikLoad ('ERF+tlrc');


length(V_ERF)

[~,Li]=ismember(LRpairs(:,1),ga.label); % ga ) grandaverage with individually kept trials, LRpairs is file from Yuval
[~,Ri]=ismember(LRpairs(:,2),ga.label);
gaLR=ga;
gaLR.individual=zeros(size(gaLR.individual)); % die mittleren Sensoren behalten Nullen
gaLR.individual(:,Li,:)=abs(ga.individual(:,Li,:))-abs(ga.individual(:,Ri,:));
gaLR.individual(:,Ri,:)=abs(ga.individual(:,Ri,:))-abs(ga.individual(:,Li,:));

aliceTtest0(gaLR, 0.4, 1);
vgenTtest0(gaLR, 0.4, 1,0.05,ga); % Input Nr. 3 = p-value
vgenTtest0(gaLR, 0.17, 1,0.05,ga);
vgenTtest0(gaLR, 0.3, 1,0.05,ga);


% compute RMS
rmsL=squeeze(sqrt(mean(ga.individual(:,Li,:).^2,2)));
rmsR=squeeze(sqrt(mean(ga.individual(:,Ri,:).^2,2)));

t=0.4;[~,p] = ttest(rmsL(:,nearest(ga.time,t)),rmsR(:,nearest(ga.time,t)));

figure;plot(ga.time,mean(rmsL),'r')
hold on
plot(ga.time,mean(rmsR))
legend('L','R')

t1=0.32;s1=nearest(ga.time,t1);
t2=0.47;s2=nearest(ga.time,t2);
areaL=trapz(ga.time(s1:s2),rmsL(:,s1:s2)');
areaR=trapz(ga.time(s1:s2),rmsR(:,s1:s2)');
[~,p]=ttest(areaL,areaR)

end

