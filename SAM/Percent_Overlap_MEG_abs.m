function patients ()
   
PatientFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';
ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

for_all (PatientFolder, 'patients')
for_all (ControlsFolder, 'controls')

end


function for_all (Folder, group)


    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    Comparison=[];
    LI_fMRI_informed_MEG=[];
    TimeInt = [.32, .47];

    for i= 1:size(nameFolds)
    
   fMRI_informed_MEG=[];  
   [Comparison, fMRI_informed_MEG]=kh_comp(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, Comparison, 'Broca_left_dil', 1, fMRI_informed_MEG)  
   [Comparison, fMRI_informed_MEG]=kh_comp(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, Comparison, 'Broca_right_dil', 2, fMRI_informed_MEG)  
   [Comparison, fMRI_informed_MEG]=kh_comp(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, Comparison, 'Wernicke_left_dil', 3, fMRI_informed_MEG)  
   [Comparison, fMRI_informed_MEG]=kh_comp(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, Comparison, 'Wernicke_right_dil', 4, fMRI_informed_MEG)  
   [fMRI_informed_MEG, LI_fMRI_informed_MEG]=LI_informed(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, fMRI_informed_MEG, LI_fMRI_informed_MEG)
   
    end

  PathComparison=strcat('/home/kh/ShareWindows/Results/Comparison_fMRI_MEG', filesep, 'Comparison_Overlap_MEG_fMRI_', group,'_', num2str(TimeInt(1)), '_', num2str(TimeInt(2)), 's_abs.mat')
  save(PathComparison, 'Comparison')
  PathLI_fMRIInformedMEG=strcat('/home/kh/ShareWindows/Results/Comparison_fMRI_MEG', filesep, 'LI_fMRI_Informed_MEG_', group,'_', num2str(TimeInt(1)), '_', num2str(TimeInt(2)), 's_abs.mat')
  save(PathLI_fMRIInformedMEG, 'LI_fMRI_informed_MEG')
    
end


function [Comparison, fMRI_informed_MEG]=kh_comp(SubjectPath, SubjectName, i, TimeInt, group, Comparison, ROI, A, fMRI_informed_MEG)

if 1==strcmp(SubjectName, 'Pat_02_13008rh_1') || 1==strcmp(SubjectName, 'Pat_02_13008rh_2') || 1==strcmp(SubjectName,'Pat_03_13014bg_1') || 1==strcmp(SubjectName,'Pat_03_13014bg_2')
    return
end

PathMask = strcat ('/home/kh/ShareWindows/data/patients/', ROI, '+tlrc');

if 1==strcmp(group, 'patients')
    Path2Subj=strcat(SubjectPath, filesep, 'TimeIntervalls');
else
    Path2Subj=strcat(SubjectPath, filesep, 'SAM');
end

cd(Path2Subj)


if 1==strcmp (SubjectName,'Pat_02_13008rh') ||  1==strcmp (SubjectName,'Pat_03_13014bg')
    PathERF = strcat('BothRuns_br_z_transf_brain01ERF_noise_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's', '+tlrc');
elseif 1==strcmp(group, 'Patients')
    PathERF = strcat('br_z_transf_brain01ERF_noise_abs_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_', SubjectName, '+tlrc');
else
    PathERF = strcat('br_z_transf_brain01ERF_noise_abs_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_', SubjectName, '+tlrc');
end

[V_ERF, Info_ERF] = BrikLoad (PathERF);

[V_Mask, Info_Mask] = BrikLoad (PathMask); 

Act_MEG = V_Mask.*V_ERF;

Path_fMRI=strcat('/media/truecrypt7/kirsten_thesis/data/', group, filesep, SubjectName,'/fMRI/statistics/Verbgeneration');

[V_fMRI, Info_fMRI] = BrikLoad (strcat(Path_fMRI, filesep, 'wspmT_0001_downsampled_z_transf+tlrc'));

Act_fMRI=V_Mask.*V_fMRI;

Act_MEG_ROI=squeeze(Act_MEG(:));
Act_fMRI_ROI=squeeze(Act_fMRI(:));
[ind_MEG]=find(Act_MEG_ROI==0);
Act_MEG_ROI(ind_MEG)=[];

[ind_fMRI]=find(Act_fMRI_ROI==0);
Act_fMRI_ROI(ind_fMRI)=[];

[ind_Act_MEG_p05]=find(Act_MEG_ROI>1.644853);
Act_MEG_p05=Act_MEG_ROI(ind_Act_MEG_p05);
Act_MEG_p05_perc=length(Act_MEG_p05)/length(Act_MEG_ROI);

[ind_Act_fMRI_p05]=find(Act_fMRI_ROI>1.644853);
Act_fMRI_p05=Act_fMRI_ROI(ind_Act_fMRI_p05);
Act_fMRI_p05_perc=length(Act_fMRI_p05)/length(Act_fMRI_ROI);

if 1==length(ind_Act_fMRI_p05)<length(ind_Act_MEG_p05) 
    
   ind_Act_fMRI_p05_samelength= zeros(length(ind_Act_MEG_p05), 1);
   ind_Act_fMRI_p05_samelength(1:length(ind_Act_fMRI_p05))=ind_Act_fMRI_p05;
end     
    
   

[ind_ismember]=ismember(ind_Act_MEG_p05,ind_Act_fMRI_p05);
sum(ind_ismember);
overlap_fMRI_MEG=sum(ind_ismember)/length(ind_Act_MEG_p05);
overlap_MEG_fMRI=sum(ind_ismember)/length(Act_fMRI_p05);
overlap_both=sum(ind_ismember)/(length(ind_Act_fMRI_p05)+length(ind_Act_MEG_p05)-sum(ind_ismember));




%3d:
[ind_Act_MEG_p05_3D]=find(Act_MEG>1.644853);
[ind_Act_fMRI_p05_3D]=find(Act_fMRI>1.644853);
[ind_Act_MEG_p05_3D_ismember]=ismember(ind_Act_MEG_p05_3D,ind_Act_fMRI_p05_3D)
[pos]=(ind_Act_MEG_p05_3D(ind_Act_MEG_p05_3D_ismember==1))

fMRI_informed_MEG.(ROI).indices_3D=pos;
rows_ind=find(ind_ismember==1)
fMRI_informed_MEG.(ROI).Activation=Act_MEG_ROI(rows_ind);


% fmri informed MEG: calculate LIs

Comparison.active_voxels05_fMRI(i,A)            =length(ind_Act_fMRI_p05);
Comparison.percOfROIactive_voxels05_fMRI(i,A)   =length(ind_Act_fMRI_p05)/length(Act_fMRI_ROI)
Comparison.active_voxels05_MEG(i,A)             =length(ind_Act_MEG_p05);
Comparison.percOfROIactive_voxels05_MEG(i,A)    =length(ind_Act_MEG_p05)/length(Act_MEG_ROI)
Comparison.Overlap.fMRIoverlapsMEG(i,A)         =overlap_fMRI_MEG;
Comparison.Overlap.MEGoverlapsfMRI(i,A)         =overlap_MEG_fMRI;
Comparison.Overlap.Overlap_both(i,A)            =overlap_both;

end



function  [fMRI_informed_MEG, LI_fMRI_informed_MEG]=LI_informed(SubjectPath, SubjectName, i, fMRI_informed_MEG, LI_fMRI_informed_MEG)

if 1==strcmp(SubjectName, 'Pat_02_13008rh_1') || 1==strcmp(SubjectName, 'Pat_02_13008rh_2') || 1==strcmp(SubjectName,'Pat_03_13014bg_1') || 1==strcmp(SubjectName,'Pat_03_13014bg_2')
    return
end

LI_fMRI_informed_MEG.Broca.Voxelvalue_p05(i,1)=(sum(fMRI_informed_MEG.Broca_left_dil.Activation)-sum(fMRI_informed_MEG.Broca_right_dil.Activation))./(sum(fMRI_informed_MEG.Broca_left_dil.Activation)+sum(fMRI_informed_MEG.Broca_right_dil.Activation))
LI_fMRI_informed_MEG.Broca.Voxelcount_p05(i,1)=(length(fMRI_informed_MEG.Broca_left_dil.Activation)-length(fMRI_informed_MEG.Broca_right_dil.Activation))./(length(fMRI_informed_MEG.Broca_left_dil.Activation)+length(fMRI_informed_MEG.Broca_right_dil.Activation))
LI_fMRI_informed_MEG.Wernicke.Voxelvalue_p05(i,1)=(sum(fMRI_informed_MEG.Wernicke_left_dil.Activation)-sum(fMRI_informed_MEG.Wernicke_right_dil.Activation))./(sum(fMRI_informed_MEG.Wernicke_left_dil.Activation)+sum(fMRI_informed_MEG.Wernicke_right_dil.Activation))
LI_fMRI_informed_MEG.Wernicke.Voxelcount_p05(i,1)=(length(fMRI_informed_MEG.Wernicke_left_dil.Activation)-length(fMRI_informed_MEG.Wernicke_right_dil.Activation))./(length(fMRI_informed_MEG.Wernicke_left_dil.Activation)+length(fMRI_informed_MEG.Wernicke_right_dil.Activation))



% außerdem: Möglichkeit, neue Maske zu erstellen mit Überlappungsbereich



end


