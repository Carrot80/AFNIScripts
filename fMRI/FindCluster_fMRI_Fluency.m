
function patients ()
   
PatientFolder = '/media/truecrypt7/kirsten_thesis/data/patients/';
ControlsFolder = '/media/truecrypt7/kirsten_thesis/data/controls/';

% for_all (PatientFolder, 'Patients', 'Verbgeneration')
% for_all (ControlsFolder, 'Controls', 'Verbgeneration')
% for_all (PatientFolder, 'Patients', 'Fluency')
% for_all (ControlsFolder, 'Controls', 'Fluency')
calcDist()

end

function for_all (Folder, group, Task)

    
    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    Table=[];

    for i= 1:size(nameFolds)
   
   [Table]=findclust(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, Table, i, Task)  
   
    end

    Path=strcat('/home/kh/ShareWindows/data/Cluster_and_Locations_fMRI_t-values_', group, '_', Task, '_ROIs.mat')
    save(Path, 'Table')

    
end


function [Table]=findclust(SubjectPath, SubjectName, Table, i, Task)



Path=strcat(SubjectPath, filesep, 'fMRI/statistics/', Task);

cd(Path)


% Maske auf T-Verteilung anpassen:
%     PathMask = strcat ('/home/kh/ShareWindows/data/patients/ROIs_complete_MEG+tlrc');
%     eval(['!3dresample -master wspmT_0001+tlrc', ' -prefix ', ' ROIs_complete_fMRI', ' -inset ', PathMask ])

% !3dcalc -a Broca_left_dil_fMRI+tlrc -b Broca_right_dil_fMRI+tlrc -exp 'a+b' -prefix Broca_both_hemi_fMRI
 
PathROI_complete = strcat ('/home/kh/ShareWindows/data/patients/ROIs_complete_fMRI+tlrc');
PathBroca = strcat ('/home/kh/ShareWindows/Masks_fMRI/Broca_both_hemi_fMRI+tlrc');
PathWernicke = strcat ('/home/kh/ShareWindows/Masks_fMRI/Wernicke_both_hemi_fMRI+tlrc');

if ~exist ('Broca_wspmT_0001+tlrc.BRIK', 'file')
    [V, Info]=BrikLoad('wspmT_0001+tlrc');
    [V_Broca, Info_Broca]=BrikLoad(PathBroca);
    V_TValues_Broca=V.*V_Broca;
    
    OptTSOut.Prefix = 'Broca_wspmT_0001'
    OptTSOut.View = '+tlrc'
    WriteBrik (V_TValues_Broca, Info, OptTSOut);
end

if 1==strcmp(Task, 'Verbgeneration')
    eval(['!3dclust -1noneg -1thresh 1.6526 3 1000 Broca_wspmT_0001+tlrc > Clust_Broca.txt'])

else
    eval(['!3dclust -1noneg -1thresh 1.658 3 1000 Broca_wspmT_0001+tlrc > Clust_Broca.txt'])
     
end


%%

if ~exist ('Wernicke_wspmT_0001+tlrc.BRIK', 'file')
    [V_Wernicke, Info_Wernicke]=BrikLoad(PathWernicke);
    V_TValues_Wernicke=V.*V_Wernicke;
    
    OptTSOut.Prefix = 'Wernicke_wspmT_0001'
    OptTSOut.View = '+tlrc'
    WriteBrik (V_TValues_Wernicke, Info, OptTSOut);
end

if 1==strcmp(Task, 'Verbgeneration')
    eval(['!3dclust -1noneg -1thresh 1.6526 3 1000 Wernicke_wspmT_0001+tlrc > Clust_Wernicke.txt'])
else
    eval(['!3dclust -1noneg -1thresh 1.658 3 1000 Wernicke_wspmT_0001+tlrc > Clust_Wernicke.txt'])
end


 
[Table, Mask]=kh_Whereami (SubjectPath, SubjectName, Table, i, 'Broca', 1)
[Table, Mask]=kh_Whereami (SubjectPath, SubjectName, Table, i, 'Wernicke', 2)


end


function [Table, Mask]=kh_Whereami (SubjectPath, SubjectName, Table, i, Mask, A)

[newData1]=importfile2 (strcat('Clust_', Mask, '.txt'));

if 0==isstruct(newData1) %kh, Patient18 hat no extreme values
    Table.(Mask).Location{i,1} = 'subject has no extreme Value';
    Table.(Mask).Coord(i,1:4)  = [NaN NaN NaN NaN];
    return
end

% Max=max(newData1.data(:,13))
% [row_Max, col]=find(newData1.data==Max)
[row_Max, col]=find(newData1.data==max(newData1.data(:,13)))

eval(['!whereami ', num2str(newData1.data(row_Max,14)),' ', num2str(newData1.data(row_Max,15)),' ', num2str(newData1.data(row_Max,16)), ' > Clust_Whereami_' Mask,'.txt' ])

fid = fopen(strcat('Clust_Whereami_',Mask, '.txt'));
Whereami=textscan(fid, '%s', 'delimiter', sprintf('\f'));
fclose(fid)

Number=strfind( Whereami{1,1}, 'Atlas CA_ML_18_MNIA: Macro Labels (N27)')
for j=1:length(Number)
    X(j)=~isempty(Number{j})
end

  X=double(X)
[row_1]=find(X==1)

if 1==isempty(row_1)
    Table.(Mask).Location{i,1} = 'no CA_ML_18_MNIA Labels';
    Table.(Mask).Coord(i,1:16) = newData1.data(row_Max,1:16);
    return
end

Table.(Mask).Location{i,1} = Whereami{1,1}{row_1+1};
Table.(Mask).Coord(i,1:16)  = newData1.data(row_Max,1:16);

end


function calcDist()

% calculates Clusterdistance
% Broca Patients

Table_Pat_Fluency=load('/home/kh/ShareWindows/data/Cluster_and_Locations_fMRI_t-values_Patients_Fluency_ROIs.mat')
Table_Pat_VG=load('/home/kh/ShareWindows/data/Cluster_and_Locations_fMRI_t-values_Patients_Verbgeneration_ROIs.mat')

d_Broca_patients=zeros(size(Table_Pat_Fluency.Table.Broca.Coord,1),1)

for j=1:length(Table_Pat_Fluency.Table.Broca.Coord)
    
    d_Broca_patients(j,1)=sqrt((Table_Pat_Fluency.Table.Broca.Coord(j,14)-Table_Pat_VG.Table.Broca.Coord(j,14))^2+(Table_Pat_Fluency.Table.Broca.Coord(j,15)-Table_Pat_VG.Table.Broca.Coord(j,15))^2+(Table_Pat_Fluency.Table.Broca.Coord(j,16)-Table_Pat_VG.Table.Broca.Coord(j,16))^2)
    
end

Path_Pat='/home/kh/ShareWindows/data/Clusterdistance_Pat_Fluency_VG_Broca_FMRI.mat';
save(Path_Pat, 'd_Broca_patients')

% Wernicke:

d_Wernicke_patients=zeros(size(Table_Pat_Fluency.Table.Wernicke.Coord,1),1)

for j=1:length(Table_Pat_Fluency.Table.Wernicke.Coord)
    
    d_Wernicke_patients(j,1)=sqrt((Table_Pat_Fluency.Table.Wernicke.Coord(j,14)-Table_Pat_VG.Table.Wernicke.Coord(j,14))^2+(Table_Pat_Fluency.Table.Wernicke.Coord(j,15)-Table_Pat_VG.Table.Wernicke.Coord(j,15))^2+(Table_Pat_Fluency.Table.Wernicke.Coord(j,16)-Table_Pat_VG.Table.Wernicke.Coord(j,16))^2)
    
end

Path_Pat='/home/kh/ShareWindows/data/Clusterdistance_Pat_Wernicke_FMRI.mat';
save(Path_Pat, 'd_Wernicke_patients')

%% controls

Table_Contr_Fluency=load('/home/kh/ShareWindows/data/Cluster_and_Locations_fMRI_t-values_Controls_Fluency_ROIs.mat')
Table_Contr_VG=load('/home/kh/ShareWindows/data/Cluster_and_Locations_fMRI_t-values_Controls_Verbgeneration_ROIs.mat')

d_Broca_contr=zeros(size(Table_Contr_Fluency.Table.Broca.Coord,1),1)

for j=1:10
    
    d_Broca_contr(j,1)=sqrt((Table_Contr_Fluency.Table.Broca.Coord(j,14)-Table_Contr_VG.Table.Broca.Coord(j,14))^2+(Table_Contr_Fluency.Table.Broca.Coord(j,15)-Table_Contr_VG.Table.Broca.Coord(j,15))^2+(Table_Contr_Fluency.Table.Broca.Coord(j,16)-Table_Contr_VG.Table.Broca.Coord(j,16))^2)
    
end


Path_Contrl='/home/kh/ShareWindows/data/Clusterdistance_Contr_Fluency_VG_Broca_FMRI.mat';
save(Path_Contrl, 'd_Broca_contr')

% Wernicke:

d_Wernicke_contr=zeros(size(Table_Contr_Fluency.Table.Wernicke.Coord,1),1)

for j=1:10
    
    d_Wernicke_contr(j,1)=sqrt((Table_Contr_Fluency.Table.Wernicke.Coord(j,14)-Table_Contr_VG.Table.Wernicke.Coord(j,14))^2+(Table_Contr_Fluency.Table.Wernicke.Coord(j,15)-Table_Contr_VG.Table.Wernicke.Coord(j,15))^2+(Table_Contr_Fluency.Table.Wernicke.Coord(j,16)-Table_Contr_VG.Table.Wernicke.Coord(j,16))^2)
    
end


Path_Contrl='/home/kh/ShareWindows/data/Clusterdistance_Contr_Wernicke_FMRI.mat';
save(Path_Contrl, 'd_Wernicke_contr')

end


function [newData1]=importfile2(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 03-Jun-2014 14:21:13

% Import the file
newData1 = importdata(fileToRead1);

if 0==isstruct(newData1) %kh, Patient18 hat no extreme values 
    return
end

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

end

