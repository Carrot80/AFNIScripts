
function for_all ()

% function created textfiles for Maxvalue in ROI, zum Ausrechnen des LI's


    ControlsFolder = '/media/truecrypt7/kirsten_thesis/data/controls/';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    Table=[];

    for i= 1:size(nameFolds)
   
      
   [Table]=findextrema(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, Table, i, 'Verbgeneration')  
    
    end

    Path=strcat('/home/kh/ShareWindows/Extrema_and_Locations_fMRI_Verbgeneration_t-values_patients.mat')
    save(Path, 'Table')
    
end


function [Table]=findextrema(SubjectPath, SubjectName, Table, i, Task)



Path=strcat(SubjectPath, filesep, 'fMRI/statistics/', Task);

cd(Path)

!3dcopy wspmT_0001.hdr wspmT_0001

% Maske auf T-Verteilung anpassen:
    % PathMask_complete = strcat ('/home/kh/ShareWindows/Masks_fMRI/Brainmask+tlrc');
    % PathMask_left = strcat ('/home/kh/ShareWindows/Masks_fMRI/Left_Brainmask+tlrc');
    % PathMask_right = strcat ('/home/kh/ShareWindows/Masks_fMRI/Right_Brainmask+tlrc');
    %  eval(['!3dresample -master wspmT_0001+tlrc', ' -prefix ', ' Brainmask_complete_fMRI', ' -inset ', PathMask_complete ]) 
    %  eval(['!3dresample -master wspmT_0001+tlrc', ' -prefix ', ' Brainmask_left_fMRI', ' -inset ', PathMask_left ]) 
    %  eval(['!3dresample -master wspmT_0001+tlrc', ' -prefix ', ' Brainmask_right_fMRI', ' -inset ', PathMask_right ]) 

 
PathMask_complete_fMRI = strcat ('/home/kh/ShareWindows/Masks_fMRI/Brainmask_complete_fMRI+tlrc');
PathMask_left_fMRI = strcat ('/home/kh/ShareWindows/Masks_fMRI/Brainmask_left_fMRI+tlrc');
PathMask_right_fMRI = strcat ('/home/kh/ShareWindows/Masks_fMRI/Brainmask_right_fMRI+tlrc');
 
    eval(['!3dExtrema -volume -mask_file ', PathMask_complete_fMRI, ' wspmT_0001+tlrc > Extrema_Mask_complete.txt'])
        eval(['!3dExtrema -volume -mask_file ', PathMask_left_fMRI, ' wspmT_0001+tlrc > Extrema_left_Brain.txt'])
            eval(['!3dExtrema -volume -mask_file ', PathMask_right_fMRI, ' wspmT_0001+tlrc > Extrema_right_Brain.txt'])

[Table, Mask]=kh_Whereami (SubjectPath, SubjectName, Table, i, 'Mask_complete', 1)
[Table, Mask]=kh_Whereami (SubjectPath, SubjectName, Table, i, 'left_Brain', 2)
[Table, Mask]=kh_Whereami (SubjectPath, SubjectName, Table, i, 'right_Brain', 3)

end


function [Table, Mask]=kh_Whereami (SubjectPath, SubjectName, Table, i, Mask, A)

[newData1]=importfile (strcat('Extrema_', Mask, '.txt'))

if 0==isstruct(newData1) %kh, Patient18 hat no extreme values
    Table.(Mask).Location{i,1} = 'subject has no extreme Value';
    Table.(Mask).Coord(i,1:4)  = [NaN NaN NaN NaN]
    return
end

eval(['!whereami ', num2str(newData1.data(1,3)),' ', num2str(newData1.data(1,4)),' ', num2str(newData1.data(1,5)), ' > Whereami_' Mask,'.txt' ])

fid = fopen(strcat('Whereami_',Mask, '.txt'));
Whereami=textscan(fid, '%s', 'delimiter', sprintf('\f'));
fclose(fid)

Number=strfind( Whereami{1,1}, 'Atlas CA_ML_18_MNIA: Macro Labels (N27)')
for j=1:length(Number)
    X(j)=~isempty(Number{j})
end

  X=double(X)
[row]=find(X==1)

Table.(Mask).Location{i,1} = Whereami{1,1}{row+1}
Table.(Mask).Coord(i,1:4)  = newData1.data(1,2:5)

end


function [newData1]=importfile(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 27-May-2014 16:09:09

DELIMITER = ' ';
HEADERLINES = 10;

% Import the file
newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);

% Create new variables in the base workspace from those fields.
if 0==isstruct(newData1) %kh, Patient18 hat no extreme values 
    return
end

vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

end