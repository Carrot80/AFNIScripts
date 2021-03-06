
function for_all ()

% function created textfiles for Maxvalue in ROI, zum Ausrechnen des LI's


    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    Table=[];

    for i= 1:size(nameFolds)
   

    TimeInt = [.32, .6];
      
   [Table]=findextrema(strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, Table, i)  
    
    end

    Path=strcat('/home/kh/ShareWindows/data', filesep, 'Extrema_and_Locations_controls.mat')
    save(Path, 'Table')
    
end


function [Table]=findextrema(SubjectPath, SubjectName, Table, i)

if 1==strcmp(SubjectName,'Pat_02_13008rh_1') || 1==strcmp(SubjectName,'Pat_02_13008rh_2') || 1==strcmp(SubjectName,'Pat_03_13014bg_1') || 1==strcmp(SubjectName,'Pat_03_13014bg_2')

    return
end

Path=strcat(SubjectPath, filesep, 'SAM');

cd(Path)

PathMask_complete = strcat ('/home/kh/ShareWindows/data/Brainmask+tlrc');
PathMask_left = strcat ('/home/kh/ShareWindows/data/Left_Brainmask+tlrc');
PathMask_right = strcat ('/home/kh/ShareWindows/data/Right_Brainmask+tlrc');

% eval(['!3dcalc -a ERF_noise_0.32-0.6s_Pat_01_13021km+tlrc -b ', PathMask_left, ' -exp a*b -prefix Left_Brainmask'])
% eval(['!3dcalc -a ERF_noise_0.32-0.6s_Pat_01_13021km+tlrc -b ', PathMask_right, ' -exp a*b -prefix Right_Brainmask'])
 eval(['!3dresample -master ', 'ERF_noise_0.32-0.6s_Pat_01_13021km+tlrc', ' -prefix ', 'Brainmask', ' -inset ', PathMask_complete ]) 
% eval(['!3dresample -master ', 'ERF_noise_0.32-0.6s_Pat_01_13021km+tlrc', ' -prefix ', 'Left_Brainmask', ' -inset ', PathMask_left ]) 


if 1==strcmp (SubjectName,'Pat_02_13008rh') ||  1==strcmp (SubjectName,'Pat_03_13014bg')
    eval(['!3dExtrema -volume -mask_file ', PathMask_complete, ' BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s+tlrc > Extrema_Mask_complete.txt'])
    eval(['!3dExtrema -volume -mask_file ', PathMask_left, ' BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s+tlrc > Extrema_left_Brain.txt'])
    eval(['!3dExtrema -volume -mask_file ', PathMask_right, ' BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s+tlrc > Extrema_right_Brain.txt'])
else
    
    eval(['!3dExtrema -volume -mask_file ', PathMask_complete, ' br_z_transf_brain01ERF_0.32-0.6s_', SubjectName,'+tlrc > Extrema_Mask_complete.txt'])
        eval(['!3dExtrema -volume -mask_file ', PathMask_left, ' br_z_transf_brain01ERF_0.32-0.6s_', SubjectName,'+tlrc > Extrema_left_Brain.txt'])
            eval(['!3dExtrema -volume -mask_file ', PathMask_right, ' br_z_transf_brain01ERF_0.32-0.6s_', SubjectName,'+tlrc > Extrema_right_Brain.txt'])
end

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