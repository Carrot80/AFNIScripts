
function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
TimeInt = [.32, .6];

for i= 5:size(nameFolds)
    
%    kh_z_transform( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .6)
%    kh_z_transform( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .31) 
%    kh_z_transform( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .6) 
   kh_convert2nifti( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, TimeInt)
%    
end


end

function kh_z_transform (SubjectPath, SubjectName, TimeBeg, TimeEnd)


    % (X-MEAN(X)) ./ STD(X)
    %  Vzscore = (V-mean(V)./std(V);

    Path = strcat(SubjectPath, filesep, 'TimeIntervalls');
    cd (Path)
    
    
    FileNameOld = strcat('ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc');

    disp(['!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ', FileNameOld, ' -prefix ', strcat('brain01', FileNameOld), ' -exp ', 'b*a'])
    eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileNameOld, ' -prefix ', strcat('brain01', FileNameOld), ' -exp ', 'b*a'])

    FileName = strcat('brain01ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc');   

    [V, Info] = BrikLoad (FileName);


    Vz=(V-mean(V(:)))/std(V(:));


    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('z_transf_', FileName);
    OptTSOut.verbose = 1;
    OptTSOut.View = '+tlrc'
    %Vsymm=double(Vlr+V>0);

    WriteBrik (Vz, Info, OptTSOut);

%     Vdif=V-Vlr;
%     Vdif(1:16,:,:)=0;
%     InfoNewTSOut = Info;
%     InfoNewTSOut.RootName = '';
%     InfoNewTSOut.BRICK_STATS = [];
%     InfoNewTSOut.BRICK_FLOAT_FACS = [];
%     InfoNewTSOut.IDCODE_STRING = '';
%     InfoNewTSOut.BRICK_TYPES=3*ones(1,1); % 1 short, 3 float.
%     InfoNewTSOut.DATASET_RANK(2)=1;


FileNameNew = OptTSOut.Prefix;   
   
eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileNameNew,  ' -prefix ', strcat('br_', FileNameNew),' -exp ' , 'b*a'])
    
    
end



function kh_convert2nifti(SubjectPath, SubjectName, TimeInt)


PathName = strcat(SubjectPath, filesep, 'TimeIntervalls');
cd (PathName)

PathERF = strcat('br_z_transf_brain01ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, '+tlrc');

% eval(['!@auto_tlrc -apar ', SubjectPath, 'orthoMNI_avg152T+tlrc',' -input ', PathERF, ' -dxyz 5']);

eval (['!3dcopy ', PathERF, ' ',strcat('br_z_transf_brain01ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, 'MNI.nii')])
% eval (['!3dcopy ', strcat('ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, '+tlrc'), ' ',strcat('ERF_', num2str(TimeInt(1,1)), '-', num2str(TimeInt(1,2)), 's_', SubjectName, '_MNI.nii')])


end