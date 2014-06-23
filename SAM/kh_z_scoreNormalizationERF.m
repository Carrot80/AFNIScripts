
function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
   kh_z_transform( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .47)
%    kh_z_transform( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .31) 
%    kh_z_transform( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .2, .6) 
%    
%    
end


end

function kh_z_transform (SubjectPath, SubjectName, TimeBeg, TimeEnd)


    % (X-MEAN(X)) ./ STD(X)
    %  Vzscore = (V-mean(V)./std(V);

    Path = strcat(SubjectPath, filesep, 'SAM');
    cd (Path)
    
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