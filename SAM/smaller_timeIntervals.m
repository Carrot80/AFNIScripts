function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
   kh_SAM_TimeInt( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .3, .6)
   kh_SAM_TimeInt( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .55) 
   kh_SAM_TimeInt( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .25, .65) 
  
   
end


end

function kh_SAM_TimeInt (SubjectPath, SubjectName, TimeBeg, TimeEnd)

Path2SAM = strcat (SubjectPath, filesep, 'SAM');
cd(Path2SAM);
Path2VS = strcat (SubjectPath, filesep, 'SAM', filesep, 'Workspace_SAM.mat');
load (Path2VS) ;

fs = 1017.25;
size_vs=size(vs);

offset_samples = 509;
vs_1_1000ms = vs(:,offset_samples:size_vs(2));
size_vs_1_1000ms=size(vs_1_1000ms);
time_samples=1:size_vs_1_1000ms(2);
time_sec=time_samples./fs;
% I could have used function nearest, too.

figure
plot(time_sec,max(abs(vs_1_1000ms)));
axis tight;
Title=strcat(num2str(TimeBeg), '_', num2str(TimeEnd), 's_', SubjectName);
title(Title)
print (Title); 
close all

sample_int_Beg=size(find(time_sec<=TimeBeg));
sample_int_End=size(find(time_sec<=TimeEnd));

vs_IntOfIn=vs_1_1000ms(:,sample_int_Beg(2):sample_int_End(2));
sum_vs_IntOfIn = abs(sum(vs_IntOfIn')); % Achtung falsch!


% Save it to load it in afni
cfg=[];
cfg.step=5;
cfg.boxSize=[-120 120 -90 90 -20 150];
str_timeInt= strcat('ERF', '_', num2str(TimeBeg), '-', num2str(TimeEnd), 's_', SubjectName);
cfg.prefix = str_timeInt; % change prefix
% cfg.torig=-500;   %  comment if you want to sum up activity of specific time intervall
% cfg.TR=1/1.01725; % comment if you want to sum up activity of specific time intervall
% VS2Brik(cfg,vs);
% max(max(vs))
VS2Brik(cfg,1e+13*(sum_vs_IntOfIn')); % =>creates ERF+orig.Brik+Head 

NewFileName = strcat(str_timeInt,'+orig');
%  disp(['!@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ', NewFileName,' -dxyz 5']) % 
eval(['!@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ', NewFileName,' -dxyz 5']) % 

kh_reduceERF2Brain (SubjectPath, SubjectName, TimeBeg, TimeEnd)
kh_z_transform (SubjectPath, SubjectName, TimeBeg, TimeEnd)

end


function kh_reduceERF2Brain (SubjectPath, SubjectName, TimeBeg, TimeEnd)

% if 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ca') || 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ht')
%     return
% end

SAMPath = strcat(SubjectPath, filesep, 'SAM');
cd (SAMPath)

FileNameOld = strcat('ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc');

disp(['!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ', FileNameOld, ' -prefix ', strcat('brain01', FileNameOld), ' -exp ', 'b*a'])
eval(['!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ', FileNameOld, ' -prefix ', strcat('brain01', FileNameOld), ' -exp ', 'b*a'])


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
   
eval(['!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ', FileNameNew,  ' -prefix ', strcat('br_', FileNameNew),' -exp ' , 'b*a'])
    
    
end
