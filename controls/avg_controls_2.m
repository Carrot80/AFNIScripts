
function forAll()


ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

 %% 

AVG = [];

for i= 1:size(nameFolds)
    
   [Info, AVG, TimeBeg, TimeEnd] = kh_get_avg( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .47, AVG, i)
   
   
end

    kh_avg_AllSubj(AVG, Info, TimeBeg, TimeEnd)
    

  
  
end




function [Info, AVG, TimeBeg, TimeEnd]=kh_get_avg (SubjectPath, SubjectName, TimeBeg, TimeEnd, AVG, i)


    Path = strcat(SubjectPath, filesep, 'SAM');
    cd (Path)
    
    FileName = strcat('br_z_transf_brain01ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc');   

    [V, Info] = BrikLoad (FileName);
    
    AVG{i}=V;
    

    %%
%     Vdif=V-Vlr;
%     Vdif(1:16,:,:)=0;
%     InfoNewTSOut = Info;
%     InfoNewTSOut.RootName = '';
%     InfoNewTSOut.BRICK_STATS = [];
%     InfoNewTSOut.BRICK_FLOAT_FACS = [];
%     InfoNewTSOut.IDCODE_STRING = '';
%     InfoNewTSOut.BRICK_TYPES=3*ones(1,1); % 1 short, 3 float.
%     InfoNewTSOut.DATASET_RANK(2)=1;
 
    
end

function kh_avg_AllSubj (AVG, Info, TimeBeg, TimeEnd)

AVG_all = (AVG{1}+AVG{2}+AVG{3}+AVG{4}+AVG{5}+AVG{6}+AVG{7}+AVG{8}+AVG{9}+AVG{10})./length(AVG);
   
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('AVG_controls_',num2str(TimeBeg),'_', num2str(TimeEnd), 's');
    OptTSOut.verbose = 1;
    OptTSOut.View = '+tlrc'
    %Vsymm=double(Vlr+V>0);
    
    cd /home/kh/ShareWindows/data/controls/AVGcontrols
    WriteBrik (AVG_all, Info, OptTSOut);
   

end