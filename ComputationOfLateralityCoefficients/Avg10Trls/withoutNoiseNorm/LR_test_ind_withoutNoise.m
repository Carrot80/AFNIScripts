 
function for_all ()

% Script written by Yuval:

ControlsFolder = '/home/kh/ShareWindows/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];


TimeIntervall (nameFolds, ControlsFolder)


end


function TimeIntervall (nameFolds, ControlsFolder)

VlrAll= [];
Vall = [];

TimeInt = [.32, .6]; % TimInt unten noch anpassen

for i= 1:size(nameFolds)
    SubjectPath = strcat(ControlsFolder, filesep, nameFolds{i,1});
    SubjectName = nameFolds{i};
    
    [VlrAll, Vall] = get_V (SubjectPath, SubjectName, VlrAll, Vall, TimeInt);
      
    UtestLR (SubjectPath, SubjectName, VlrAll, Vall, TimeInt)
    
end   
     

end

%%


function [VlrAll, Vall] = get_V (SubjectPath, SubjectName, VlrAll, Vall, TimeInt)

% SubjectToAnalyse = strcat('/media/truecrypt7/Linux_data/KeepTrialsNoNoiseEstimation/', 'zzz_ca');
% if 1 == strcmp (SubjectPath, SubjectToAnalyse)
% return
% end

  PathName = strcat(SubjectPath, filesep, 'UTestNoNoise');
  cd (PathName)

for i= 1:8


    FileName = strcat('ERF_avgTrials_', num2str(i), '+orig');

    [V, Info] = BrikLoad (FileName);
    
    % Sum Samples in Time Interval
    fs = 1017.25;

    offset_samples = 509;
      
    
    V_1_1000ms = V(:,:,:, offset_samples:length(V)); 
    
    time_samples=1:size(V_1_1000ms, 4);
    time_sec=time_samples./fs;
    
    
    sample_int_Beg=size(find(time_sec<=TimeInt(1,1)));
    sample_int_End=size(find(time_sec<=TimeInt(1,2)));
    
    % sum forth dimension:
    V_SumTime=sum(V_1_1000ms(:,:,:, sample_int_Beg(2):sample_int_End(2)),4);
    
    clear V V_1_1000ms       
    
    Vlr=flipdim(V_SumTime,2);
    if i==1
        Vall=V_SumTime;
        VlrAll=Vlr;
    else
        Vall(:,:,:,i)=V_SumTime;
        VlrAll(:,:,:,i)=Vlr;
    end
    
    clear V_SumTime
    
end

end

%%

function UtestLR (SubjectPath, SubjectName, VlrAll, Vall, TimeInt)
 
% SubjectToAnalyse = strcat('/media/truecrypt7/Linux_data/KeepTrialsNoNoiseEstimation/', 'zzz_ca');
% 
% if 1 == strcmp (SubjectPath, SubjectToAnalyse)
% return
% end

    U=zeros(size(Vall,1),size(Vall,2),size(Vall,3));
    
    for i=1:size(Vall,1)
        for j=1:size(Vall,2)
            for k=1:size(Vall,3)
                if Vall(i,j,k)>0
                    p = ranksum(squeeze(Vall(i,j,k,:)),squeeze(VlrAll(i,j,k,:)));
                    u=1-p;
                    dif=mean(squeeze(Vall(i,j,k,:)))-mean(squeeze(VlrAll(i,j,k,:)));
                    if dif>0
                        R=1;% R is 1 when current is larger than other side
                    else
                        R=-1;
                    end
                    U(i,j,k)=R*u;
                end
            end
        end
        i
    end

    [~, Info] = BrikLoad (strcat('/home/kh/ShareWindows/data/controls_SAM/', SubjectName, filesep, 'SAM/', 'ERF_0.4-0.6s_', SubjectName, '+orig')); % Info would be 

      PathName = strcat(SubjectPath, filesep, 'UTestNoNoise');
      cd (PathName)
    
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('Utest', '_', 'LR', '_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's');
    OptTSOut.verbose = 1;
    OptTSOut.View = '+orig' ;
    %Vsymm=double(Vlr+V>0);
    WriteBrik (U, Info, OptTSOut);
    
end


