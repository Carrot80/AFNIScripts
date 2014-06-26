% berechnet LIs für AVG je 10 Trials; funktioniert noch nicht

function LR_TTST_AVG10()
    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';
    PatientFolder  = '/home/kh/ShareWindows/data/patients/patients_SAM';
    TimeInt = [.32, .6];
    for_all( ControlsFolder, 'Controls', TimeInt )
%     for_all( PatientFolder,  'Patients', TimeInt )
end 


function for_all (Folder, group, TimeInt)
    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    TimeIntervall (nameFolds, Folder, group, TimeInt)
end


function TimeIntervall (nameFolds, ControlsFolder, group, TimeInt)

    VlrAll= [];
    Vall = [];

    for i= 1:size(nameFolds)
        SubjectPath = strcat(ControlsFolder, filesep, nameFolds{i,1});
        SubjectName = nameFolds{i};
        kh_SAM(SubjectPath, SubjectName, group)
        [VlrAll, Vall] = get_V (SubjectPath, SubjectName, group, TimeInt, VlrAll, Vall);
        TTestLR (SubjectPath, SubjectName, group, TimeInt, VlrAll, Vall)
    end
end

%%

function kh_SAM(SubjectPath, SubjectName, group)

    if 1 == strcmp (SubjectName, 'Pat_03_13014bg') || 1 == strcmp (SubjectName, 'Pat_02_13008rh')
        return
    end

     % ToDo kann später entfernt werden:
     for j=1:10
         if exist(strcat(SubjectPath, filesep, 'keptTrials', filesep, 'ERF_avgTrials_', num2str(j),'+orig.BRIK'))
             eval (['! rm ', strcat(SubjectPath, filesep, 'keptTrials', filesep, 'ERF_avgTrials_', num2str(j), '+orig.BRIK')])
         end
         if exist(strcat(SubjectPath, filesep, 'keptTrials', filesep, 'ERF_avgTrials_', num2str(j),'+orig.HEAD'))
             eval (['! rm ', strcat(SubjectPath, filesep, 'keptTrials', filesep, 'ERF_avgTrials_', num2str(j), '+orig.HEAD')])
         end
     end

    NewFolder = strcat(SubjectPath, filesep, 'keptTrials');
    if ~exist(NewFolder, 'dir')
        mkdir(NewFolder)
    end

    switch group
        case 'Controls'
            load(strcat(SubjectPath, filesep, 'SAM', filesep, 'Workspace_SAM.mat'));
            save(strcat(SubjectPath, filesep, 'SAM', filesep, 'Workspace_SAM.mat'), 'ActIndex', 'ActWgts', 'SAMHeader', 'avgBL') % ToDo: kann nach einem Run entfernt werden 

        otherwise 
            PathSAM = strcat(SubjectPath, filesep, 'SAM', filesep);
            cd (PathSAM) % change matlab dir to Path of weights (SAMdir)
            [SAMHeader, ActIndex, ActWgts] = readWeights('M400,1-50Hz,VGa.wts');
    end        
            
    % load avg:
    if exist (strcat(SubjectPath, filesep, 'avgTrials_', SubjectName, '.mat'));
        movefile(strcat(SubjectPath, filesep, 'avgTrials_', SubjectName, '.mat'), strcat(SubjectPath, filesep, 'keptTrials'))
    end
    load(strcat(SubjectPath, filesep, 'keptTrials', filesep,'avgTrials_', SubjectName, '.mat' ))
    
    for i=1:length(avgTrials.trial(:,1,1))

        vs = [];
        ns = [];
        vs=ActWgts*squeeze(avgTrials.trial(i,:,:));
        ns=mean(abs(ActWgts),2);
        vs=vs./repmat(ns,1,size(vs,2));

        cd(NewFolder)

        cfg         = [];
        cfg.step    = 5;
        cfg.boxSize = [-120 120 -90 90 -20 150];
        cfg.prefix  = strcat('ERF_avgTrials_', num2str(i)); % change prefix
        cfg.torig   = -500;   %  comment if you want to sum up activity of specific time intervall
        cfg.TR      = 1/1.01725; % comment if you want to sum up activity of specific time intervall
        VS2Brik(cfg,1e+13*abs(vs)); % =>creates ERF+orig.Brik+Head

    %     MNIFile = strcat (SubjectPath, filesep, 'SAM', filesep, 'orthoMNI_avg152T+tlrc');
%     eval ([ '!@auto_tlrc -apar ', MNIFile, ' -input ', 'ERF_avgTrials_', num2str(i), '+orig -dxyz 5'  ]);  
     
    end
end



function [VlrAll, Vall] = get_V (SubjectPath, SubjectName, group, TimeInt, VlrAll, Vall)

    if 1 == strcmp (SubjectName, 'Pat_03_13014bg') || 1 == strcmp (SubjectName, 'Pat_02_13008rh')
        return
    end

    Path = strcat(SubjectPath, filesep, 'keptTrials');

for i= 1:10

    FileName = strcat(Path, filesep, 'ERF_avgTrials_', num2str(i), '+orig');

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
    V_SumTime=sum(V_1_1000ms(:,:,:, sample_int_Beg(2):sample_int_End(2)),4); % stimmt so, da nur positive Werte summiert werden (siehe kh_SAM_Beamforming_keeptrials => abs(vs))
    
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

function TTestLR (SubjectPath, SubjectName,  group, TimeInt, VlrAll, Vall)
 
    if 1 == strcmp (SubjectName, 'Pat_03_13014bg') || 1 == strcmp (SubjectName, 'Pat_02_13008rh')
        return
    end

    T=zeros(size(Vall,1),size(Vall,2),size(Vall,3));
    
    for i=1:size(Vall,1)
        for j=1:size(Vall,2)
            for k=1:size(Vall,3)
                if Vall(i,j,k)>0
                    [~, p] = signrank(squeeze(Vall(i,j,k,:)),squeeze(VlrAll(i,j,k,:)));
                    u=1-p;
                    dif=mean(squeeze(Vall(i,j,k,:)))-mean(squeeze(VlrAll(i,j,k,:)));
                    if dif>0
                        R=1;% R is 1 when current is larger than other side
                    else
                        R=-1;
                    end
                    T(i,j,k)=R*u;
                end
            end
        end
    end

    [~, Info] = BrikLoad (strcat(SubjectPath, filesep, 'TimeIntervalls', filesep, 'ERF_0.32-0.6s_', SubjectName, '+orig')); % Info-File wird benötigt  
    
    Path = strcat(SubjectPath, filesep, 'keptTrials');    
    cd (Path) 
    
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('Signranktest', '_', 'LR', '_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's');
    OptTSOut.verbose = 1;
    OptTSOut.View = '+orig' ;
    %Vsymm=double(Vlr+V>0);
    WriteBrik (T, Info, OptTSOut);
    
    % remove files to clear space:
    for j=1:10
        delete (strcat(SubjectPath, filesep, 'keptTrials', filesep, 'ERF_avgTrials_', num2str(j), '+orig.BRIK'))
        delete (strcat(SubjectPath, filesep, 'keptTrials', filesep, 'ERF_avgTrials_', num2str(j), '+orig.HEAD'))
    end
    
end


