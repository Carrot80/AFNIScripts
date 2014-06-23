function LR_ttest_dep()
    controlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';
    PatientFolder  = '/home/kh/ShareWindows/data/patients/patients_SAM';
    TimeInt = [.32, .6, 320, 600; .32 .47, 320, 470; .4, .6, 400, 600 ];
    for_all( controlsFolder, 'controls', TimeInt )
    for_all( PatientFolder,  'patients', TimeInt )
end


function for_all (Folder, group, TimeInt)
    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    Main (nameFolds, Folder, group, TimeInt)
end


function Main (nameFolds, Folder, group, TimeInt)

    VlrAll = [];
    Vall   = [];
    LI_All = [];
    for i = 1:size(nameFolds)
        SubjectPath = strcat(Folder, filesep, nameFolds{i,1});
        SubjectName = nameFolds{i};

        %       [avg] = AVG_CleanData (SubjectPath, SubjectName, group)
        %       kh_SAM(SubjectPath, SubjectName, group, avg)

        %       for ztime=1:size(TimeInt,1)
        %           [VlrAll, Vall] = get_V (SubjectPath, SubjectName, VlrAll, Vall, TimeInt(ztime,:));
        %           TTestLR (SubjectPath, SubjectName, VlrAll, Vall, TimeInt(ztime,:))
        %       end

        %       deleteFiles(SubjectPath, SubjectName) % belongs to TTestLR

        for ztime=1:size(TimeInt,1)
            %           kh_TTest_normalize (SubjectPath, SubjectName, TimeInt(ztime,:))
            %           kh_extractActROI (SubjectPath, SubjectName, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt(ztime,:))
            %           kh_extractActROI (SubjectPath, SubjectName, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt(ztime,:))
          
%             [LI_All] = collect_LI (SubjectPath, SubjectName, strcat(num2str(TimeInt(ztime,1)),'_', num2str(TimeInt(ztime,2)),'_s' ), LI_All, strcat('Int',num2str(TimeInt(ztime, 3)),'To',num2str(TimeInt(ztime, 4)), 'ms'))
        end
    end
%     Path = strcat( '/home/kh/ShareWindows/data/', group, filesep, 'LIs_signVoxels_singleTrials.mat'); 
%     save (Path, 'LI_All')
    for ztime=1:size(TimeInt,1)
        collect_LI_excel (group, strcat('Int',num2str(TimeInt(ztime, 3)),'To',num2str(TimeInt(ztime, 4)), 'ms'))
    end
end


function [avg]=AVG_CleanData (SubjectPath, SubjectName, group)
    % Hauptordner ist für Patienten mit 2 Runs uninteressant
    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg')
        avg=[];
        return
    end
    
    % falls bereits Dateien existieren, muss dieser Schritt nicht mehr durchgeführt werden 
    if exist(strcat(SubjectPath, filesep, 'TTest', filesep, 'ERF_1000ms_Trial_1+orig.BRIK'))
        avg=[];
        return 
    end
    
    switch SubjectName
        case 'Pat_03_13014bg_1'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanDataSNR.mat');
        case 'Pat_03_13014bg_2'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanDataSNR_best.mat');
        case  'Pat_07_13033gc'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanData_nobadTrls.mat');
        case  'Pat_08_13026pj'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanData_rejcomp.mat');
        case  'Pat_11_13030rs'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanData_rejcomp.mat');
        case  'Pat_14_13039sg'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanDataSNR2.mat');
        case  'Pat_17_13060ec'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanData_rejectcomp2.mat');
        case  'Pat_19_13055eg'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanData_rejcomp.mat');
        case  'Pat_21_13056hz'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanData_rejcomp108Trials.mat');
        case  'Pat_22_13059oc'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanData_rejcomp2.mat');
        case  'Pat_24_13067pj'
            File_CleanData= strcat (SubjectPath, filesep, 'TTest', filesep, 'CleanDataSNR1_5.mat');
        otherwise
            File_CleanData = strcat(SubjectPath, filesep, 'TTest', filesep, 'CleanData.mat');
    end

    % alte Zwischenergebnis löschen. Kann nach einen run entfernt werden (ToDo)
    if exist(strcat(SubjectPath, filesep, 'TTest', filesep, 'AvgKeepTrials.mat'), 'file')
        delete(strcat(SubjectPath, filesep, 'TTest', filesep, 'AvgKeepTrials.mat'))
    end    

    load(File_CleanData)
    CleanDataBL = correctBL(CleanData,[-0.32 -.02]);
    cfg = [];
    cfg.keeptrials='yes'; 
    
    % hier gehts weiter!
    avg = ft_timelockanalysis(cfg, CleanDataBL); % eigentlich nicht notwendig, wäre auch möglich, trials aus CleanData zu nehmen.

end

%%


function kh_SAM(SubjectPath, SubjectName, group, avg)

    % Hauptordner ist für Patienten mit 2 Runs uninteressant
    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg')
        return
    end

    % falls bereits Dateien existieren, muss dieser Schritt nicht mehr durchgeführt werden: 
    if exist(strcat(SubjectPath, filesep, 'TTest', filesep, 'ERF_1000ms_Trial_1+orig.BRIK'))
        return 
    end
    
    switch group
        case 'controls'
            load(strcat(SubjectPath, filesep, 'SAM', filesep, 'Workspace_SAM.mat'));
            save(strcat(SubjectPath, filesep, 'SAM', filesep, 'Workspace_SAM.mat'), 'ActIndex', 'ActWgts', 'SAMHeader', 'avgBL') % ToDo: kann nach einem Run entfernt werden 

        otherwise 
            PathSAM = strcat(SubjectPath, filesep, 'SAM', filesep);
            cd (PathSAM) % change matlab dir to Path of weights (SAMdir)
            [SAMHeader, ActIndex, ActWgts] = readWeights('M400,1-50Hz,VGa.wts');
    end        
            
    cd(strcat(SubjectPath, filesep, 'TTest'))
        
    if exist('ERF_1_800ms_Trial_1+orig.BRIK', 'file') % ToDo: kann nach löschen rausgenommen werden
        return
    end

    for i=1:length(avg.trial(:,1,1))
        vs = [];
        ns = [];
        vs=ActWgts*squeeze(avg.trial(i,:,:));
        ns=mean(abs(ActWgts),2);
        vs=vs./repmat(ns,1,size(vs,2));

        fs = 1017.25;
        offset_samples = 509;
        vs_1_1000ms = vs(:,offset_samples:size(vs,2));
        vs_IntOfIn=vs_1_1000ms; % von 0 bis 1000ms

        cfg         = [];
        cfg.step    = 5;
        cfg.boxSize = [-120 120 -90 90 -20 150];
        cfg.prefix  = strcat('ERF_1000ms_Trial_', num2str(i)); % change prefix
        cfg.torig   = 1;   %  comment if you want to sum up activity of specific time intervall
        cfg.TR      = 1/1.01725; % comment if you want to sum up activity of specific time intervall
        VS2Brik(cfg,1e+12*abs(vs_IntOfIn)); % =>creates ERF+orig.Brik + Head
    end   
 end


function [VlrAll, Vall] = get_V (SubjectPath, SubjectName, VlrAll, Vall, TimeInt)

    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg')  % ToDo: zzz_ht später entfernen
        return
    end

    Path = strcat(SubjectPath, filesep, 'TTest');
    CleanData=dir(fullfile(Path, 'CleanData*.mat'))
    load (strcat(SubjectPath,  filesep, 'TTest', filesep, CleanData.name))
    
    for i= 1:length(CleanData.trial)
        
        FileName=strcat(Path, filesep, 'ERF_1000ms_Trial_', num2str(i), '+orig');
        [V, Info] = BrikLoad (FileName);

        % Sum Samples in Time Interval
        fs        = 1017.25;             % Samplingrate vom MEG
        time      = (1:1000)/fs;         %
        TimeBeg   = nearest(time, TimeInt(1))
        TimeEnd   = nearest(time, TimeInt(2))
        V_TimeInt = V(:,:,:, TimeBeg:TimeEnd); 
        % sum forth dimension:
        V_TimeInt_sum=sum(V_TimeInt(:,:,:,:),4);

        clear V V_TimeInt       

        Vlr=flipdim(V_TimeInt_sum,2);
        if i==1
            Vall   = V_TimeInt_sum;
            VlrAll = Vlr;
        else
            Vall(:,:,:,i)   = V_TimeInt_sum;
            VlrAll(:,:,:,i) = Vlr;
        end
        clear V_TimeInt_sum V
    end
end

%%

function TTestLR (SubjectPath, SubjectName, VlrAll, Vall, TimeInt)
 
    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg')  % ToDo: zzz_ht später entfernen
        return
    end
    T=zeros(size(Vall,1),size(Vall,2),size(Vall,3));
    
    for i=1:size(Vall,1)
        for j=1:size(Vall,2)
            for k=1:size(Vall,3)
                if Vall(i,j,k)>0
                    [~, p]  = ttest(squeeze(Vall(i,j,k,:)),squeeze(VlrAll(i,j,k,:)));
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
    
    Path = strcat(SubjectPath, filesep, 'TTest'); 
    cd (Path)
    [~, Info] = BrikLoad (strcat(SubjectPath, filesep, 'TimeIntervalls', filesep, 'ERF_noise_abs_0.32-0.6s_', SubjectName, '+orig')); % Info would be 
    
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('ttest', '_', 'LR', '_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's');
    OptTSOut.verbose = 1;
    OptTSOut.View = '+orig' ;
    %Vsymm=double(Vlr+V>0);
    WriteBrik (T, Info, OptTSOut);   
end

function deleteFiles(SubjectPath, SubjectName)
    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg') || 1 == strcmp (SubjectName, 'zzz_ht') % ToDo: zzz_ht später entfernen
        return
    end

    ERFFiles=dir(fullfile(SubjectPath, filesep, 'TTest', filesep, 'ERF_*.BRIK'))

    for j=1:length(ERFFiles)
        delete (strcat(SubjectPath, filesep, 'TTest', filesep, 'ERF_1000ms_Trial_', num2str(j), '+orig.BRIK'))
        delete (strcat(SubjectPath, filesep, 'TTest', filesep, 'ERF_1000ms_Trial_', num2str(j), '+orig.HEAD'))
    end
end

function kh_TTest_normalize (SubjectPath, SubjectName, TimeInt)

    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg')
        return
    end
    
    PathName = strcat(SubjectPath, filesep, 'TTest');
    cd (PathName)

    FileTTest = strcat('ttest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's', '+orig');
    File_OrthoMNI=strcat(SubjectPath, filesep, 'keptTrials', filesep, 'orthoMNI_avg152T+tlrc');
    
    disp(['!@auto_tlrc -apar ', File_OrthoMNI,' -suffix MNI -input ', FileTTest, ' -dxyz 5 -ok_notice']);
    eval(['!@auto_tlrc -apar ', File_OrthoMNI,' -suffix MNI -input ', FileTTest, ' -dxyz 5 -ok_notice']);

    FileName = strcat('ttest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 'sMNI+tlrc'); 

    eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileName, ' -prefix ', strcat('brain01_', FileName), ' -exp ', 'b*a']) % um Aktivität auf Hirn einzuschränken
    
    if exist(strcat(SubjectPath, filesep, 'TTest', filesep, 'ttest_LR_', num2str(TimeInt(1)), '_',num2str(TimeInt(2)), 'sMNI+tlrc.BRIK'), 'file')
        delete (strcat(SubjectPath, filesep, 'TTest', filesep, 'ttest_LR_', num2str(TimeInt(1)), '_',num2str(TimeInt(2)), 'sMNI+tlrc.BRIK'))
        delete (strcat(SubjectPath, filesep, 'TTest', filesep, 'ttest_LR_', num2str(TimeInt(1)), '_',num2str(TimeInt(2)), 'sMNI+tlrc.HEAD'))
    end
end

function kh_extractActROI (SubjectPath, SubjectName,  ROI_left, ROI_right, ROIName, TimeInt, Time)

    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg')
        return
    end

    PathTTest = strcat (SubjectPath, filesep, 'TTest', filesep);
    cd (PathTTest)
    Path2TValues = strcat (SubjectPath, filesep, 'TTest', filesep, 'brain01_ttest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 'sMNI+tlrc' );
    [V_TValues, Info_TValues] = BrikLoad (Path2TValues);

    PathMask_left = strcat ('/home/kh/ShareWindows/data/patients/', ROI_left, '+tlrc');
    [Mask_left, Info_MASK_left] = BrikLoad (PathMask_left);
    PathMask_right = strcat ('/home/kh/ShareWindows/data/patients/', ROI_right, '+tlrc');
    [Mask_right, Info_MASK_right] = BrikLoad (PathMask_right);

    Left_Voxels = find(Mask_left==1);
    leftAct=V_TValues(Left_Voxels);
    [signLeftVoxels]=find(leftAct>=.95);
    Right_Voxels = find(Mask_right==1);
    rightAct=V_TValues(Right_Voxels);
    [signRightVoxels]=find(rightAct>=.95);

    % calculate LI:
     LI.signRightVoxels  = length(signRightVoxels);
     LI.signLeftVoxels = length(signLeftVoxels);
     LI.SizeROI = size(Left_Voxels, 1);
     LI.relActLeft = length(signLeftVoxels)./size(Left_Voxels, 1);
     LI.relActRight = length(signRightVoxels)./size(Right_Voxels, 1);
     LI.Max_percchange = (length(signLeftVoxels)-length(signRightVoxels))./length(signRightVoxels);
     LI.Classic = (length(signLeftVoxels)-length(signRightVoxels))./(length(signLeftVoxels)+length(signRightVoxels));
     if LI.Classic >= .2
          LI.Lateralization = 'left';
     elseif LI.Classic <= -.2 
        LI.Lateralization = 'right';
     elseif LI.Classic > -.2 && LI.Classic < .2
     LI.Lateralization = 'bilateral';
     end

     PathFile = strcat (SubjectPath, filesep, 'TTest', filesep, 'LI_', ROIName, 'dil_SumOfSignVoxels_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)),  '_s.mat' );
     save (PathFile, 'LI')
end


function [LI_All] = collect_LI (SubjectPath, SubjectName, TimeInt, LI_All, Time )

    if 1 == strcmp (SubjectName, 'Pat_02_13008rh') || 1 == strcmp (SubjectName, 'Pat_03_13014bg')
        return
    end

    PathLIBroca = strcat (SubjectPath, filesep, 'TTest/', 'LI_', 'Brocadil_SumOfSignVoxels_', TimeInt, '.mat');
    load (PathLIBroca)
    LI_All.(SubjectName).(Time).Broca = LI;
    clear LI

    PathLIWernicke = strcat (SubjectPath, filesep, 'TTest/', 'LI_', 'Wernickedil_SumOfSignVoxels_', TimeInt, '.mat');
    load (PathLIWernicke)
    LI_All.(SubjectName).(Time).Wernicke = LI;
    clear LI
end

function collect_LI_excel (group, Time)

Files = strcat('/home/kh/ShareWindows/data/', group, filesep, group, '_SAM');
Dir = dir(Files);
isub = [Dir(:).isdir];
nameFolds = {Dir(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
load(strcat('/home/kh/ShareWindows/data/', group, filesep, 'LIs_signVoxels_singleTrials.mat'));



for i= 1:length(nameFolds)
LI.(Time).Broca(i,1) = LI_All.(nameFolds{i}).(Time).Broca.Classic;
end


for i= 1:length(nameFolds)
LI.(Time).Broca(i,2) = LI_All.(nameFolds{i}).(Time).Broca.Max_percchange;
end

for i= 1:length(nameFolds)
LI.(Time).Broca(i,3) = LI_All.(nameFolds{i}).(Time).Broca.signLeftVoxels;
end

for i= 1:length(nameFolds)
LI.(Time).Broca(i,4) = LI_All.(nameFolds{i}).(Time).Broca.signRightVoxels;
end

for i= 1:length(nameFolds)
LI.(Time).Broca(i,5) = LI_All.(nameFolds{i}).(Time).Broca.SizeROI;
end

for i= 1:length(nameFolds)
LI.(Time).Broca(i,6) = LI_All.(nameFolds{i}).(Time).Broca.relActLeft;
end


for i= 1:length(nameFolds)
LI.(Time).Broca(i,7) = LI_All.(nameFolds{i}).(Time).Broca.relActRight;

end


for i= 1:length(nameFolds)
LI.(Time).Wernicke(i,1) = LI_All.(nameFolds{i}).(Time).Wernicke.Classic;
end


for i= 1:length(nameFolds)
LI.(Time).Wernicke(i,2) = LI_All.(nameFolds{i}).(Time).Wernicke.Max_percchange;
end

for i= 1:length(nameFolds)
LI.(Time).Wernicke(i,3) = LI_All.(nameFolds{i}).(Time).Wernicke.signLeftVoxels;
end

for i= 1:length(nameFolds)
LI.(Time).Wernicke(i,4) = LI_All.(nameFolds{i}).(Time).Wernicke.signRightVoxels;
end

for i= 1:length(nameFolds)
LI.(Time).Wernicke(i,5) = LI_All.(nameFolds{i}).(Time).Wernicke.SizeROI;
end

for i= 1:length(nameFolds)
LI.(Time).Wernicke(i,6) = LI_All.(nameFolds{i}).(Time).Wernicke.relActLeft;
end


for i= 1:length(nameFolds)
LI.(Time).Wernicke(i,7) = LI_All.(nameFolds{i}).(Time).Wernicke.relActRight;

end




end