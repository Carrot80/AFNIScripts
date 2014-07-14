function patients ()
   
    PatientFolder = '/home/kh/ShareWindows/data/SAM_BL_350ms/patients';
    ControlsFolder = '/home/kh/ShareWindows/data/SAM_BL_350ms/controls';

    forAll (PatientFolder, 'Patients')
    forAll (ControlsFolder, 'Controls')

end


function forAll(Folder, group)
    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    LI_All_RelAct=zeros(size(nameFolds), 4);
    TimeInt=[.32, .6];
    LI_ALL_zscores=zeros(24, 24);
    
    for i= 1:size(nameFolds)
       kh_SAM_RelAct( strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, TimeInt, group)
       [LI_ALL_zscores]=kh_extractActROI(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt, group, i, LI_ALL_zscores, 1:4, 9:12, 17:20)
       [LI_ALL_zscores]=kh_extractActROI(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt, group, i , LI_ALL_zscores, 5:8, 13:16, 21:24)
    end
    
    Path_LI_All=strcat('/home/kh/ShareWindows/Results/LI_RelativeActSubstraction/zscores_LI_RelAct_', num2str(TimeInt(1)), '_', num2str(TimeInt(2)), 's_', group, '.mat')
    save (Path_LI_All, 'LI_ALL_zscores')

end

function kh_SAM_RelAct (SubjectPath, SubjectName, TimeInt, group)

    if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg') 
        return
    end
    
    SAMPath = strcat(SubjectPath, filesep, 'SAM');
    cd (SAMPath)
    %read weights:
    [SAMHeader, ActIndex, ActWgts]=readWeights('M400,1-50Hz,VGa.wts');
    % load avg:
    OldPath = strcat('/home/kh/ShareWindows/data/', group, filesep, group, '_SAM', filesep, SubjectName, filesep); 
    
    switch SubjectName
        case 'Pat_03_13014bg_1'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanDataSNR.mat');
        case 'Pat_03_13014bg_2'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanDataSNR_best.mat');
        case  'Pat_07_13033gc'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanData_nobadTrls.mat');
        case  'Pat_08_13026pj'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanData_rejcomp.mat');
        case  'Pat_11_13030rs'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanData_rejcomp.mat');
        case  'Pat_14_13039sg'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanDataSNR2.mat');
        case  'Pat_17_13060ec'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanData_rejectcomp2.mat');
        case  'Pat_19_13055eg'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanData_rejcomp.mat');
        case  'Pat_21_13056hz'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanData_rejcomp108Trials.mat');
        case  'Pat_22_13059oc'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanData_rejcomp2.mat');
        case  'Pat_24_13067pj'
            File_CleanData= strcat (OldPath,'TTest/', 'CleanDataSNR1_5.mat');
        otherwise
            File_CleanData = strcat(OldPath,'TTest/', 'CleanData.mat');
    end

    % load avg:
    load(File_CleanData)
    CleanData_BL=correctBL(CleanData, [-0.35 -0.03]);
    cfg=[];
    avgBL=ft_timelockanalysis(cfg, CleanData_BL)

    fs = 1017.25;
    % Baselineintervall 350-50ms prestim:
    avgBaseline=avgBL.avg(:,158:490); % in etwa
    avgVG_1_1000=avgBL.avg(:,509:size(avgBL.avg,2)); %in etwa
    time_samples=(1:size(avgVG_1_1000,2))./fs;
    avgVG_TimeInt=avgVG_1_1000(:, nearest(time_samples, TimeInt(1)):nearest(time_samples, TimeInt(2)));
    VS_Baseline=ActWgts*avgBaseline;
    VS_VG=ActWgts*avgVG_TimeInt;
    ns=mean(abs(ActWgts),2); 
    VS_VG_ns=VS_VG./repmat(ns,1,size(VS_VG,2));
    VS_Baseline_ns=VS_Baseline./repmat(ns,1,size(VS_Baseline,2));
    VS_RelActSub=sum(abs(VS_VG_ns'))-sum(abs(VS_Baseline_ns'));

    NewDir=strcat(SubjectPath, filesep, 'RelativeActSubstraction')
    if ~exist(NewDir, 'dir')
        mkdir (NewDir)
    end
    cd(NewDir)

    % Save it to load it in afni
    cfg=[];
    cfg.step=5;
    cfg.boxSize=[-120 120 -90 90 -20 150];
    str_timeInt= strcat('RelAct', '_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's');
    cfg.prefix = str_timeInt; % change prefix
    VS2Brik(cfg,VS_RelActSub'); % =>creates ERF+orig.Brik+Head 

    NewFileName = strcat(str_timeInt,'+orig');
    eval(['!@auto_tlrc -apar ', strcat(OldPath, 'keptTrials', filesep, 'orthoMNI_avg152T+tlrc'), ' -input ', NewFileName,' -dxyz 5']) % 

    kh_reduceERF2Brain (SubjectPath, SubjectName, TimeInt)
    kh_z_transform (SubjectPath, SubjectName, TimeInt)
end


function kh_reduceERF2Brain (SubjectPath, SubjectName, TimeInt)

    if exist(strcat('br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.BRIK'))
        delete(strcat('br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.BRIK'))
        delete(strcat('br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.HEAD'))
    end
    FileNameOld = strcat('RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc');
    eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileNameOld, ' -prefix ', strcat('br01', FileNameOld), ' -exp ', 'b*a'])
    if exist(strcat('RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.BRIK'))
        delete (strcat('RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.BRIK'))
        delete (strcat('RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.HEAD'))
    end
    
    if exist(strcat('RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+orig.BRIK'))
        delete (strcat('RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+orig.BRIK'))
        delete (strcat('RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+orig.HEAD'))
    end
    
end


function kh_z_transform (SubjectPath, SubjectName, TimeInt)
    
    FileName = strcat('br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc');   

    [V, Info] = BrikLoad (FileName);
    Vz=(V-mean(V(:)))/std(V(:));   
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('z_transf_', FileName);
    OptTSOut.verbose = 1;
    OptTSOut.View = '+tlrc'
    WriteBrik (Vz, Info, OptTSOut);
  
    FileNameNew = OptTSOut.Prefix;     
    eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileNameNew,  ' -prefix ', strcat('br_', FileNameNew),' -exp ' , 'b*a'])
    delete (strcat('z_transf_br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.BRIK'))
    delete (strcat('z_transf_br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc.HEAD'))     
end



function [LI_ALL_zscores]=kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt, group, i, LI_ALL_zscores, A, B, C)

    if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg') 
        return
    end
    cd (strcat(SubjectPath, filesep, 'RelativeActSubstraction'))
    FileName = strcat('br_z_transf_br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc');
    [V_ERF, Info_ERF] = BrikLoad (FileName);

    PathMask_left = strcat ('/home/kh/ShareWindows/data/patients/', ROI_left, '+tlrc');
    [Mask_left, Info_MASK_left] = BrikLoad (PathMask_left); 

    PathMask_right = strcat ('/home/kh/ShareWindows/data/patients/', ROI_right, '+tlrc');
    [Mask_right, Info_MASK_right] = BrikLoad (PathMask_right);

    for j= 1:length(V_ERF)
        LeftAct(:,:,:) = Mask_left.*V_ERF(:,:,:);
    end

    for j= 1:length(V_ERF)
        RightAct(:,:,:) = Mask_right.*V_ERF(:,:,:);
    end
    
    %% p>.01
    
    [ind_LeftAct]=find(LeftAct>2.33)
    zscores_LeftAct=LeftAct(ind_LeftAct)

    [ind_RightAct]=find(RightAct>2.33)
    zscores_RightAct=RightAct(ind_RightAct)

    LI.LI_Voxelvalue_p01=(sum(zscores_LeftAct)-sum(zscores_RightAct))./(sum(zscores_LeftAct)+sum(zscores_RightAct))
    LI.LI_Voxelcount_p01=(length(zscores_LeftAct)-length(zscores_RightAct))./(length(zscores_LeftAct)+length(zscores_RightAct))
    LI.Voxelcount_p01_LeftVox=length(zscores_LeftAct);
    LI.Voxelcount_p01_RightVox=length(zscores_RightAct);

    %% p>.001

    [ind_LeftAct_p001]=find(LeftAct>3.090232)
    zscores_LeftAct_p001=LeftAct(ind_LeftAct_p001)

    [ind_RightAct_p001]=find(RightAct>3.090232)
    zscores_RightAct_p001=RightAct(ind_RightAct_p001)

    LI.LI_Voxelvalue_p001=(sum(zscores_LeftAct_p001)-sum(zscores_RightAct_p001))./(sum(zscores_LeftAct_p001)+sum(zscores_RightAct_p001))
    LI.LI_Voxelcount_p001=(length(zscores_LeftAct_p001)-length(zscores_RightAct_p001))./(length(zscores_LeftAct_p001)+length(zscores_RightAct_p001))
    LI.Voxelcount_p001_LeftVox=length(zscores_LeftAct_p001);
    LI.Voxelcount_p001_RightVox=length(zscores_RightAct_p001);

LI_Path=strcat(SubjectPath, filesep, 'LI_', ROI, '_Voxelvalue_noise_abs_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_Path, 'LI') 

%% p>.05

[ind_LeftAct_p05]=find(LeftAct>1.644853)
zscores_LeftAct_p05=LeftAct(ind_LeftAct_p05)

[ind_RightAct_p05]=find(RightAct>1.644853)
zscores_RightAct_p05=RightAct(ind_RightAct_p05)

LI.LI_Voxelvalue_p05=(sum(zscores_LeftAct_p05)-sum(zscores_RightAct_p05))./(sum(zscores_LeftAct_p05)+sum(zscores_RightAct_p05));
LI.LI_Voxelcount_p05=(length(zscores_LeftAct_p05)-length(zscores_RightAct_p05))./(length(zscores_LeftAct_p05)+length(zscores_RightAct_p05));
LI.Voxelcount_p05_LeftVox=length(zscores_LeftAct_p05);
LI.Voxelcount_p05_RightVox=length(zscores_RightAct_p05);

LI_Path=strcat(SubjectPath, filesep, 'LI_', ROI, '_zscores_', num2str(TimeInt(1,1)),'_', num2str(TimeInt(1,2)),'ms.mat')
save (LI_Path, 'LI') 

LI_ALL_zscores(i,A)= [LI.LI_Voxelvalue_p05 LI.LI_Voxelcount_p05 LI.Voxelcount_p05_LeftVox LI.Voxelcount_p05_RightVox ]
LI_ALL_zscores(i,B)= [LI.LI_Voxelvalue_p01 LI.LI_Voxelcount_p01 LI.Voxelcount_p01_LeftVox LI.Voxelcount_p01_RightVox ]
LI_ALL_zscores(i,C)= [LI.LI_Voxelvalue_p001 LI.LI_Voxelcount_p001 LI.Voxelcount_p001_LeftVox LI.Voxelcount_p001_RightVox ]
    
    LI_Path=strcat(SubjectPath, filesep, 'RelativeActSubstraction', filesep, 'LI_Max', ROI, '_', num2str(TimeInt(1)),'-', num2str(TimeInt(2)),'s.mat');
    save (LI_Path, 'LI') 
end

