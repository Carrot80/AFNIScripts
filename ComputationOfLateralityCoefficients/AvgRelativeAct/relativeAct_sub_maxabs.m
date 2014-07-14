% computes  LI (sqrt) of maximum of absolute values
function patients ()
   % funkioniert zur Lateralisierung nicht gut
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
    TimeInt=[.32, .47];
    for i= 1:size(nameFolds)
       kh_SAM_RelAct( strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, TimeInt, group)
       kh_extractActROI(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, 'Broca_left_dil', 'Broca_right_dil', 'Broca', TimeInt, group)
       kh_extractActROI(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, 'Wernicke_left_dil', 'Wernicke_right_dil', 'Wernicke', TimeInt, group)
      [LI_All_RelAct]= collect_LI (strcat(Folder, filesep, nameFolds{i,1}), i, LI_All_RelAct, nameFolds{i,1}, TimeInt)
    end
    
    Path_LI_All=strcat('/home/kh/ShareWindows/Results/LI_RelativeActSubmaxabs/LI_All_RelAct_', num2str(TimeInt(1)), '_', num2str(TimeInt(2)), 's_', group, '.mat')
    save (Path_LI_All, 'LI_All_RelAct')

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
   

    NewDir=strcat(SubjectPath, filesep, 'RelativeActSubmaxabs')
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
%     kh_z_transform (SubjectPath, SubjectName, TimeBeg, TimeEnd)
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


function kh_z_transform (SubjectPath, SubjectName, TimeBeg, TimeEnd)
    
    FileName = strcat('br01RelAct_', num2str(TimeBeg), '-', num2str(TimeEnd), 's+tlrc');   

    [V, Info] = BrikLoad (FileName);
    Vz=(V-mean(V(:)))/std(V(:));   
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('z_transf_', FileName);
    OptTSOut.verbose = 1;
    OptTSOut.View = '+tlrc'
    WriteBrik (Vz, Info, OptTSOut);
  
    FileNameNew = OptTSOut.Prefix;     
    eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileNameNew,  ' -prefix ', strcat('br_', FileNameNew),' -exp ' , 'b*a'])

    PathERF = strcat('br_z_transf_brain01ERF_noise_abs_', num2str(TimeBeg), '-', num2str(TimeEnd), 's_', SubjectName, '+tlrc');
      
end



function kh_extractActROI (SubjectPath, SubjectName, ROI_left, ROI_right, ROI, TimeInt, group)

    if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg') 
        return
    end
    cd (strcat(SubjectPath, filesep, 'RelativeActSubmaxabs'))
    FileName = strcat('br01RelAct_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's+tlrc');
    [V_ERF, Info_ERF] = BrikLoad (FileName);

    PathMask_left = strcat ('/home/kh/ShareWindows/data/patients/', ROI_left, '+tlrc');
    [Mask_left, Info_MASK_left] = BrikLoad (PathMask_left); 

    PathMask_right = strcat ('/home/kh/ShareWindows/data/patients/', ROI_right, '+tlrc');
    [Mask_right, Info_MASK_right] = BrikLoad (PathMask_right);

    for i= 1:length(V_ERF)
        LeftAct(:,:,:) = Mask_left.*V_ERF(:,:,:);
    end

    for i= 1:length(V_ERF)
        RightAct(:,:,:) = Mask_right.*V_ERF(:,:,:);
    end

    %%
    maxabs_RightAct = max(abs(RightAct(:)));
    maxabs_LeftAct = max(abs(LeftAct(:)));
    
    LI_maxabs=(maxabs_LeftAct-maxabs_RightAct)/(maxabs_LeftAct+maxabs_RightAct);
    LI_squared=(maxabs_LeftAct^2-maxabs_RightAct^2)/(maxabs_LeftAct^2+maxabs_RightAct^2);
    LI_sqrt=sqrt(abs(LI_squared));
    if LI_squared <0
        LI_maxabs_squared=LI_sqrt*(-1);
    else LI_maxabs_squared=LI_sqrt;
    end

    LI.LI_maxabs=LI_maxabs
    LI.LI_maxabs_squared=LI_maxabs_squared
    
    LI_Path=strcat(SubjectPath, filesep, 'RelativeActSubmaxabs', filesep, 'LI_maxabs', ROI, '_', num2str(TimeInt(1)),'-', num2str(TimeInt(2)),'s.mat');
    save (LI_Path, 'LI') 
end


function [LI_All_RelAct]=collect_LI (SubjectPath, i, LI_All_RelAct, SubjectName, TimeInt )

if 1==strcmp(SubjectName,'Pat_02_13008rh') || 1==strcmp(SubjectName,'Pat_03_13014bg') 
    return
end

load(strcat(SubjectPath, filesep,  'RelativeActSubmaxabs', filesep, 'LI_maxabsBroca_', num2str(TimeInt(1)), '-',num2str(TimeInt(2)), 's.mat'))
LI_All(1,1)=LI.LI_maxabs;
clear LI.LI_maxabs
LI_All(1,2)=LI.LI_maxabs_squared;
clear LI.LI_maxabs_squared

load(strcat(SubjectPath, filesep,  'RelativeActSubmaxabs', filesep, 'LI_maxabsWernicke_', num2str(TimeInt(1)), '-',num2str(TimeInt(2)), 's.mat'))
LI_All(1,3)=LI.LI_maxabs;
clear LI.LI_maxabs

LI_All(1,4)=LI.LI_maxabs_squared;
clear LI.LI_maxabs_squared

LI_All_RelAct(i,:)=LI_All
end




