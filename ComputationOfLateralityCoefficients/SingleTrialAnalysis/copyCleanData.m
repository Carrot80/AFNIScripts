% kopiert CleanData files in den ShareWindowsOrdner

function DicomImportMain  

    PatientFolder = 'D:\kirsten_thesis\data\patients\';
    ControlsFolder = 'L:\kirsten_thesis\data\controls\' % auf Rechner aus Platzgründen gelöscht
    
%     DicomImport (PatientFolder, 'patients' )
    DicomImport (ControlsFolder, 'controls')
    
end


function DicomImport(Mainfolder, group)
    
     List = dir( Mainfolder );
      
    for i = 1:size (List)
        if ( 0 == strcmp( List(i,1).name, '.') && 0 == strcmp( List(i,1).name, '..') )
            CopyFiles( strcat(Mainfolder, List(i,1).name), List(i,1).name, group  );
        end
    end
end

function CopyFiles(SubjectPath, SubjectName, group)


if 1==strcmp(SubjectName, 'Pat_02_13008rh')
    File_CleanData1= strcat (SubjectPath, filesep, 'MEG\01_Input_no_noisereduction\1', filesep, 'CleanData.mat');
    NewDir=strcat('D:\Arbeit\LinuxExchange\data\', group, filesep, group, '_SAM', filesep, SubjectName, '_1', filesep, 'TTest')
    if ~exist(NewDir, 'dir')
        mkdir (NewDir)
    end
    copyfile(File_CleanData1, NewDir)
    
    File_CleanData2= strcat (SubjectPath, filesep, 'MEG\01_Input_no_noisereduction\2', filesep, 'CleanData.mat');
    NewDir2=strcat('D:\Arbeit\LinuxExchange\data\', group, filesep, group, '_SAM', filesep, SubjectName, '_2', filesep, 'TTest')
    if ~exist(NewDir2, 'dir')
        mkdir (NewDir2)
    end
    copyfile(File_CleanData2, NewDir2)
    
    return
    
elseif 1==strcmp(SubjectName, 'Pat_03_13014bg')
    File_CleanData1= strcat (SubjectPath, filesep, 'MEG\02_PreProcessing_1', filesep, 'CleanDataSNR.mat');
    NewDir=strcat('D:\Arbeit\LinuxExchange\data\', group, filesep, group, '_SAM', filesep, SubjectName, '_1', filesep, 'TTest')
    if ~exist(NewDir, 'dir')
        mkdir (NewDir)
    end
    copyfile(File_CleanData1, NewDir)
    
    File_CleanData2= strcat (SubjectPath, filesep, 'MEG\02_PreProcessing_2', filesep, 'CleanDataSNR_best.mat');
    NewDir2=strcat('D:\Arbeit\LinuxExchange\data\', group, filesep, group, '_SAM', filesep, SubjectName, '_2', filesep, 'TTest')
    if ~exist(NewDir2, 'dir')
        mkdir (NewDir2)
    end
    copyfile(File_CleanData2, NewDir2)
    
    return
    
elseif 1==strcmp(SubjectName, 'Pat_07_13033gc')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData_nobadTrls.mat');
elseif 1==strcmp(SubjectName, 'Pat_08_13026pj')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData_rejcomp.mat');
elseif 1==strcmp(SubjectName, 'Pat_11_13030rs')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData_rejcomp.mat');
elseif 1==strcmp(SubjectName, 'Pat_14_13039sg')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanDataSNR2.mat');
elseif 1==strcmp(SubjectName, 'Pat_17_13060ec')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData_rejectcomp2.mat');
elseif 1==strcmp(SubjectName, 'Pat_19_13055eg')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData_rejcomp.mat');
elseif 1==strcmp(SubjectName, 'Pat_21_13056hz')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData_rejcomp108Trials.mat');
elseif 1==strcmp(SubjectName, 'Pat_22_13059oc')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData_rejcomp2.mat');
elseif 1==strcmp(SubjectName, 'Pat_24_13067pj')
    File_CleanData= strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanDataSNR1_5.mat');
else
    File_CleanData = strcat (SubjectPath, filesep, 'MEG\', '02_PreProcessing', filesep, 'CleanData.mat');
    
end
    
NewDir=strcat('D:\Arbeit\LinuxExchange\data\', group, filesep, group, '_SAM', filesep, SubjectName, filesep, 'TTest')
if ~exist(NewDir, 'dir')
    mkdir (NewDir)
end
copyfile(File_CleanData, NewDir)


end