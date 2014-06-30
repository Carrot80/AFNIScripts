% computes SAM covariance and weights
% MarkerFile and MEG data have to be directly in Subject Folder 
% createPARAM('M400','ERF','VG',[0.3 0.47],'VG',[-0.35 0],[1 50],[-0.35 0.8]);
function patients ()
   
    PatientFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';
    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

    ForAll (PatientFolder, 'patients')
    ForAll (ControlsFolder, 'controls')

end

function ForAll(Folder, group)

DIR = dir (Folder)
isub = [DIR(:).isdir]; 
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
kh_SAM( strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, group)

end

end


function kh_SAM(SubjectPath, SubjectName, group)
 
    if 1==strcmp(SubjectName, 'Pat_02_13008rh') || 1==strcmp(SubjectName, 'Pat_03_13014bg') || 1==strcmp(SubjectName, 'Pat_02_13008rh_1')
        return
    end

    SAMDir=strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName, filesep, 'SAM');
    if exist(SAMDir, 'dir')
       return
    end
    
    NewDir=strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName);
    if ~exist(NewDir, 'dir')
       mkdir(NewDir)
    end
    
    if ~exist (strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName, filesep, 'MarkerFile.mrk'))
        copyfile(strcat(SubjectPath, filesep, 'MarkerFile.mrk'), strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName))
    end

    [PathFileName, FileName]=findFileName(SubjectPath,SubjectName)

    if ~exist(strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName, filesep, FileName))
        movefile(PathFileName, strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName))
    end
    
    if ~exist(strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName, filesep, 'hs_file'))
        movefile(strcat(SubjectPath, filesep, 'hs_file'), strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName))
    end
    if ~exist(strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName, filesep, 'config'))
        movefile(strcat(SubjectPath, filesep, 'config'), strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName))
    end
    if ~exist(strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName, filesep, SubjectName, '.rtw'))
        movefile(strcat(SubjectPath, filesep, SubjectName,'.rtw'), strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName))
    end

    if ~exist(strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName, filesep, 'hull.shape'))
        movefile(strcat(SubjectPath, filesep, 'hull.shape'), strcat('/home/kh/ShareWindows/data/SAM_BL_350ms', filesep, group, filesep, SubjectName))
    end
        
    PathfileNameHBcor=strcat('/home/kh/ShareWindows/data/SAM_BL_350ms/', group, filesep, SubjectName, filesep, 'hb_tr_lf_c,rfhp0.1Hz') ;
    PathfileNameHBcorNoTr=strcat('/home/kh/ShareWindows/data/SAM_BL_350ms/', group,filesep, SubjectName, filesep,'hb_lf_c,rfhp0.1Hz') ;
    PathfileNameNoHBcor=strcat('/home/kh/ShareWindows/data/SAM_BL_350ms/', group,filesep, SubjectName, filesep,'tr_lf_c,rfhp0.1Hz') ;
    PathfileNameLfcor=strcat('/home/kh/ShareWindows/data/SAM_BL_350ms/', group,filesep, SubjectName, filesep,'lf_c,rfhp0.1Hz') ;
    if exist (PathfileNameHBcor, 'file')
        fileName='hb_tr_lf_c,rfhp0.1Hz';
    elseif exist (PathfileNameHBcorNoTr, 'file')
        fileName='hb_lf_c,rfhp0.1Hz';
    elseif exist(PathfileNameNoHBcor, 'file')
        fileName='tr_lf_c,rfhp0.1Hz';
    elseif exist(PathfileNameLfcor, 'file')
        fileName='lf_c,rfhp0.1Hz';
    else
        PathFileName = [];
        fileName=[];
    end

    cd (strcat('/home/kh/ShareWindows/data/SAM_BL_350ms/', filesep, group))
   
    % compute covariance
    eval(['!SAMcov64 -r ', SubjectName, ' -d ',fileName,' -m M400 -v'])

    % then compute weights:
    eval(['!SAMwts64 -r  ', SubjectName, ' -d ', fileName,' -m M400 -c VGa -v'])
end

function [PathFileName, FileName]=findFileName(SubjectPath, SubjectName)
   
     PathfileNameHBcor=strcat(SubjectPath, filesep,'hb_tr_lf_c,rfhp0.1Hz') ;
     PathfileNameHBcorNoTr=strcat(SubjectPath, filesep,'hb_lf_c,rfhp0.1Hz') ;
     PathfileNameNoHBcor=strcat(SubjectPath, filesep,'tr_lf_c,rfhp0.1Hz') ;
     PathfileNameLfcor=strcat(SubjectPath, filesep,'lf_c,rfhp0.1Hz') ;

     if exist (PathfileNameHBcor, 'file')
         PathFileName=PathfileNameHBcor;
         FileName='hb_tr_lf_c,rfhp0.1Hz';
     elseif exist (PathfileNameHBcorNoTr, 'file')
         PathFileName=PathfileNameHBcorNoTr;
         FileName='hb_lf_c,rfhp0.1Hz';
     elseif exist(PathfileNameNoHBcor, 'file')
         PathFileName=PathfileNameNoHBcor;
         FileName='tr_lf_c,rfhp0.1Hz';
     elseif exist(PathfileNameLfcor, 'file')
         PathFileName=PathfileNameLfcor;
         FileName='lf_c,rfhp0.1Hz';
     else         
         PathFileName = [];
         FileName=[];
     end

end
