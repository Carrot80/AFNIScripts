% SAM


function for_all_subjects  

      ControlsFolder = 'D:\kirsten_thesis\data\controls\';
    
     SelectSubjects (ControlsFolder)
%      SelectSubjects (PatientFolder)
end


function SelectSubjects (Mainfolder)

    List = dir( Mainfolder );

 for i = 1 : size (List)
      if ( 0 == strcmp( List(i,1).name, '.') && 0 == strcmp( List(i,1).name, '..') )
          SubjectPath = strcat(Mainfolder, List(i,1).name) ;
          SubjectName = List(i,1).name  
          
        
          [Path] = MakePath(SubjectPath, SubjectName)
          CopyFiles (Path, SubjectName)
          
           
          
          
          
%           GetAVG (SubjectName, Path)
          
      end

 end
 
%            statsMEGPlanar (SubjectName, Path)
 
end

function  CopyFiles (Path, SubjectName)

configFileSAMPath = strcat('D:\kirsten_thesis\data\controls_SAM', filesep, SubjectName,filesep, 'config');

if exist (configFileSAMPath, 'file')
    return
end

[fileName]=PathForFileName(SubjectName, Path);
hs= strcat(Path.DataInput, filesep, 'hs_file');
configFile = strcat(Path.DataInput, filesep, 'config');
BadTrials = strcat(Path.Preprocessing, filesep, 'BadTrials.mat');
BadTrials_New = strcat('D:\kirsten_thesis\data\controls_SAM', filesep, SubjectName, filesep, 'BadTrials_Jumps.mat');

trlsel = strcat(Path.Preprocessing, filesep, 'trlsel.mat');

SAMPATH= (strcat('D:\kirsten_thesis\data\controls_SAM', filesep, SubjectName, filesep));
copyfile(fileName, SAMPATH)
copyfile(hs, SAMPATH)
copyfile(configFile, SAMPATH)
copyfile(BadTrials, BadTrials_New)
copyfile(trlsel, SAMPATH)



end

function [Path] = MakePath(SubjectPath, SubjectName)
        
        Path                     = [];
        Path.Subject             = SubjectPath ; 
        Path.DataInput           = strcat ( SubjectPath, '\MEG\01_Input_no_noisereduction')                 ;
        Path.Preprocessing       = strcat ( SubjectPath, '\MEG\02_PreProcessing')          ;
        Path.StatisticsSensorLevel = strcat (SubjectPath, '\MEG\SensorLevelAnalysis')          ;
end

   

    function      [fileName]=PathForFileName(SubjectName, Path)
        
   
     PathfileNameHBcor=strcat(Path.DataInput, filesep,'hb_tr_lf_c,rfhp0.1Hz') ;
     PathfileNameHBcorNoTr=strcat(Path.DataInput, filesep,'hb_lf_c,rfhp0.1Hz') ;
     PathfileNameNoHBcor=strcat(Path.DataInput, filesep,'tr_lf_c,rfhp0.1Hz') ;
     PathfileNameLfcor=strcat(Path.DataInput, filesep,'lf_c,rfhp0.1Hz') ;

     if exist (PathfileNameHBcor, 'file')
         fileName=PathfileNameHBcor;
     elseif exist (PathfileNameHBcorNoTr, 'file')
         fileName=PathfileNameHBcorNoTr;
     elseif exist(PathfileNameNoHBcor, 'file')
         fileName=PathfileNameNoHBcor;
     elseif exist(PathfileNameLfcor, 'file')
         fileName=PathfileNameLfcor
         
     else 
         
         fileName = [];
 
                 
         
    end
    end