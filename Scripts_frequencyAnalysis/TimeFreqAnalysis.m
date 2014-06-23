function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];


for i= 1:size(nameFolds)
    
   FreqAnalysis( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i})
    
end


end


function FreqAnalysis (SubjectPath, SubjectName)


% SubjectToAnalyse = '/home/kh/data/controls_SAM/zzz_ca';
% if 1 == strcmp (SubjectPath, SubjectToAnalyse)
%     return
% end

[fileName] = findFileName(SubjectPath, SubjectName)
load (strcat (SubjectPath, filesep, 'CleanTrl_1_50Hz.mat'));

cfg         = [];
cfg.dataset = fileName;
cfg.trl     = trl ;
cfg.demean  = 'yes';

data = ft_preprocessing(cfg);


cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 5:100;                            % freq of interest 3 to 100Hz
cfg.t_ftimwin    = 1./cfg.foi;
cfg.toi          = -0.3:0.01:0.7;                    % time window "slides" from -0.1 to 0.5 sec in steps of 0.02 sec (20 ms)
cfg.channel='MEG';
cfg.tapsmofrq =ones(1,length(cfg.foi));
%cfg.feedback='no';
TF = ft_freqanalysis(cfg, data);
% now plot one channel


Path = strcat(SubjectPath, filesep, 'FreqAnalysis', filesep, 'TF');
save (Path, 'TF')

cfg=[];
cfg.layout       = '4D248.lay';
cfg.interactive='yes';
cfg.baseline         = [-.3 0];%'yes';
cfg.baselinetype     = 'relative';
ft_multiplotTFR(cfg, TF);


cfg=[];
cfg.layout       = '4D248.lay';
cfg.interactive='yes';
cfg.xlim=[0.3 .5];
cfg.ylim = [10 10];
ft_topoplotTFR(cfg, TF);


end


function [fileName]=findFileName(SubjectPath, SubjectName)

   
     PathfileNameHBcor=strcat(SubjectPath, filesep,'hb_tr_lf_c,rfhp0.1Hz') ;
     PathfileNameHBcorNoTr=strcat(SubjectPath, filesep,'hb_lf_c,rfhp0.1Hz') ;
     PathfileNameNoHBcor=strcat(SubjectPath, filesep,'tr_lf_c,rfhp0.1Hz') ;
     PathfileNameLfcor=strcat(SubjectPath, filesep,'lf_c,rfhp0.1Hz') ;

     if exist (PathfileNameHBcor, 'file')
         fileName=PathfileNameHBcor;
     elseif exist (PathfileNameHBcorNoTr, 'file')
         fileName=PathfileNameHBcorNoTr;
     elseif exist(PathfileNameNoHBcor, 'file')
         fileName=PathfileNameNoHBcor;
     elseif exist(PathfileNameLfcor, 'file')
         fileName=PathfileNameLfcor;
         
     else 
         
         fileName = [];
     end

end