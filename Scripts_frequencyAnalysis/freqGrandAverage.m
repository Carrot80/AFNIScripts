function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
TF_All = [];

for i= 1:size(nameFolds)
    
   [TF_All] = FreqAnalysis( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, i, TF_All)
    

     
end



GrandAVGTF = freqgrandaverage (TF_All.('zzz_ca'), TF_All.('zzz_ht'), TF_All.('zzz_ka'), TF_All.('zzz_md'), TF_All.('zzz_mf'), TF_All.('zzz_ms'), TF_All.('zzz_sc'), TF_All.('zzz_sf'), TF_All.('zzz_si'), TF_All.('zzz_wi'))
save GrandAVGTF GrandAVGTF

    cfg=[];
    cfg.layout       = '4D248.lay';
    cfg.interactive='yes';
    cfg.baseline         = [-.3 0];%'yes';
    cfg.baselinetype     = 'relative';
    ft_multiplotTFR(cfg, GrandAVGTF);


    cfg=[];
    cfg.layout       = '4D248.lay';
    cfg.interactive='yes';
    cfg.xlim=[0.3 .5];
    cfg.ylim = [10 10];
    cfg.baseline         = [-.3 0];%'yes';
    cfg.baselinetype     = 'relative';
    ft_topoplotTFR(cfg, GrandAVGTF);


end


function [TF_All] = FreqAnalysis (SubjectPath, SubjectName, i, TF_All)


load ( strcat(SubjectPath, filesep, 'FreqAnalysis', filesep, 'TF'));

TF_All.(SubjectName) = TF;
clear TF

end

