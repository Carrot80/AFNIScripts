function ForAll()

% multiply weigths with avgtrials and normalize
ControlsFolder = '/home/kh/ShareWindows/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %# returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 2:size(nameFolds, 1)
    
   kh_SAM( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i})

end


end




function kh_SAM(SubjectPath, SubjectName)

% if exist(Path_vs, 'file')
%     return
% end


% change matlab dir to Path of weights (SAMdir):
PathSAM = strcat(SubjectPath, filesep, 'UTestNoNoise', filesep);
cd (PathSAM)

load (strcat (SubjectPath, filesep, 'SAM/Workspace_SAM.mat'));


% load avg:
PathAVG = strcat('/home/kh/ShareWindows/data/AVG_keeptrials', filesep, 'avgTrials_', SubjectName, '.mat');
load(PathAVG)


for i=1:length(avgTrials.trial(:,1,1))
    
    vs = [];
    ns = [];
    vs=ActWgts*squeeze(avgTrials.trial(i,:,:));
%     ns=mean(abs(ActWgts),2);
%     vs=vs./repmat(ns,1,size(vs,2));
    
    cfg         = [];
    cfg.step    = 5;
    cfg.boxSize = [-120 120 -90 90 -20 150];
    cfg.prefix  = strcat('ERF_avgTrials_', num2str(i)); % change prefix
    cfg.torig   = -500;   %  comment if you want to sum up activity of specific time intervall
    cfg.TR      = 1/1.01725; % comment if you want to sum up activity of specific time intervall
    % VS2Brik(cfg,vs);
    % max(max(vs))
%     VS2Brik(cfg,1e+13*abs(vs)); % =>creates ERF+orig.Brik+Head

    VS2Brik(cfg,1e+13*abs(vs));   
    
    
%     MNIFile = strcat (SubjectPath, filesep, 'SAM', filesep, 'orthoMNI_avg152T+tlrc');
    
%     eval ([ '!@auto_tlrc -apar ', MNIFile, ' -input ', 'ERF_avgTrials_', num2str(i), '+orig -dxyz 5'  ]);
%     
  
    
end

% oder nur einzelne files speichern

end

