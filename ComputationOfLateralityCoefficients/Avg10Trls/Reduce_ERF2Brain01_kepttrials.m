

function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/data/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
   kh_reduceERF2Brain( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i})
  
         
end

end


function kh_reduceERF2Brain (SubjectPath, SubjectName)

% if 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ca') || 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ht')
%     return
% end

SAMPath = strcat(SubjectPath, filesep, 'keptTrials/');
cd (SAMPath)

for i = 1 : 10
    
    FileName = strcat('ERF_avgTrials_', num2str(i),'+tlrc');
    
    disp(['!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ', FileName, ' -prefix ', strcat('brain01_', FileName), ' -exp ', 'b*a'])
    eval(['!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ', FileName, ' -prefix ', strcat('brain01_', FileName), ' -exp ', 'b*a'])
    
end

end
