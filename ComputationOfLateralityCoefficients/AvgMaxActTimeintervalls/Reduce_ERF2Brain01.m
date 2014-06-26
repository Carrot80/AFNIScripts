

function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
   kh_reduceERF2Brain( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .32, .47)
  
         
end

end


function kh_reduceERF2Brain (SubjectPath, SubjectName, TimeBeg, TimeEnd)

% if 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ca') || 1 == strcmp (SubjectPath, '/home/kh/data/controls_SAM/zzz_ht')
%     return
% end

SAMPath = strcat(SubjectPath, filesep, 'SAM');
cd (SAMPath)

FileNameOld = strcat('ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc');

disp(['!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ', FileNameOld, ' -prefix ', strcat('brain01', FileNameOld), ' -exp ', 'b*a'])
eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileNameOld, ' -prefix ', strcat('brain01', FileNameOld), ' -exp ', 'b*a'])


end
