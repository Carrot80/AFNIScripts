function forAll()

% dies für alle Patienten nutzen und umbauen
ControlsFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i= 1:size(nameFolds)
    
   kh_SAM_normalize( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i})
    
end


end



function kh_SAM_normalize (SubjectPath, SubjectName)


% SubjectToAnalyse = '/home/kh/data/controls_SAM/zzz_wi';
% if 1 == strcmp (SubjectPath, SubjectToAnalyse)
%     return
% end


Path = strcat (SubjectPath);
cd (Path)

%  Path_Ortho = strcat(SubjectPath, filesep, 'SAM', filesep, 'ortho+orig');

eval (['!@auto_tlrc -base /home/kh/linux_openmp_64/MNI_avg152T1+tlrc -suffix MNI_avg152T -input ortho+orig -ok_notice']);

% !@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ERF+orig -dxyz 5  

end
