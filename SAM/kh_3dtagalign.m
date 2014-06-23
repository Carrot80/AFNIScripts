
function for_all_pat()

Patientpath = '/home/kh/ShareWindows/data/patients/patients_SAM';

Dir = dir(Patientpath);
isub = [Dir(:).isdir]; 
nameFolds = {Dir(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

for i=1:length(nameFolds)
    
kh_3dTagalign (Patientpath, nameFolds{i} )

end


end


function kh_3dTagalign (PatientPath, PatientName)

cd ( strcat(PatientPath, filesep, PatientName))

% 1. 
%  !3dTagalign -master ~/brainhull/master+orig -prefix ./ortho anat+orig 

% 2.:
% hs2afni()

% 3.:
% !3dSkullStrip -input ortho+orig -prefix mask -mask_vol -skulls -o_ply ortho

% 4.:
% !meshnorm ortho_brainhull.ply > hull.shape





 
end