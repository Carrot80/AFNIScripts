function forAll()
    % dies für alle Patienten nutzen und umbauen
    ControlsFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    TimeInt = [.32, .60];

    for i= 1:size(nameFolds)
       kh_UTest_normalize( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, TimeInt)
    end
end



function kh_UTest_normalize (SubjectPath, SubjectName, TimeInt)

    % if 0 == strcmp (SubjectName, 'Pat_02_13008rh')
    %     return
    % end


    % Path = strcat (SubjectPath, filesep, 'keptTrials');

    PathName = strcat(SubjectPath, filesep, 'UTest');
    cd (PathName)

    if exist (strcat(SubjectPath, filesep, 'orthoMNI_avg152T+tlrc.BRIK'))
        movefile(strcat(SubjectPath, filesep, 'orthoMNI_avg152T+tlrc.BRIK'), strcat(SubjectPath, filesep,'keptTrials', filesep, 'orthoMNI_avg152T+tlrc.BRIK'));
    end

    if exist(strcat(SubjectPath, filesep, 'orthoMNI_avg152T+tlrc.HEAD'))
        movefile(strcat(SubjectPath, filesep, 'orthoMNI_avg152T+tlrc.HEAD'), strcat(SubjectPath, filesep,'keptTrials', filesep, 'orthoMNI_avg152T+tlrc.HEAD'));
    end

    PathUtest = strcat('Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 's', '+orig');

    disp(['!@auto_tlrc -apar ', '/home/kh/data/patients_SAM/', SubjectName, '/UTest/orthoMNI_avg152T+tlrc',' -suffix MNI -input ', PathUtest, ' -ok_notice']);
    % eval(['!@auto_tlrc -apar ', 'orthoMNI_avg152T+tlrc',' -suffix MNI -input ', PathUtest, ' -dxyz 5']);
    eval(['!@auto_tlrc -apar ', 'orthoMNI_avg152T+tlrc',' -suffix MNI -input ', PathUtest, ' -dxyz 5']);


    eval(['!@auto_tlrc -apar ', SubjectPath, '/UTest/orthoMNI_avg152T+tlrc',' -suffix MNI -input ', PathUtest, ' -dxyz 5']);
    % eval(['!@auto_tlrc -apar ', 'orthoMNI_avg152T+tlrc',' -suffix MNI -input ', PathUtest, ' -dxyz 5']);


    % !@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ERF+orig -dxyz 5  

    FileName = strcat('Utest_LR_', num2str(TimeInt(1,1)), '_', num2str(TimeInt(1,2)), 'sMNI+tlrc'); 

    eval(['!3dcalc -a /home/kh/ShareWindows/data/mniBrain01+tlrc -b ', FileName, ' -prefix ', strcat('brain01_', FileName), ' -exp ', 'b*a'])

end
