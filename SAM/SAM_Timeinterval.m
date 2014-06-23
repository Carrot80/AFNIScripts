function forAll()

    % dies für alle Patienten nutzen und umbauen
    ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

    DIR = dir (ControlsFolder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    for i= 1:size(nameFolds)
        kh_SAM_TimeInt( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .320, .470)
        %    kh_SAM_TimeInt( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .200, .310) 
        %    kh_SAM_TimeInt( strcat(ControlsFolder, filesep, nameFolds{i,1}), nameFolds{i}, .200, .600) 
    end
end

function kh_SAM_TimeInt (SubjectPath, SubjectName, IntBeg, IntEnd)
    Path2SAM = strcat (SubjectPath, filesep, 'SAM');
    cd(Path2SAM);
    Path2VS = strcat (SubjectPath, filesep, 'SAM', filesep, 'Workspace_SAM.mat');
    load (Path2VS) ;

    fs = 1017.25;
    size_vs=size(vs);

    offset_samples = 509;
    vs_1_1000ms = vs(:,offset_samples:size_vs(2));
    size_vs_1_1000ms=size(vs_1_1000ms);
    time_samples=1:size_vs_1_1000ms(2);
    time_sec=time_samples./fs;
    % I could have used function nearest, too.

    figure
    plot(time_sec,max(abs(vs_1_1000ms)));
    axis tight;
    Title=strcat(num2str(IntBeg), '_', num2str(IntEnd), 's_', SubjectName);
    title(Title)
    % print (Title, 'fig'); 
    close all

    sample_int_Beg=size(find(time_sec<=IntBeg));
    sample_int_End=size(find(time_sec<=IntEnd));

    vs_IntOfIn=vs_1_1000ms(:,sample_int_Beg(2):sample_int_End(2));
    sum_vs_IntOfIn = abs(sum(vs_IntOfIn')); %falsch! 

    % Save it to see it in afni
    cfg=[];
    cfg.step=5;
    cfg.boxSize=[-120 120 -90 90 -20 150];
    str_timeInt= strcat('ERF_noise_abs_', '_', num2str(IntBeg), '-', num2str(IntEnd), 's_', SubjectName);
    cfg.prefix = str_timeInt; % change prefix
    % cfg.torig=-500;   %  comment if you want to sum up activity of specific time intervall
    % cfg.TR=1/1.01725; % comment if you want to sum up activity of specific time intervall
    % VS2Brik(cfg,vs);
    % max(max(vs))
    VS2Brik(cfg,1e+13*(sum_vs_IntOfIn')); % =>creates ERF+orig.Brik+Head 

    NewFileName = strcat(str_timeInt,'+orig');
    %  disp(['!@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ', NewFileName,' -dxyz 5']) % 
    eval(['!@auto_tlrc -apar orthoMNI_avg152T+tlrc -input ', NewFileName,' -dxyz 5']) % 
end
