% Create ROI from Atlas:

% see data zzz_wi as example

!3dcalc -a TT_N27_EZ_ML+tlrc -exp 'equals(a,13)' -prefix test
! ~/abin$ 3dcalc -a TT_N27_EZ_ML+tlrc -exp 'or(equals(a,13),equals(a,11))' -prefix test
% ~/abin$ => Ort, in dem Befehl ausgeführt wird, dort wird auch Maske
% abgespeichert
% TT_N27_EZ_ML => Name of Atlas
% 13 => Index of ROI



%
!3dExtrema -volume -sep_dist 30 -mask_file /home/kh/abin/test+tlrc ERF_0.14-0.23s_zzz_wi+tlrc
%/home/kh/abin/test+tlrc = mask file => call it differently (see
%left_Wernicke zzz.wi)

% creates mask:
!3dresample -master ERF_0.32-0.6s_zzz_wi+tlrc -prefix Broca -inset /home/kh/abin/test+tlrc  
% !3dresample -master ERF_0.32-0.6s_zzz_wi+tlrc -prefix Wernicke -inset /home/kh/data/controls_SAM/zzz_wi/SAM/left_wernicke  

% closure: also takes extrema at the border of mask:
!3dExtrema -volume -sep_dist 30 -closure -mask_file Broca+tlrc ERF_0.32-0.6s_zzz_wi+tlrc
!3dExtrema -volume -sep_dist 30 -closure -mask_file Broca+tlrc ERF_0.32-0.6s_zzz_wi+tlrc > brocaMax.txt % write it to textfile


% flip ROI to get mirror of right side:

[V, Info] = BrikLoad ('Broca+tlrc');
    Vlr=flipdim(V,1);
%     Vdif=V-Vlr;
%     Vdif(1:16,:,:)=0;
%     InfoNewTSOut = Info;
%     InfoNewTSOut.RootName = '';
%     InfoNewTSOut.BRICK_STATS = [];
%     InfoNewTSOut.BRICK_FLOAT_FACS = [];
%     InfoNewTSOut.IDCODE_STRING = '';
%     InfoNewTSOut.BRICK_TYPES=3*ones(1,1); % 1 short, 3 float.
%     InfoNewTSOut.DATASET_RANK(2)=1;
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = ['BrocaR'];
    OptTSOut.verbose = 1;
    OptTSOut.View = '+tlrc'
    %Vsymm=double(Vlr+V>0);
    WriteBrik (Vlr, Info, OptTSOut);
    
    % ROIfor Wernicke: make smaller (play around) or other atlas?