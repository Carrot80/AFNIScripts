% flip ROI to get mirror of right side:

!3dcopy AAL_left_AngSuprTempSup_dil_3vox3D.nii AAL_left_AngSuprTempSup_dil_3vox3D
 
[V, Info] = BrikLoad ('Wernicke_left_dil_fMRI+tlrc');
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
    OptTSOut.Prefix = ['Wernicke_right_dil_fMRI'];
    OptTSOut.verbose = 1;
    OptTSOut.View = '+tlrc'
    %Vsymm=double(Vlr+V>0);
    WriteBrik (Vlr, Info, OptTSOut);
    
    

%     !3dcalc -a Wernicke_right_dil_fMRI+tlrc -exp 'ispositive(a)' -prefix Wernicke_right_dil_fMRI_recalc 