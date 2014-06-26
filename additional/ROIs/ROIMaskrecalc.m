function for_all ()

   
%     kh_recalc('R_Broca_3voxdil3D')
    kh_recalc( 'L_Broca_3voxdil3D')
%     kh_recalc('R_AngSuprTempSup_3voxdil3D')
%     kh_recalc( 'L_AngSuprTempSup_3voxdil3D')

    

end

function kh_recalc(ROIold)

Path2ROI = '/home/kh/data/ROIsForAfni/AAL';
cd (Path2ROI)

% rechtshemisphärische ROIs müssen einsen erhalten:
eval(['!3dcalc -a ', ROIold, '+tlrc', ' -exp ''ispositive(a)'' -prefix ', ROIold, '_recalc' ])
disp(['!3dcalc -a ', ROIold, '+tlrc', ' -exp ''ispositive(a)'' -prefix ', ROIold, '_recalc' ])

% dann Anpassung auf größe von MNI_avg152T1, so dass Maske nicht übersteht:
disp(['!3dcalc -a ', ROIold, '_recalc+tlrc', ' -b ', 'MNI_avg152T1+tlrc', ' -exp ''a*b'' -prefix ', ROIold, '_recalc_recalc' ])
eval(['!3dcalc -a ', ROIold, '_recalc+tlrc', ' -b ', 'MNI_avg152T1+tlrc', ' -exp ''a*b'' -prefix ', ROIold, '_recalc_recalc' ])

% nun besteht Maske wieder aus Werten zwischen o und 1, deshalb wieder
% positive Werte in einsen umkehren:
eval(['!3dcalc -a ', ROIold, '_recalc_recalc+tlrc', ' -exp ''ispositive(a)'' -prefix ', ROIold, 'recalc_recalc_recalc' ])

eval(['!3dcalc -a ', ROIold, '_recalc+tlrc', ' -b ', 'MNI_avg152T1+tlrc', ' -exp ''a*b'' -prefix ', ROIold, '_recalc_recalc' ])

!3dresample -master ERF_noise_0.32-0.6s_Pat_01_13021km+tlrc -prefix mniBrain_1mm -inset /home/kh/abin/MNI_avg152T1+tlrc
!3dcalc -a mniBrain_1mm+tlrc -exp 'ispositive(a)' -prefix mniBrain_1mm_01

!3dcalc -a /home/kh/data/mniBrain01+tlrc -b ERF_0.2-0.6s_zzz_ca+tlrc -prefix mERF_0.2-0.6s_zzz_ca -exp 'b*a'

end