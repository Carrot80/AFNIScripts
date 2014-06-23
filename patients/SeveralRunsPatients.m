% wenn 2 Runs von einem Patienten vorhanden:  

!3dcalc -a '/home/kh/ShareWindows/data/patients/patients_SAM/Pat_02_13008rh_1/TimeIntervalls/br_z_transf_brain01ERF_noise_0.32-0.6s_Pat_02_13008rh_1+tlrc' -b '/home/kh/ShareWindows/data/patients/patients_SAM/Pat_02_13008rh_2/TimeIntervalls/br_z_transf_brain01ERF_noise_0.32-0.6s_Pat_02_13008rh_2+tlrc' -exp '(a+b)/2' -prefix 'BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s' 
!3dcopy BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s+tlrc BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s.nii 

!3dcalc -a '/home/kh/ShareWindows/data/patients/patients_SAM/Pat_02_13008rh_2/TimeIntervalls/br_z_transf_brain01ERF_noise_0.32-0.6s_Pat_02_13008rh_2+tlrc' -b '/home/kh/ShareWindows/data/patients/patients_SAM/Pat_02_13008rh_1/TimeIntervalls/br_z_transf_brain01ERF_noise_0.32-0.6s_Pat_02_13008rh_1+tlrc' -exp 'a+b' -prefix 'Sum_BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s' 
!3dcopy Sum_BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s+tlrc Sum_BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s.nii 


!3dcalc -a '/home/kh/ShareWindows/data/patients/patients_SAM/Pat_03_13014bg_1/TimeIntervalls/br_z_transf_brain01ERF_noise_0.32-0.6s_Pat_03_13014bg_1+tlrc' -b '/home/kh/ShareWindows/data/patients/patients_SAM/Pat_03_13014bg_2/TimeIntervalls/br_z_transf_brain01ERF_noise_0.32-0.6s_Pat_03_13014bg_2+tlrc' -exp '(a+b)/2' -prefix 'BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s' 
!3dcopy BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s+tlrc BothRuns_br_z_transf_brain01ERF_noise_0.32-0.6s.nii 
