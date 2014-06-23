% Cleaning-Skripte:

% preprocessing Scripte (außerhalb von Linux):

AverageData.m
AVG_keepTrials.m

%

% Übersicht über Linux-Skripte
% Beispielskripte von Yuval:
Script_SAM_Yuval % Yuvals Vorlage für SAM_Beamforming und für Berechnung von virtuellen Sensoren und Multiplikation von Weights mit AVG 


% eigene Skripte in chronologischer Reihenfolge:
kh_MRI2Afni.m % 
kh_SAM_Beamformer % eigenes Skript, um SAM_Beamforming durchzuführen
SAM_normalize % um SAM_results und Orthogonale Anatomie zu normalisieren in MNI avgT152
SAM_Timeinterval % calculates mean activity for time intervals
(ERF_transform2MNIspace.m) % nachträgliche transformation der funktionellen Daten zum MNI space, da vorher falsches template verwendet
Reduce_ERF2Brain01.m % reduces die ERF activity which expands the brain after MNI Normalization back to Brain Volume size
kh_z_scoreNormalizationERF.m 


% U-Tst for all controls (group level):
kh_avg % creating AVG for  all controls (und visuelle Darstellung) => vorher kh_z_scoreNormalizationERF.m 
LR_ttest % U-Test for Left and Right Hemisphere created by Yuval

smaller_timeintervals.m % includes Reduce_ERF2Brain01.m and kh_z_scoreNormalizationERF.m 
kh_avg
LRttest

% Erstellung von ROI-Masken:

Create_ROI.m
FlipROI % flips left ROI to right Hemisphere
ROIMaskrecalc.m % after flipping ROI, one has to recalc right ROI, as it does not consist of values of ones after sampling, but of 0,..7. Recalculated to resample values of 1
ReMasterROI.m % ROI muss resampled werden, damit sie gleiche Auflösung hat

% ROI Analysis
ROI_Analysis.m % looks for extrem values for each ROI and saves to Textfile
loadTextComputeLI.m % computes LIs of maximum values
extractActROI.m % extrahiert Aktivität aller Voxel aus ROI
extractfullActROI


% U-Test für einzelne Subjects: WICHTIG
kh_SAM_Beamforming_keepttrials.m % Wiederholung für AVGTrials, mitteln von 10er Gruppen;
LR_UTest_ind.m % U-Test für verschiedene Zeitintervalle; Linke vs. rechte Hemisphere, hier nohc keine ROIextraction
UTest_normalize.m % NOrmalisierung in MNI space
extract_UValuesROI.m % Berechnung des LI-Wertes für ROIs
kh_collect_LI.m % sammeln der Werte für alle Subjects in einer Matrix

LR_UTest_ind_UValues.m % the same as LR_ttest_ind, except for that I wanted to get stats (z-value)


% Sanity Check without NoiseEstimation: U-Test für einzelne Subjects
kh_SAM_Beamforming_keepttrials_withoutNoiseEst.m % Wiederholung für AVGTrials, mitteln von 10er Gruppen;
LR_test_ind_withoutNoise.m % U-Test für verschiedene Zeitintervalle; Linke vs. rechte Hemisphere, hier nohc keine ROIextraction
UTest_normalize.m % NOrmalisierung in MNI space
extract_UValuesROI.m % Berechnung des LI-Wertes für ROIs
kh_collect_LI % sammeln der Werte für alle Subjects in einer Matrix



% patients without noise estimation:

kh_3dtagalign.m
kh_MRI2Afni.m
SAM_patients.m
kh_SAM_Beamforming_keeptrials.m
patients_avg.m
Timeintervall_patients.m % wird nur für Infofile benötigt, um später Brik zu speichern
LR_UTest_ind_UValues.m % funktioniert nur mit Infofile aus Timeintervall_patients.m
SAM_normalize.m % normalisiert ortho
UTest_normalize.m % NOrmalisierung in MNI space der UTest-Ergebnisse
extract_UValuesROI_patients.m % Berechnung des LI-Wertes für ROIs
kh_collect_LI % sammeln der Werte für alle Subjects in einer Matrix

Extract_Max_Activity_no_noise_voxelcomparison.m % Y. meinte in Mail, dass es helfen könnte, keine noisenormalisierung durchzuführen und dann LI auf Voxelbasis zu berechnen


% außerdem für Diss:
extract_Act_patients2.m % erstellt Abbildungen aus maximaler Activität pro Sample für linke und rechte ROI
Extract_Max_Activity.m % extrahiert Maxima in ROI und berechnet LIs für Summe aus Maxima pro sample und aus höchsten Maximum aller samples



% patients with noise normalization
patients_avg.m
Timeintervall_patients.m % mit z-transformation und normalisierung


% 
LI_zsores.m % für Patienten
LI_zsores_controls.m % für Kontrollen, Achtung: Summe wurde durch Anzahl dividiert, funktioniert nicht, noch rückgängig machen! (siehe Patienten)
FindExtrema.m % Findet pro Maske Extremwerte und gibt Atlas-Region aus => gesamtes Gehirn
FindExtrema_Controls.m  % Findet pro Maske Extremwerte und gibt Atlas-Region aus  => gesamtes Gehirn

%
Timeintervall_test_abs.m % gleicht Fehler aus (sum(abs)) statt abs(sum)
M400.m % extrahiert nur diesen Zeitbereich
Li_z_Scores_405ms.m % berechnet LI für 405ms für Kontrollen, hat allerdings nicht funktioniert (Li's teilweise rechtslateral)

% fMRI Fluency VG comparison:

FindExtrema_fMRI_Fluency.m % complete Brain
FindExtrema_fMRI_FluencyVG_comp_ROIs.m % only ROIS
FindClust_fMRI_Fluency.m %finds Clusters with peak activity for Broca and Wernicke 

% nützlich:

merge_matrices.m % für vektoren
ForExcel.m

%MEG-fMRI-comparison:
FindExtrema_fMRI.m % finds extrema in Verbgeneration task, converting wspmT_001 to z-normalisation and downsampling to MEG resolution
FindExtrema_fMRI_downs_ROI.m % Findet pro Maske Extremwerte und gibt Atlas-Region aus =>> nur ROI downsampled
FindExtrema_MEG_ROI_patients.m % Findet pro Maske Extremwerte und gibt Atlas-Region aus =>> nur ROI
IPN_ccc % Pearson R and CCC for MEG and fMRI Verbgeneration
Percent_Overlap_MEG_abs.m %inklusive fMRI informed MEG => noch ausbauen; % außerdem: Möglichkeit, neue Maske zu erstellen mit Überlappungsbereich
merge_Clusters.m

FindClust_VG.m % Findet pro Maske Clusterwerte  

% TTest für jeweils 10 gemittelte Trials:
LR_TTest_AVG10_noise_abs.m % noch unvollständig, Prüfung an einem Probanden führte zu keinen wesentlichen Änderungen gegenüber dem U-Test

% Versuch, ttest linke gegen rechte hemisphäre zu berechnen bei ungemittelten Daten:
CopyCleanData.m % kopiert CleanData.mat von Windows zum Linuxordner 
LR_ttest_singleTrials.m  % => hier auch Beschreibung für Patienten, welche gesäuberten Daten genommen worden sidn

% zum ausprobieren:
kh_reshape.m: % aus 2D => 4D ohne abspeichern
