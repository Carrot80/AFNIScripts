
Fluency_Pat=load(strcat('D:\Arbeit\LinuxExchange\data\Cluster_and_Locations_fMRI_t-values_', 'Patients','_Fluency_ROIs.mat'));
VG_Pat=load(strcat('D:\Arbeit\LinuxExchange\data\Cluster_and_Locations_fMRI_t-values_', 'Patients','_Verbgeneration_ROIs.mat'));
Fluency_Con=load(strcat('D:\Arbeit\LinuxExchange\data\Cluster_and_Locations_fMRI_t-values_', 'Controls','_Fluency_ROIs.mat'));
VG_Con=load(strcat('D:\Arbeit\LinuxExchange\data\Cluster_and_Locations_fMRI_t-values_', 'Controls','_Verbgeneration_ROIs.mat'));

t_VG_Broca=VG_Pat.Table.Broca.Coord(:,13)
t_VG_Broca(25:34)=VG_Con.Table.Broca.Coord(:,13)

t_VG_Wernicke=VG_Pat.Table.Wernicke.Coord(:,13)
t_VG_Wernicke(25:34)=VG_Con.Table.Wernicke.Coord(:,13)

t_Fluency_Broca=VG_Pat.Table.Broca.Coord(:,13)
t_Fluency_Broca(25:34)=VG_Con.Table.Broca.Coord(:,13)

t_Fluency_Wernicke=VG_Pat.Table.Wernicke.Coord(:,13)
t_Fluency_Wernicke(25:34)=VG_Con.Table.Wernicke.Coord(:,13)

t_values(:,1)=t_VG_Broca
t_values(:,2)=t_VG_Wernicke
t_values(:,3)=t_Fluency_Broca
t_values(:,4)=t_Fluency_Wernicke

xlswrite('t_values', t_values)


%% MEG fMRI comparison



Pat=load(strcat('D:\Arbeit\LinuxExchange\Results\Comparison_fMRI_MEG\Clusters_and_Locations_', 'Patients','_fMRI_MEG_comparison.mat'));
Con=load(strcat('D:\Arbeit\LinuxExchange\Results\Comparison_fMRI_MEG\Clusters_and_Locations_', 'Controls','_fMRI_MEG_comparison.mat'));

z_MEG_Broca=Pat.Table.MEG.Broca.Coord(:,13)
z_MEG_Broca(25:34)=Con.Table.MEG.Broca.Coord(:,13)

z_MEG_Wernicke=Pat.Table.MEG.Wernicke.Coord(:,13)
z_MEG_Wernicke(25:34)=Con.Table.MEG.Wernicke.Coord(:,13)

z_fMRI_Broca=Pat.Table.fMRI.Broca.Coord(:,13)
z_fMRI_Broca(25:34)=Con.Table.fMRI.Broca.Coord(:,13)

z_fMRI_Wernicke=Pat.Table.fMRI.Wernicke.Coord(:,13)
z_fMRI_Wernicke(25:34)=Con.Table.fMRI.Wernicke.Coord(:,13)

z_values(:,1)=z_MEG_Broca
z_values(:,2)=z_MEG_Wernicke
z_values(:,3)=z_fMRI_Broca
z_values(:,4)=z_fMRI_Wernicke

xlswrite('z_values', z_values)


ccc=[2*corr(z_MEG_Broca,z_fMRI_Broca).*std(z_MEG_Broca).*std(z_fMRI_Broca)]./[std(z_MEG_Broca)^2+std(z_fMRI_Broca)^2+(mean(z_MEG_Broca)-mean(z_fMRI_Broca))^2 ];

ccc_Wernicke=[2*corr(z_MEG_Wernicke,z_fMRI_Wernicke).*std(z_MEG_Wernicke).*std(z_fMRI_Wernicke)]./[std(z_MEG_Wernicke)^2+std(z_fMRI_Wernicke)^2+(mean(z_MEG_Wernicke)-mean(z_fMRI_Wernicke))^2 ];

distance(:,1)=d.Broca
distance(:,2)=d.Wernicke
xlswrite('distance', distance)


d_Broca_merged = zeros(34,1)

j=1

for i=1:34
   
    d_Broca_merged(j,:)=d.Broca(i, :);
    j=j+2;
end

xlswrite('d_Broca', d_Broca_merged)

 
d_Wernicke_merged = zeros(34,1)

j=1

for i=1:34
   
    d_Wernicke_merged(j,:)=d.Wernicke(i, :);
    j=j+2;
end

xlswrite('d_Wernicke', d_Wernicke_merged)

