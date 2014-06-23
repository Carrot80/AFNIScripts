
Files = 'D:\Arbeit\LinuxExchange\data\patients\patients_SAM';
Dir = dir(Files);
isub = [Dir(:).isdir];
nameFolds = {Dir(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];


for i= 1:length(nameFolds)
LI_patients.Broca(i,1) = LI_All.(nameFolds{i}).Broca.Classic;
end


for i= 1:length(nameFolds)
LI_patients.Broca(i,2) = LI_All.(nameFolds{i}).Broca.Max_percchange;
end

for i= 1:length(nameFolds)
LI_patients.Broca(i,3) = LI_All.(nameFolds{i}).Broca.signLeftVoxels;
end


for i= 1:length(nameFolds)
LI_patients.Broca(i,4) = LI_All.(nameFolds{i}).Broca.signRightVoxels;
end

for i= 1:length(nameFolds)
LI_patients.Broca(i,5) = LI_All.(nameFolds{i}).Broca.SizeROI;
end

for i= 1:length(nameFolds)
LI_patients.Broca(i,6) = LI_All.(nameFolds{i}).Broca.relActLeft;
end


for i= 1:length(nameFolds)
LI_patients.Broca(i,7) = LI_All.(nameFolds{i}).Broca.relActRight;

end

%%


for i= 1:length(nameFolds)
LI_patients.Broca(i,1) = LI_All.(nameFolds{i}).Broca.Classic;
end


for i= 1:length(nameFolds)
LI_patients.Wernicke(i,2) = LI_All.(nameFolds{i}).Wernicke.Max_percchange;
end

for i= 1:length(nameFolds)
LI_patients.Wernicke(i,3) = LI_All.(nameFolds{i}).Wernicke.signLeftVoxels;
end


for i= 1:length(nameFolds)
LI_patients.Wernicke(i,4) = LI_All.(nameFolds{i}).Wernicke.signRightVoxels;
end

for i= 1:length(nameFolds)
LI_patients.Wernicke(i,5) = LI_All.(nameFolds{i}).Wernicke.SizeROI;
end

for i= 1:length(nameFolds)
LI_patients.Wernicke(i,6) = LI_All.(nameFolds{i}).Wernicke.relActLeft;
end


for i= 1:length(nameFolds)
LI_patients.Wernicke(i,7) = LI_All.(nameFolds{i}).Wernicke.relActRight;

end

