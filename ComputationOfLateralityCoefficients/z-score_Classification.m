function ()

forall('zscores_LI_RelAct_0.32_0.47s_Controls.mat')
end


function forall(FileName) 
load('zscores_LI_RelAct_0.32_0.47s_Controls.mat')
x=LI_ALL_zscores
Broca=x([1:24],[1,9, 17])
Wernicke=x([1:24],[5,13,21])



end

function leftlat()

for i=1:size(Broca,2)
Broca_15_right(i)=sum(Broca(:,i)<=-.15)
Broca_15_left(i)=sum(Broca(:,i)>=.15)
Broca_15_bil(i)=sum(Broca(:,i)<.15 & Broca(:,i)>-.15)

Broca_all(i,1:3)=[Broca_15_left(i) Broca_15_right(i) Broca_15_bil(i)]
end

for i=1:size(Broca,2)
Broca_1_right(i)=sum(Broca(:,i)<=-.1)
Broca_1_left(i)=sum(Broca(:,i)>=.1)
Broca_1_bil(i)=sum(Broca(:,i)<.1 & Broca(:,i)>-.1)

Broca_all(i,1:3)=[Broca_1_left(i) Broca_1_right(i) Broca_1_bil(i)]
end


for i=1:size(Broca,2)
Broca_2_right(i)=sum(Broca(:,i)<=-.2)
Broca_2_left(i)=sum(Broca(:,i)>=.2)
Broca_2_bil(i)=sum(Broca(:,i)<.2 & Broca(:,i)>-.2)

Broca_all(i,1:3)=[Broca_2_left(i) Broca_2_right(i) Broca_2_bil(i)]
end



for i=1:size(Wernicke,2)
Wernicke_15_right(i)=sum(Wernicke(:,i)<=-.15)
Wernicke_15_left(i)=sum(Wernicke(:,i)>=.15)
Wernicke_15_bil(i)=sum(Wernicke(:,i)<.15 & Wernicke(:,i)>-.15)

Wernicke_all(i,1:3)=[Wernicke_15_left(i) Wernicke_15_right(i) Wernicke_15_bil(i)]
end

for i=1:size(Wernicke,2)
Wernicke_1_right(i)=sum(Wernicke(:,i)<=-.1)
Wernicke_1_left(i)=sum(Wernicke(:,i)>=.1)
Wernicke_1_bil(i)=sum(Wernicke(:,i)<.1 & Wernicke(:,i)>-.1)

Wernicke_all(i,1:3)=[Wernicke_1_left(i) Wernicke_1_right(i) Wernicke_1_bil(i)]
end


for i=1:size(Wernicke,2)
Wernicke_2_right(i)=sum(Wernicke(:,i)<=-.2)
Wernicke_2_left(i)=sum(Wernicke(:,i)>=.2)
Wernicke_2_bil(i)=sum(Wernicke(:,i)<.2 & Wernicke(:,i)>-.2)

Wernicke_all(i,1:3)=[Wernicke_2_left(i) Wernicke_2_right(i) Wernicke_2_bil(i)]
end

end