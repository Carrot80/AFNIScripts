plot(trig)
trig=readTrig_BIU;
plot(trig)
unique(trig)
trig4224=bitand(uint16(trig),4224);
hold
on
hold on
plot(trig4224,'r')
trig=bitand(uint16(trig),4224);
trigSh=zeros(size(trig));
trigSh(2:end)=trig(1:end-1);
trigOnset=find((trig-trigSh)>0);
trig=double(bitand(uint16(trig),4224));
trigSh=zeros(size(trig));
trigSh(2:end)=trig(1:end-1);
trigOnset=find((trig-trigSh)>0);
plot(trig,'r')
hold on
plot(trigSh,'g')
plot(trigOnset,trig(trigOnset),'.k')
plot(trig,'r')
hold on
plot(trigSh,'g')
plot(trigOnset,trig(trigOnset),'.k')