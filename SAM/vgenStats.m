% Compute Diff between abs left and right values
[~,Li]=ismember(LRpairs(:,1),ga.label); % ga ) grandaverage with individually kept trials, LRpairs is file from Yuval
[~,Ri]=ismember(LRpairs(:,2),ga.label);
gaLR=ga;
gaLR.individual=zeros(size(gaLR.individual)); % die mittleren Sensoren behalten Nullen
gaLR.individual(:,Li,:)=abs(ga.individual(:,Li,:))-abs(ga.individual(:,Ri,:));
gaLR.individual(:,Ri,:)=abs(ga.individual(:,Ri,:))-abs(ga.individual(:,Li,:));

aliceTtest0(gaLR, 0.4, 1);
vgenTtest0(gaLR, 0.4, 1,0.05,ga); % Input Nr. 3 = p-value
vgenTtest0(gaLR, 0.17, 1,0.05,ga);
vgenTtest0(gaLR, 0.3, 1,0.05,ga);


% compute RMS
rmsL=squeeze(sqrt(mean(ga.individual(:,Li,:).^2,2)));
rmsR=squeeze(sqrt(mean(ga.individual(:,Ri,:).^2,2)));

t=0.4;[~,p] = ttest(rmsL(:,nearest(ga.time,t)),rmsR(:,nearest(ga.time,t)));

figure;plot(ga.time,mean(rmsL),'r')
hold on
plot(ga.time,mean(rmsR))
legend('L','R')

t1=0.32;s1=nearest(ga.time,t1);
t2=0.47;s2=nearest(ga.time,t2);
areaL=trapz(ga.time(s1:s2),rmsL(:,s1:s2)');
areaR=trapz(ga.time(s1:s2),rmsR(:,s1:s2)');
[~,p]=ttest(areaL,areaR)