f = linspace(500,4000,2000);
sav1=tdm26(0,1); pa1=sav1.pa; pa1.gam=1; [~,~,~,~,~,Zc1,~]=fdmod23(pa1,f); Ze1=midear_imp(pa1,f,Zc1);
sav3=tdm26(0,3); pa3=sav3.pa; pa3.gam=1; [~,~,~,~,~,Zc3,~]=fdmod23(pa3,f); Ze3=midear_imp(pa3,f,Zc3);
db1 = 20*log10(abs(Ze1)); rip1 = db1 - movmean(db1,50);
db3 = 20*log10(abs(Ze3)); rip3 = db3 - movmean(db3,50);
fprintf('Max Ripple m=1: %g dB\n', max(abs(rip1)));
fprintf('Max Ripple m=3: %g dB\n', max(abs(rip3)));
