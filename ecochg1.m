% ECochG protocol #1 - simulate extra-cochlear ECochG "click" response
function ecochg1
pr.fr=4;        % tone frequency (kHz)
pr.lv=100;      % tone level (dB SPL)
pr.td=0.6;      % stimulus duration (ms)
pr.tr=0.1;      % ramp duration (ms)
pr.ts=0.2;      % stimulus start (ms)
pr.tp=5;        % plotted time (ms)
pr.sd=0.25;     % synaptic delay (ms)
pr.xe=2;        % electrode positions (mm)
pr.wvc=1;       % volume-conduction width (mm)
pr.hca=0.999;   % high-pass-filter coefficient
pr.hcc=80;      % hair-cell currrent
pr.nrc=1.8e-5;  % neural-rate currrent
pr.sumlim=[-3 1];      % summed range
pr.diflim=pr.sumlim/2; % difference range
diary ecochg.log
tdm25('ecochg',pr);
diary off
print('-dpng',1,'','ecochg1.png');
print('-dpng',2,'','ecochg1a.png');
end