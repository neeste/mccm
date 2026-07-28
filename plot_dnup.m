function plot_dnup
tic; tdm24('swpdn.wav',-1);    toc;print(figure(1),'-dpng', 'swpdn-ft.png');
tic; tdm24('swpdn.wav',1,1,3); toc;print(figure(1),'-dpng', 'swpdn-ld.png');
tic; tdm24('swpup.wav',-1);    toc;print(figure(1),'-dpng', 'swpup-ft.png');
tic; tdm24('swpup.wav',1,1,3); toc;print(figure(1),'-dpng', 'swpup-ld.png');
end
