function plot_aba
tic; tdm24('aba.wav',-1);    toc;print(figure(1),'-dpng', 'aba-ft.png');
tic; tdm24('aba.wav',1,1,3); toc;print(figure(1),'-dpng', 'aba-ld.png');
end
