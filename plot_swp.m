function plot_swp
tic; tdm24('swp1.wav',-1);    toc;print(figure(1),'-dpng', 'swp1-ft.png');
tic; tdm24('swp1.wav',1,1,3); toc;print(figure(1),'-dpng', 'swp1-ld.png');
tic; tdm24('swp2.wav',-1);    toc;print(figure(1),'-dpng', 'swp2-ft.png');
tic; tdm24('swp2.wav',1,1,3); toc;print(figure(1),'-dpng', 'swp2-ld.png');
end
