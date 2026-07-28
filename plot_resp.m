function plot_resp(stcfg)
tdm24(stcfg,-1); prneps(1);    movefile('prneps01.eps', 'aba-ft.eps');
tdm24(stcfg,1); prneps(1);     movefile('prneps01.eps', 'aba-nr.eps');
tdm24(stcfg,1,1,3); prneps(1); movefile('prneps01.eps', 'aba-ld.eps');
end
