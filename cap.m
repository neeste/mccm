global DBG_W; DBG_W={};
tdm26('tbabr',1,0,0); n1=numel(DBG_W);
tdm26('tbabr',3,0,0);
save('traces.mat','DBG_W','n1');
fprintf('captured %d traces (nch1=%d, nch3=%d)\n', numel(DBG_W), n1, numel(DBG_W)-n1);
