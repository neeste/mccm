for fn={'tdm26.m','fdm26.m','parfit26.m'}
  e=checkcode(fn{1},'-string');
  err=regexp(e,'\n','split'); nerr=sum(~cellfun(@isempty,regexp(err,'(Error|Parse error|Invalid)','once')));
  fprintf('%-12s : %d serious mlint issues\n', fn{1}, nerr);
end
disp('INTEGRITY_DONE');
