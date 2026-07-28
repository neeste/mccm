for f={'tdm26.m','modpar26c3b.m','m3b_gate.m'}
  e=checkcode(f{1},'-string');
  n=sum(~cellfun(@isempty,regexp(regexp(e,'\n','split'),'(Error|Parse error|Invalid)','once')));
  fprintf('%-16s serious mlint issues: %d\n', f{1}, n);
end
disp('SYNCHK3_DONE');
