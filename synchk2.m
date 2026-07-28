e=checkcode('parfit26.m','-string');
n=sum(~cellfun(@isempty,regexp(regexp(e,'\n','split'),'(Error|Parse error|Invalid)','once')));
fprintf('parfit26.m serious mlint issues: %d\n', n);
e2=checkcode('launch_c1shape.m','-string');
n2=sum(~cellfun(@isempty,regexp(regexp(e2,'\n','split'),'(Error|Parse error|Invalid)','once')));
fprintf('launch_c1shape.m serious mlint issues: %d\n', n2);
disp('SYNCHK2_DONE');
