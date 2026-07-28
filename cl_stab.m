% Is the gain>=10 blow-up PHYSICAL or NUMERICAL?
%   coupeig Re>0  => coupled model genuinely unstable; no integrator can help
%   coupeig Re<0 but explicit scheme diverges => NUMERICAL; implicit would help
% This decides whether "no stable amplifying regime" is a fact about the
% mechanism or an artifact of the explicit time-stepping.
k5s=0.3; base=modpar26(4);
fprintf('\n  gain |  coupeig maxRe_osc   verdict\n');
fprintf('%s\n',repmat('-',1,50));
for g=[0 1 3 10 30]
    pa=modpar26(4); pa.ohcsgn=+1; pa.ohcgain=g; pa.k5o=base.k5o*k5s;
    try
        evalc('S=tdm26(''coupeig'',struct(''pa'',pa));');
        v='sub-critical (stable)'; if (S.maxRe_osc>=0), v='UNSTABLE (physical)'; end
        fprintf('  %4.0f |  %+14.1f   %s\n', g, S.maxRe_osc, v);
    catch e
        fprintf('  %4.0f |  FAILED: %s\n', g, e.message);
    end
end
fprintf('\nIf gain=10 is sub-critical here but diverged in the time domain,\n');
fprintf('the limit is the explicit integrator, NOT the mechanism.\n');
disp('CL_STAB_DONE');
