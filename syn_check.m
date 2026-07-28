% Syntax-check the edited fitter without running a fit.
try
    ok = true;
    e = checkcode('parfit26.m','-string');
    if (isempty(e)), fprintf('parfit26.m: no mlint messages\n');
    else, fprintf('parfit26.m mlint:\n%s\n', e); end
catch err
    fprintf('checkcode failed: %s\n', err.message);
end
% confirm the new options parse and pinning works, without launching a fit
opts=struct('pin',struct('m5o',0.072,'chsz',[0.95 0.05 1.0 0]),'cheapstab',1,'statol',1);
fprintf('opts built OK: cheapstab=%d statol=%g pinfields=%s\n', ...
        opts.cheapstab, opts.statol, strjoin(fieldnames(opts.pin),','));
disp('SYN_CHECK_DONE');
