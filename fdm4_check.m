% Does fdm26's DEFAULT path (the one parfit26 uses for maperr) work for m=4?
for nch=[3 4]
    try
        R=fdm26(struct('pa',modpar26(nch)));
        f={}; if(isstruct(R)), f=fieldnames(R); end
        fprintf('nch=%d: OK, maperr=%s\n', nch, ...
            char(string(  (isstruct(R)&&isfield(R,'maperr'))*1  )));
        if (isstruct(R)&&isfield(R,'maperr')), fprintf('   maperr=%.2f\n', R.maperr); end
    catch e
        fprintf('nch=%d: THREW -> %s\n', nch, e.message);
    end
end
disp('FDM4_CHECK_DONE');
