PROJ='/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('refit_m4_full.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
H=parfit26('handles');
% which fields must be pinned? the guard names them
b=modpar26(4); lost=H.lost(pa,H.setpar(b,H.getpar(pa),[]),[]);
fprintf('\n  unpinned, fields lost: %s\n', strjoin(lost,', '));
b2=modpar26(4); for i=1:numel(lost), b2.(lost{i})=pa.(lost{i}); end
lost2=H.lost(pa,H.setpar(b2,H.getpar(pa),[]),[]);
fprintf('  with those pinned    : %s\n', ternx(isempty(lost2),'(none) OK',strjoin(lost2,', ')));
fprintf('  chsz = [%s] (%d entries; fitidx 31..%d available)\n', ...
        strtrim(sprintf('%.3f ',pa.chsz)), numel(pa.chsz), 30+numel(pa.chsz));
t=tic; m=abr_metric(pa,false); t1=toc(t);
t=tic; evalc('tdm26(''coupeig'',struct(''pa'',pa));'); t2=toc(t);
t=tic; S=score26(pa,'fast',false); t3=toc(t);
ev=t1+t2+t3;
fprintf('  eval cost: abr %.1f + coupeig %.1f + score26 %.1f = %.1f s\n', t1,t2,t3,ev);
fprintf('  6 h -> %d evals | 7 h -> %d evals\n', floor(6*3600/ev), floor(7*3600/ev));
fprintf('  start: maperr %.2f amp %+.2f osc %+.2f | shoulder %.4f (%d/16>0.5)\n', ...
        S.maperr,S.amp_gain,S.maxRe_osc,m.shoulder,sum(m.sho(isfinite(m.sho))>0.5));
disp('M4PRE_DONE');
function s=ternx(c,a,b), if c, s=a; else, s=b; end, end
