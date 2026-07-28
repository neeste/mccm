% Is the CL compartment simply too small?  d1 pushes ST<->CL; with chsz(CL)=0.05
% the BM works into a tiny, stiff volume.  Sweep the CL size and watch whether
% the partition displacements (and hence the compression) recover.
% 3-chamber reference @2kHz/40dB: d1=2.58e-08, d2=4.73e-08, d2/hbmx=7.89
fr=2; lvl=40;
fprintf('\n chsz(CL) |    max|d1|     max|d2|   d2/d1   d2/hbmx   lat(ms)\n');
fprintf('%s\n',repmat('-',1,64));
for c4=[0.05 0.2 0.5 1.0 2.0]
    p.fr=fr; p.lv=lvl; p.pa=modpar26(4); p.pa.hbmode='bm';
    p.pa.chsz=[0.95 0.05 1.0 c4];
    try
        S=tdm26('wnr1',p,0,0); d=S.dgn;
        fprintf('%9.2f | %11.3e %11.3e %7.3f %9.2f %9.2f\n', ...
                c4, d.d1mx, d.d2mx, d.d2mx/d.d1mx, d.ratio, S.tpk);
    catch e
        fprintf('%9.2f | FAILED: %s\n', c4, e.message);
    end
end
fprintf('\n3-chamber ref: d1=2.581e-08 d2=4.733e-08 d2/d1=1.83 d2/hbmx=7.89\n');
disp('CHSZ4_SWEEP_DONE');
