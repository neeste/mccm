L=load('traces.mat'); D=L.DBG_W; frq=[0.5 1 2 4]; lvl=[20 40 60 80];
det=@(w,dt,b0) det_impl(w,dt,b0);
labels={'nch=1','nch=3'};
for blk=1:2
    off=(blk-1)*16;
    fprintf('\n==== %s ====  tpk for baseline = {2%%mean, p10, p25, median}, thr=3x\n', labels{blk});
    fprintf('%4s %3s %8s | %7s %7s %7s %7s\n','fr','lv','global','b2pct','p10','p25','med');
    for i=1:16
        d=D{off+i}; w=d.w; dt=d.dtms; n=numel(w);
        fr=frq(ceil(i/4)); lv=lvl(mod(i-1,4)+1);
        [~,ig]=max(w); tg=(ig-1)*dt;
        b2=mean(w(1:max(1,round(0.02*n))));
        b10=prctile(w,10); b25=prctile(w,25); bmed=median(w);
        fprintf('%4.1f %3d %8.2f | %7.2f %7.2f %7.2f %7.2f\n', fr,lv,tg, ...
            det(w,dt,b2),det(w,dt,b10),det(w,dt,b25),det(w,dt,bmed));
    end
end
function tpk=det_impl(w,dt,b0)
n=numel(w); j=find(w>=3*b0,1,'first');
if isempty(j), tpk=NaN; return; end
while (j<n && w(j+1)>=w(j)), j=j+1; end
tpk=(j-1)*dt;
end
