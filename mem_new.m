% set middle-ear mass
me.smem(iep,iep) = pa.mdi / pa.adi^2;
me.smem(ied,ied) = pa.mcp / pa.acp^2 + pa.med / pa.aed^2;
me.smem(ied,ima) = -pa.med / pa.aed^2;
me.smem(ima,ied) = -pa.med / pa.aed^2;
me.smem(ima,ima) = pa.mma / pa.ama^2 + pa.med / pa.aed^2;
me.smem(ist,ist) = pa.mst;
% set middle-ear damping
me.smer(iep,iep) = (pa.rdi + rep + ra) / pa.adi^2;
me.smer(iep,ied) = -ra / (pa.acp * pa.adi);
me.smer(ied,iep) = -ra / (pa.acp * pa.aed);
me.smer(ied,ied) = (pa.rcp + ra) / pa.acp^2 + pa.red / pa.aed^2;
me.smer(ied,ima) = -pa.red / pa.aed^2;
me.smer(ima,ied) = -pa.red / pa.aed^2;
me.smer(ima,ima) = (pa.rma +  gme^2 * pa.rim) / pa.ama^2 + pa.red / pa.aed^2;
me.smer(ima,ist) = -gme * pa.rim / (pa.ama * pa.ast);
me.smer(ist,ima) = -gme * pa.rim / (pa.ama * pa.ast);
me.smer(ist,ist) = (pa.rst + pa.rim) / pa.ast^2;
% set middle-ear stiffness
me.smek(iep,iep) = (pa.kdi + ka) / pa.adi^2;
me.smek(iep,ied) = -ka / (pa.acp * pa.adi);
me.smek(ied,iep) = -ka / (pa.acp * pa.adi);
me.smek(ied,ied) = ka / pa.acp^2 + pa.red / pa.aed^2;
me.smek(ied,ima) = -pa.ked / pa.aed^2;
me.smek(ima,ied) = -pa.ked / pa.aed^2;
me.smek(ima,ima) = (pa.kma +  gme^2 * pa.kim) / pa.ama^2;
me.smek(ima,ist) = -(gme * pa.kim) / (pa.ama * pa.ast);
me.smek(ist,ima) = -(gme * pa.kim) / (pa.ama * pa.ast);
me.smek(ist,ist) = pa.kst + pa.kim / pa.ast^2;
%
me.smer(1,:)=me.smer(1,:)/me.ame(1);
me.smer(2,:)=me.smer(2,:)/me.ame(2);
me.smer(3,:)=me.smer(3,:)/me.ame(3);
me.smer(4,:)=me.smer(4,:)/me.ame(4);
% stapes BC
me.alfx = (pa.mst / pa.ast^2) / (2 * pa.rho * dx);
