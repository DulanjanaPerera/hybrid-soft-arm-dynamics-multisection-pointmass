function validateRlReferences()
%VALIDATERLREFERENCES Offline schedule, continuity and cycle checks.
functions={@rlCircleSweepReferenceBlock,@rlHelixReferenceBlock};
for k=1:2
    f=functions{k};[th,ph,~,~,~,~,done,finish]=f(0);
    assert(th==0 && ph==0 && ~done);
    previous=[th;ph];maxThetaRate=0;
    for t=.05:.05:finish+.1
        [th,ph,~,phase,~,direction,done]=f(t);
        assert(isfinite(th) && isfinite(ph) && ph>=-1e-12 && ph<=1.5+1e-12);
        dth=(th-previous(1))/.05;maxThetaRate=max(maxThetaRate,abs(dth));
        assert(abs(dth)<=2*pi/20+1e-8);
        assert(abs(ph-previous(2))/.05<=.451);
        if phase==5 || phase==6 || phase==7
            assert(abs(direction)==1);
        end
        if done,assert(ph==0);end
        previous=[th;ph];
    end
    [~,ph,~,phase,~,~,done]=f(finish);
    assert(done && ph==0 && phase==9);
    [th,ph,~,phase]=f(NaN);assert(th==0 && ph==0 && phase==-1);
    fprintf('%s: duration %.2fs, max theta speed %.6f rad/s\n',func2str(f),finish,maxThetaRate);
end
levels=[.2 .8 1 1.2 1.5];
for j=1:5
    start=5+(j-1)*123;
    [~,ph,trial,phase]=rlCircleSweepReferenceBlock(start+6);
    assert(trial==j && phase==2 && abs(ph-levels(j))<1e-12);
    for lap=1:5
        t=start+8+(lap-1)*20+1;
        [th,ph,trial,phase,cycle]=rlCircleSweepReferenceBlock(t);
        [th2,~,~,phase2]=rlCircleSweepReferenceBlock(t+.05);
        assert(trial==j && phase==5 && phase2==5 && cycle==lap);
        assert(abs(ph-levels(j))<1e-12 && abs((th2-th)/.05-2*pi/20)<1e-9);
    end
end
for j=1:4
    start=5+(j-1)*130;
    [~,ph,trial,phase,~,direction]=rlHelixReferenceBlock(start+60);
    assert(trial==j && abs(ph-1.5)<1e-12 && phase==7 && direction==1-2*mod(j-1,2));
    [~,ph,~,phase]=rlHelixReferenceBlock(start+125);assert(ph==0 && phase==4);
end
disp('PASS: circle sweep and loading/unloading helix references.');
end
