function validateRlPressureBlocks()
%VALIDATERLPRESSUREBLOCKS Offline waveform checks only; never invokes hardware.
functions={@rlPressureStepsBlock,@rlPressureSineBlock,@rlPressureCombinedBlock};
for k=1:numel(functions)
    f=functions{k};[p,~,~,~,~,done,finish]=f(0);
    assert(all(p==0) && ~done && finish>0);
    previous=p;maxRate=0;
    for t=.05:.05:finish+.1
        [p,trial,phase,channel,cycle,done]=f(t);
        assert(all(isfinite(p)) && all(p>=-1e-12 & p<=3+1e-12));
        assert(trial>=0 && phase>=0 && channel>=0 && channel<=3 && cycle>=0);
        maxRate=max(maxRate,max(abs(p-previous))/.05);previous=p;
        if k<3
            assert(nnz(abs(p)>1e-12)<=1);
            if channel>0,other=[1 2 3]~=channel;assert(all(p(other)==0));end
        end
        if done,assert(all(p==0));end
    end
    [p,~,phase,~,~,done]=f(finish);
    assert(done && phase==9 && all(p==0));
    [p,~,phase]=f(NaN);assert(phase==-1 && all(p==0));
    fprintf('%s: duration %.3fs, max sampled request slew %.4f bar/s\n',func2str(f),finish,maxRate);
    assert(maxRate<1.0); % Default waveforms are deliberately smooth/slow.
end
% Raw sine matches the previous single-P1 test after the initial baseline.
for tau=[0 1 5 10 15 20 30]
    p=rlPressureSineBlock(5+tau);
    assert(abs(p(1)-1.5*(1-cos(.2*tau)))<1e-12);
end
disp('PASS: RL open-loop pressure source blocks; no hardware/model run.');
end
