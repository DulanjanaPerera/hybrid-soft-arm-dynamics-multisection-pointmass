function validateResidualScreening(dataset,screen)
%VALIDATERESIDUALSCREENING Check fitted statistics and target-free inference.
p=load(fullfile(screen,'protocol.mat'));assert(strcmp(p.cfg.status,'complete'));
constraints=table;solver=table;
for task={'forward','inverse'}
    task=task{1};f=load(fullfile(screen,task,'frozen_selection.mat'));
    T=residualBenchmarkData(dataset,task,'train',p.cfg.trainingStride);
    V=residualBenchmarkData(dataset,task,'validation',1);
    assert(all(startsWith(T.Run,'day1_'))&&all(startsWith(V.Run,'day1_')));
    for i=1:numel(f.bank)
        m=f.bank(i);assert(max(abs(m.contextMean-mean(T.Context)))<1e-10);
        assert(max(abs(m.targetMean-mean(T.Y)))<1e-10);
        % A small inference subset; target arrays cannot influence predictions.
        D=V;D.Seq=D.Seq(:,:,1:4);D.Context=D.Context(1:4,:);D.Base=D.Base(1:4,:);
        D.Y=D.Y(1:4,:);D.Truth=D.Truth(1:4,:);D.p0=D.p0(1:4,:);
        a=residualPredict(m,D);D.Y(:)=1e8;D.Truth(:)=-1e8;
        b=residualPredict(m,D);assert(isequaln(a,b)&&all(isfinite(a(:))));
    end
    m=f.bank(1);assert(isa(m.W,'double'));
    X=[double(reshape(T.Seq,[],size(T.Seq,3))') (T.Context-m.contextMean)./m.contextStd];
    A=[X ones(size(X,1),1)];P=eye(size(A,2));P(end,end)=0;
    Z=double(single((T.Y-m.targetMean)./m.targetStd));dist=[];
    for lambda=p.cfg.ridgeLambda
        W=[A;sqrt(lambda*size(A,1))*P]\[Z;zeros(size(A,2),2)];dist(end+1)=norm(W-m.W,'fro'); %#ok<AGROW>
    end
    [error,j]=min(dist);assert(error<1e-8);
    solver=[solver;table(string(task),p.cfg.ridgeLambda(j),error,'VariableNames', ...
        {'task','selectedLambda','reproductionError'})]; %#ok<AGROW>
    if strcmp(task,'inverse')
        for partition={'validation','external_test'}
            D=residualBenchmarkData(dataset,task,partition{1},1);B=[-1/3 -1/3;2/3 -1/3;-1/3 2/3];
            for i=1:numel(f.bank)
                R=residualPredict(f.bank(i),D);pred=D.p0+(D.Base+R)*B';
                gap=max(max(-pred,pred-3),0);
                row=table(string(partition{1}),string(f.bank(i).name),max(gap,[],'all'), ...
                    mean(any(gap>.01,2)),'VariableNames',{'partition','model', ...
                    'maximumBoundViolation_bar','violationAbovePoint01Fraction'});
                constraints=[constraints;row]; %#ok<AGROW>
            end
        end
    end
end
writetable(solver,fullfile(screen,'ridge_reproduction.csv'));
writetable(constraints,fullfile(screen,'inverse_constraints.csv'));
fprintf('PASS: train-only fitted statistics, target-free frozen inference, exact double-QR ridge reproduction and unclipped pressure-constraint audit.\n');
end
