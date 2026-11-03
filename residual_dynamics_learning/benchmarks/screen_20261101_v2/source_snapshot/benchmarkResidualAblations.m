function benchmarkResidualAblations(dataset,screen)
%BENCHMARKRESIDUALABLATIONS Validation-motivated shortcut checks, fixed protocol.
folder=fullfile(screen,'ablations');assert(~isfolder(folder));mkdir(folder);
protocol='Ridge grid unchanged; inverse pressure-only/motion-only and forward kinematic residual; validation-only fitting/selection.';
save(fullfile(folder,'protocol.mat'),'protocol');
for task={'inverse','forward'}
    task=task{1};d=load(fullfile(screen,task,'frozen_selection.mat'));
    T=residualBenchmarkData(dataset,task,'train',d.cfg.trainingStride);
    V=residualBenchmarkData(dataset,task,'validation',1);bank=d.bank([]);
    if strcmp(task,'inverse'),kinds={'pressure_ar','motion_ridge'};
    else,kinds={'kinematic_ridge'};T=kinematic(T);V=kinematic(V);end
    for kind=kinds
        kind=kind{1};proto=d.bank(1);proto.name=kind;proto.net=[];
        proto.contextMean=mean(T.Context);proto.contextStd=std(T.Context);
        proto.contextStd(proto.contextStd<1e-8)=1;
        Y=T.Y;if strcmp(kind,'pressure_ar'),Y=T.Base+T.Y;end
        proto.targetMean=mean(Y);proto.targetStd=std(Y);proto.targetStd(proto.targetStd<1e-8)=1;
        proto.kind=kind;if strcmp(kind,'kinematic_ridge'),proto.kind='ridge';end
        X=features(T,proto);best=Inf;
        for lambda=d.cfg.ridgeLambda
            m=proto;A=[X ones(size(X,1),1)];P=eye(size(A,2));P(end,end)=0;
            Z=(Y-m.targetMean)./m.targetStd;
            m.W=[A;sqrt(lambda*size(A,1))*P]\[Z;zeros(size(A,2),2)];
            row=metric(task,m,V);if row.macroRMSE<best,best=row.macroRMSE;chosen=m;end
        end
        bank(end+1)=chosen;
    end
    save(fullfile(folder,[task '_frozen.mat']),'bank','protocol'); % Before external loading.
    report=table;
    for k=1:numel(bank),r=metric(task,bank(k),V);r.partition="validation";report=[report;r];end %#ok<AGROW>
    E=residualBenchmarkData(dataset,task,'external_test',1);
    if strcmp(task,'forward'),E=kinematic(E);end
    for k=1:numel(bank),r=metric(task,bank(k),E);r.partition="external_test";report=[report;r];end %#ok<AGROW>
    writetable(report,fullfile(folder,[task '_summary.csv']));disp(report);
end
end
function D=kinematic(D)
D.Base=D.Current+D.Velocity.*D.dt;D.Context=[D.dt D.Base];
D.Y=1000*(D.Truth-D.Base);
end
function X=features(D,m)
c=(D.Context-m.contextMean)./m.contextStd;
X=[double(reshape(D.Seq,[],size(D.Seq,3))') c];
if strcmp(m.kind,'pressure_ar')
    X=[double(reshape(D.Seq(5:7,:,:),[],size(D.Seq,3))') c(:,end-2)];
elseif strcmp(m.kind,'motion_ridge')
    for j=0:size(D.Seq,2)-1,X(:,8*j+(5:7))=0;end
end
end
function row=metric(task,m,D)
R=residualPredict(m,D);
if strcmp(task,'inverse')
    B=[-1/3 -1/3;2/3 -1/3;-1/3 2/3];pred=D.p0+(D.Base+R)*B';scale=1;
    gap=max(max(-pred,pred-3),0);maxGap=max(gap,[],'all');bad=mean(any(gap>.01,2));
else,pred=D.Base+R/1000;scale=1000;maxGap=NaN;bad=NaN;end
err=scale*(pred-D.Truth);scores=[];
for id=unique(D.Run)',scores(end+1)=sqrt(mean(sum(err(D.Run==id,:).^2,2)));end %#ok<AGROW>
row=table(string(m.name),mean(scores),sqrt(mean(sum(err.^2,2))),maxGap,bad, ...
    'VariableNames',{'model','macroRMSE','pooledRMSE','maximumBoundViolation_bar','violationAbovePoint01Fraction'});
end
