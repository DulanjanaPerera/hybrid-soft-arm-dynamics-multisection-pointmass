function benchmarkResidualModels(dataset,outFolder)
%BENCHMARKRESIDUALMODELS Offline first screening; no Simulink or hardware.
assert(~isfolder(outFolder));mkdir(outFolder);
before=rng;cleanup=onCleanup(@()rng(before)); %#ok<NASGU>
cfg=struct('seed',73,'trainingStride',3,'epochs',20,'batchSize',256, ...
    'learningRate',.001,'MLPwidths',[32 16],'GRUwidth',16, ...
    'ridgeLambda',[1e-6 1e-3 .1],'selection','macro whole-run validation error', ...
    'dataset',dataset,'status','running','ridgeSolver','double augmented QR least squares');
save(fullfile(outFolder,'protocol.mat'),'cfg');
for task={'forward','inverse'}
    task=task{1};folder=fullfile(outFolder,task);mkdir(folder);
    fprintf('Preparing %s development partitions only...\n',task);
    T=residualBenchmarkData(dataset,task,'train',cfg.trainingStride);
    V=residualBenchmarkData(dataset,task,'validation',1);
    assert(all(startsWith(T.Run,'day1_'))&&all(startsWith(V.Run,'day1_')));
    proto=struct('name','','kind','','W',[],'net',[], ...
        'contextMean',mean(T.Context),'contextStd',std(T.Context), ...
        'targetMean',mean(T.Y),'targetStd',std(T.Y));
    proto.contextStd(proto.contextStd<1e-8)=1;proto.targetStd(proto.targetStd<1e-8)=1;
    X=[double(reshape(T.Seq,[],size(T.Seq,3))') (T.Context-proto.contextMean)./proto.contextStd];
    XV=[double(reshape(V.Seq,[],size(V.Seq,3))') (V.Context-proto.contextMean)./proto.contextStd];
    Y=single((T.Y-proto.targetMean)./proto.targetStd);
    YV=single((V.Y-proto.targetMean)./proto.targetStd);
    bank=repmat(proto,0,1);validation=table;details=table;
    [s,p]=score(task,'nominal',V,zeros(size(V.Y)));validation=[validation;s];details=[details;p];
    [s,p]=score(task,'persistence',V,baselineResidual(task,V,'persistence'));validation=[validation;s];details=[details;p];
    if strcmp(task,'forward')
        [s,p]=score(task,'constant_velocity',V,baselineResidual(task,V,'velocity'));validation=[validation;s];details=[details;p];
    end
    best=Inf;bestRidge=proto;
    for lambda=cfg.ridgeLambda
        model=proto;model.kind='ridge';model.name=sprintf('ridge_%g',lambda);
        A=[X ones(size(X,1),1)];pen=eye(size(A,2));pen(end,end)=0;
        model.W=[A;sqrt(lambda*size(A,1))*pen]\[double(Y);zeros(size(A,2),2)];
        [s,~]=score(task,model.name,V,residualPredict(model,V));
        if s.macroRMSE<best,best=s.macroRMSE;bestRidge=model;end
    end
    bestRidge.name='ridge';bank(end+1)=bestRidge;
    for kind={'mlp','gru'}
        kind=kind{1};rng(cfg.seed);model=proto;model.kind=kind;model.name=kind;
        if strcmp(kind,'mlp')
            layers=[featureInputLayer(size(X,2),'Normalization','none'); ...
                fullyConnectedLayer(32);tanhLayer;fullyConnectedLayer(16);tanhLayer; ...
                fullyConnectedLayer(2);regressionLayer];
            inputs=single(X);validInputs=single(XV);
        else
            layers=[sequenceInputLayer(8+size(T.Context,2),'Normalization','none'); ...
                gruLayer(16,'OutputMode','last');fullyConnectedLayer(2);regressionLayer];
            inputs=seqCells(T,proto);validInputs=seqCells(V,proto);
        end
        opts=trainingOptions('adam','MaxEpochs',cfg.epochs,'MiniBatchSize',cfg.batchSize, ...
            'InitialLearnRate',cfg.learningRate,'L2Regularization',1e-4, ...
            'Shuffle','every-epoch','ValidationData',{validInputs,YV}, ...
            'ValidationFrequency',ceil(size(Y,1)/cfg.batchSize),'ValidationPatience',4, ...
            'OutputNetwork','best-validation-loss','ExecutionEnvironment','cpu', ...
            'Verbose',false,'Plots','none');
        fprintf('Training %s %s (%d train anchors)...\n',task,kind,size(Y,1));ticId=tic;
        [model.net,info]=trainNetwork(inputs,Y,layers,opts);elapsed=toc(ticId);
        save(fullfile(folder,[kind '_training.mat']),'info','elapsed','cfg');
        bank(end+1)=model;
    end
    for i=1:numel(bank)
        [s,p]=score(task,bank(i).name,V,residualPredict(bank(i),V));
        validation=[validation;s];details=[details;p];
    end
    [~,i]=min(validation.macroRMSE);selected=validation.model(i);
    % Freeze all candidates and selection BEFORE any day2 sequence is loaded.
    save(fullfile(folder,'frozen_selection.mat'),'bank','selected','validation','cfg');
    writetable(validation,fullfile(folder,'validation_summary.csv'));
    writetable(details,fullfile(folder,'validation_by_run.csv'));
    fprintf('%s validation selected %s. Frozen before external evaluation.\n',task,selected);
    E=residualBenchmarkData(dataset,task,'external_test',1);
    external=table;perRun=table;
    names={'nominal','persistence'};if strcmp(task,'forward'),names{end+1}='constant_velocity';end
    for name=names
        if strcmp(name{1},'nominal'),R=zeros(size(E.Y));
        elseif strcmp(name{1},'constant_velocity'),R=baselineResidual(task,E,'velocity');
        else,R=baselineResidual(task,E,'persistence');end
        [s,p]=score(task,name{1},E,R);external=[external;s];perRun=[perRun;p];
    end
    for i=1:numel(bank)
        [s,p]=score(task,bank(i).name,E,residualPredict(bank(i),E));external=[external;s];perRun=[perRun;p];
    end
    writetable(external,fullfile(folder,'external_summary.csv'));
    writetable(perRun,fullfile(folder,'external_by_run.csv'));
    save(fullfile(folder,'evaluation.mat'),'validation','external','selected','perRun');
    if strcmp(task,'forward'),residualBenchmarkRollouts(dataset,folder,bank,proto);end
end
cfg.status='complete';save(fullfile(outFolder,'protocol.mat'),'cfg');
end
function cells=seqCells(D,m)
c=(D.Context-m.contextMean)./m.contextStd;
seq=cat(1,D.Seq,permute(repmat(single(c),1,1,size(D.Seq,2)),[2 3 1]));
cells=cell(size(seq,3),1);for j=1:numel(cells),cells{j}=seq(:,:,j);end
end
function R=baselineResidual(task,D,kind)
if strcmp(task,'inverse'),R=D.Previous-D.Base;
elseif strcmp(kind,'velocity'),R=1000*(D.Current+D.Velocity.*D.dt-D.Base);
else,R=1000*(D.Current-D.Base);end
end
function [summary,detail]=score(task,name,D,R)
if strcmp(task,'inverse')
    B=[-1/3 -1/3;2/3 -1/3;-1/3 2/3];pred=D.p0+(D.Base+R)*B';unit="bar";
    violation=mean(any(pred< -1e-8|pred>3+1e-8,2));scale=1;
else,pred=D.Base+R/1000;unit="mm";violation=NaN;scale=1000;end
err=scale*(pred-D.Truth);detail=table;values=[];
for id=unique(D.Run)'
    m=D.Run==id;rmse=sqrt(mean(sum(err(m,:).^2,2)));values(end+1)=rmse; %#ok<AGROW>
    detail=[detail;table(string(name),id,nnz(m),rmse,'VariableNames',{'model','runId','anchors','vectorRMSE'})]; %#ok<AGROW>
end
summary=table(string(name),unit,size(err,1),mean(values),sqrt(mean(sum(err.^2,2))), ...
    violation,'VariableNames',{'model','units','anchors','macroRMSE','pooledRMSE','pressureViolationFraction'});
end
