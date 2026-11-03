function y=residualPredict(model,D)
%RESIDUALPREDICT Frozen transforms/net, outputs residual in task units.
c=(D.Context-model.contextMean)./model.contextStd;
X=[double(reshape(D.Seq,[],size(D.Seq,3))') c];
if strcmp(model.kind,'pressure_ar')
    X=[double(reshape(D.Seq(5:7,:,:),[],size(D.Seq,3))') c(:,end-2)];
    z=[X ones(size(X,1),1)]*model.W;
elseif strcmp(model.kind,'motion_ridge')
    for j=0:size(D.Seq,2)-1,X(:,8*j+(5:7))=0;end
    z=[X ones(size(X,1),1)]*model.W;
elseif strcmp(model.kind,'ridge'),z=[X ones(size(X,1),1)]*model.W;
elseif strcmp(model.kind,'mlp'),z=predict(model.net,single(X),'MiniBatchSize',512);
else
    seq=cat(1,D.Seq,permute(repmat(single(c),1,1,size(D.Seq,2)),[2 3 1]));
    cells=cell(size(seq,3),1);for j=1:numel(cells),cells{j}=seq(:,:,j);end
    z=predict(model.net,cells,'MiniBatchSize',512);
end
y=double(z).*model.targetStd+model.targetMean;
if strcmp(model.kind,'pressure_ar'),y=y-D.Base;end % Convert predicted d to residual.
end
