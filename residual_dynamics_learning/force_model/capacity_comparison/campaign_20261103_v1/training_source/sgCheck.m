function sgCheck(root)
base=fullfile(root,'residual_dynamics_learning','force_model');addpath(base,fullfile(base,'streaming'));
cfg=sgConfig();data=fullfile(base,'streaming','data_20261103_v1');
if ~isfolder(data),sgPrepare(fullfile(base,'dataset_20261103_v1'),data,cfg);end
D=load(fullfile(data,'data.mat'));
jobs=zeros(4,2);train=find(arrayfun(@(s)strcmp(s.partition,'train'),D.segments));
for i=1:4,jobs(i,:)=[train(i) 1];end
B=sgBatch(D,jobs,cfg);
assert(all(sum(B.auxDt.*B.auxMask,3)<=cfg.horizon_s+1e-9));
assert(all(sum(B.dt.*B.mask,3)<=cfg.chunk_s+1e-9));
for widths={32,64,[32 32]}
    cfg.widths=widths{1};params=sgInit(cfg);numeric=sgPlain(params);
    % Nonzero head makes recurrent and flow gradient checks informative.
    numeric.output=.02*randn(size(numeric.output));
    params=numeric;for f=fieldnames(params)',params.(f{1})=dlarray(params.(f{1}));end
    H=sgZero(cfg,4);tic;[loss,g,carry]=dlfeval(@sgLoss,params,B,H,cfg,D.normalizer);elapsed=toc;
    assert(isfinite(extractdata(loss)));
    for f=fieldnames(g)',assert(all(isfinite(extractdata(g.(f{1}))),'all'));end
    plus=numeric;minus=numeric;epsilon=1e-5;
    plus.skip(1,1)=plus.skip(1,1)+epsilon;minus.skip(1,1)=minus.skip(1,1)-epsilon;
    for f=fieldnames(plus)',plus.(f{1})=dlarray(plus.(f{1}));minus.(f{1})=dlarray(minus.(f{1}));end
    a=dlfeval(@sgLoss,plus,B,H,cfg,D.normalizer);b=dlfeval(@sgLoss,minus,B,H,cfg,D.normalizer);
    fd=(extractdata(a)-extractdata(b))/(2*epsilon);ad=extractdata(g.skip(1,1));
    assert(abs(fd-ad)<1e-4*max(1,abs(fd)),'Gradient mismatch.');
    h=H;
    for j=1:size(B.x,3)
        old=h;h=sgCell(numeric,(B.x(:,:,j)-D.normalizer.mean)./D.normalizer.std,h,cfg);
        for l=1:numel(h),h{l}=h{l}.*B.mask(:,:,j)+old{l}.*(1-B.mask(:,:,j));end
    end
    for l=1:numel(h),assert(max(abs(h{l}-extractdata(carry{l})),[],'all')<1e-12);end
    B2=B;B2.auxTruth=B.auxTruth+1;B2.truth=B.truth-1;
    [~,~,other]=dlfeval(@sgLoss,params,B2,H,cfg,D.normalizer);
    for l=1:numel(h),assert(isequal(extractdata(other{l}),extractdata(carry{l})));end
    fprintf('PASS widths %s: gradient AD %.8g / FD %.8g, measured-state purity, physical horizon; 4-batch gradient %.3fs.\n',mat2str(cfg.widths),ad,fd,elapsed);
end
fprintf('ALL STREAMING CHECKS PASS.\n');
end
