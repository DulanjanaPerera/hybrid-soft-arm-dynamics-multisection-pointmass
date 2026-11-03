function sgAdditionalChecks(root)
base=fullfile(root,'residual_dynamics_learning','force_model');src=fullfile(base,'streaming');addpath(base,src);
D=load(fullfile(src,'data_20261103_v1','data.mat'));cfg=sgConfig();cfg.widths=[32 32];p=sgInit(cfg,false);
p.output=.03*randn(size(p.output));jobs=[1 1;2 1;3 1];B=sgBatch(D,jobs,cfg);H=sgZero(cfg,3);
[a,r]=sgFree(p,B,H,cfg,D.normalizer);copy=B;copy.truth=randn(size(B.truth));copy.auxTruth=randn(size(B.auxTruth));
[b,s]=sgFree(p,copy,H,cfg,D.normalizer);assert(isequal(a,b)&&isequal(r,s),'Future target leaked into prediction.');
% Splitting a measured prefix at an arbitrary chunk boundary preserves both layers.
x=(B.x-D.normalizer.mean)./D.normalizer.std;state=H;
for k=1:size(x,3),state=sgCell(p,x(:,:,k),state,cfg);end
split=floor(size(x,3)/2);other=H;
for k=1:split,other=sgCell(p,x(:,:,k),other,cfg);end
other=cellfun(@double,other,'UniformOutput',false);
for k=split+1:size(x,3),other=sgCell(p,x(:,:,k),other,cfg);end
assert(isequal(state,other));
for hz=[.5 1]
    cfg.horizon_s=hz;
    for sid=1:numel(D.segments)
        seg=D.segments(sid);S=D.sequences{seg.run};
        for k=seg.starts
            last=find(S.physicalTime<=S.physicalTime(k)+hz+1e-10,1,'last');
            rows=k:last;assert(all(S.mature(rows))&&all(S.segment(rows)==seg.sourceSegment));
            assert(all(S.split(rows)==string(seg.partition))&&S.physicalTime(last)-S.physicalTime(k)<=hz+1e-9);
        end
    end
end
fid=fopen(fullfile(src,'additional_checks_20261103.txt'),'w');
fprintf(fid,'PASS future-target perturbation leaves free predictions/forces identical.\nPASS stacked state replay invariant to chunk splitting.\nPASS every common .5/1s train/selection anchor stays within physical horizon, mature valid segment and original whole-trial partition.\n');fclose(fid);
fprintf('ADDITIONAL CHECKS PASS.\n');
end
