function state=sgZero(cfg,N)
state=cell(numel(cfg.widths),1);
for j=1:numel(state),state{j}=zeros(cfg.widths(j),N);end
end
