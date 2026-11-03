function params=sgInit(cfg,deep)
if nargin<2,deep=true;end
rng(cfg.seed);params=struct;input=8;
for layer=1:numel(cfg.widths)
    H=cfg.widths(layer);tag=sprintf('l%d',layer);
    params.([tag 'Input'])=randn(3*H,input)*sqrt(2/(input+H));
    params.([tag 'Recurrent'])=[orth(randn(H));orth(randn(H));orth(randn(H))];
    params.([tag 'Bias'])=[zeros(H,1);ones(H,1);zeros(H,1)];input=H;
end
params.output=zeros(2,input);params.skip=zeros(2,8);params.outputBias=zeros(2,1);
if deep,for f=fieldnames(params)',params.(f{1})=dlarray(params.(f{1}));end,end
end
