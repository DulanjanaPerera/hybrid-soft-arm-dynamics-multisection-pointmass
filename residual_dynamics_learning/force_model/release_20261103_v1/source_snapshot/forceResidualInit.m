function params=forceResidualInit(cfg,deep)
%FORCERESIDUALINIT New force-network weights; never reuse position corrections.
if nargin<2,deep=true;end
rng(cfg.seed);H=cfg.hiddenUnits;
params=struct('input',randn(3*H,8)*sqrt(2/(8+H)), ...
    'recurrent',[orth(randn(H));orth(randn(H));orth(randn(H))], ...
    'bias',[zeros(H,1);ones(H,1);zeros(H,1)], ...
    'output',zeros(2,H),'skip',zeros(2,8),'outputBias',zeros(2,1));
if deep,for name=fieldnames(params)',params.(name{1})=dlarray(params.(name{1}));end,end
end
