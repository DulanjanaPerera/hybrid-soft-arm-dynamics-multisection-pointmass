function info=prepareResidualControllerNetworks(root)
%PREPARERESIDUALCONTROLLERNETWORKS Offline packing only, no model/hardware run.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_dynamics_learning','force_model');folder=fullfile(root,'residual_controller');
addpath(root,folder,fullfile(root,'passivity_adaptive_control'),fullfile(base,'streaming'));
paths={fullfile(base,'fine_tuning','campaign_20261103_v1','release','force_model.mat'), ...
    fullfile(base,'capacity_comparison','campaign_20261103_v1','candidate_force_model.mat')};
hashes={'1c38edc36ba6ff4a41da8f974a88ca1c5ea775655426da91124a35b255a9c870', ...
    '8f2e0bb55ea9399b62abd5abe0b8f6dd15ab0470a400ad3e593e1edcd8f9f545'};
info=table;
for j=1:2
    assert(strcmp(sgHash(paths{j}),hashes{j}),'Frozen model hash mismatch.');a=load(paths{j});m=a.model;
    assert(numel(m.cfg.widths)==1);H=m.cfg.widths;assert(ismember(H,[64 128]));
    assert(isequal(m.cfg.nominal.K,650*[2 1;1 2])&&isequal(m.cfg.nominal.D,40*eye(2)));
    assert(m.cfg.nominal.boundMode==1&&isequal(m.cfg.nominal.g,[0;0;-9.81]));
    p=m.params;network=struct('input',zeros(384,8),'recurrent',zeros(384,128), ...
        'bias',zeros(384,1),'output',zeros(2,128),'skip',p.skip,'outputBias',p.outputBias, ...
        'mean',m.normalizer.mean,'std',m.normalizer.std,'forceScale',m.cfg.forceScale_N, ...
        'forceBound',m.cfg.forceBound_N,'units',double(H),'forcePerBar',pi*.0065^2*1e5);
    for gate=0:2
        rows=gate*128+(1:H);source=gate*H+(1:H);
        network.input(rows,:)=p.l1Input(source,:);network.recurrent(rows,1:H)=p.l1Recurrent(source,:);
        network.bias(rows)=p.l1Bias(source);
    end
    network.output(:,1:H)=p.output;
    for field=fieldnames(network)',assert(all(isfinite(network.(field{1})),'all'));end
    assert(all(network.std>0));
    path=fullfile(folder,sprintf('residualNetwork%d_v1.mat',H));
    if isfile(path)
        old=load(path);assert(isequaln(old.network,network),'Never overwrite a different packed network.');
    else
        save(path,'network','-v7');
    end
    info=[info;table(H,string(paths{j}),string(hashes{j}),string(path),string(sgHash(path)), ...
        'VariableNames',{'units','source','sourceSHA256','packed','packedSHA256'})]; %#ok<AGROW>
end
writetable(info,fullfile(folder,'NETWORKS.csv'));disp(info(:,[1 3 5]));
end
