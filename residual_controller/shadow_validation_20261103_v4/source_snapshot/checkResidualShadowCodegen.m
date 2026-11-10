function results=checkResidualShadowCodegen(root,folder)
%CHECKRESIDUALSHADOWCODEGEN Source emission only; never loads a Simulink model.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
assert(~isfolder(folder),'Choose a fresh output folder.');mkdir(folder);
setupResidualShadow(root);cfg=coder.config('lib');cfg.GenCodeOnly=true;
cfg.GenerateReport=false;cfg.EnableOpenMP=false;results=struct();
cases={'rlResidualShadowBlock','rlResidualInverseShadowBlock'};
arguments={{zeros(3,1),false,false,0,0,zeros(3,1),0,0,false,zeros(11,1)}, ...
    {residualShadowContext(),zeros(2,1),zeros(10,1),false,zeros(8,1)}};
for k=1:numel(cases)
    item=struct('entry',cases{k},'success',false,'message','');
    try
        codegen('-config',cfg,cases{k},'-args',arguments{k},'-d',fullfile(folder,cases{k}));
        item.success=true;item.message='MATLAB Coder C source generation passed; no target or SLX run.';
    catch err
        item.message=getReport(err,'extended','hyperlinks','off');
    end
    results.(cases{k})=item;disp(item);
end
save(fullfile(folder,'codegen_results.mat'),'results');
f=fopen(fullfile(folder,'codegen_results.json'),'w');fwrite(f,jsonencode(results,PrettyPrint=true));fclose(f);
end
