function result=checkResidualSinglePassShadowCodegen(root,folder,entry)
%CHECKRESIDUALSINGLEPASSSHADOWCODEGEN Fixed-size C source emission, no SLX run.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<3,entry='rlResidualSinglePassShadowBlock';end
assert(ismember(entry,{'rlResidualSinglePassShadowBlock','rlResidualSinglePassCollectionShadowBlock'}));
assert(~isfolder(folder),'Choose a fresh output folder.');mkdir(folder);
setupResidualSinglePassShadow(root);cfg=coder.config('lib');cfg.GenCodeOnly=true;
cfg.GenerateReport=false;cfg.EnableOpenMP=false;
result=struct('entry',entry,'sourceGenerationPassed',false,'cSourceEmitted',false,'message','');
try
    codegen('-config',cfg,entry,'-args', ...
        {residualShadowContext(),zeros(2,1),zeros(10,1),false,zeros(2,1)}, ...
        '-d',fullfile(folder,'generated'));
    result.sourceGenerationPassed=true;result.message='C source emission passed; no target/model run.';
catch err
    result.message=getReport(err,'extended','hyperlinks','off');
end
result.cSourceEmitted=isfile(fullfile(folder,'generated',[entry '.c']));
save(fullfile(folder,'codegen_result.mat'),'result');
f=fopen(fullfile(folder,'codegen_result.json'),'w');fwrite(f,jsonencode(result,PrettyPrint=true));fclose(f);
disp(result);
end
