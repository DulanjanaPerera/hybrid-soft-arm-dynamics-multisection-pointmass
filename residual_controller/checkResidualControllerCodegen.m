function results=checkResidualControllerCodegen(root,folder)
%CHECKRESIDUALCONTROLLERCODEGEN Standalone C emission only; no SLX/build/run.
assert(~isfolder(folder));mkdir(folder);addpath(root,fullfile(root,'residual_controller'),fullfile(root,'passivity_adaptive_control'));
cfg=coder.config('lib');cfg.GenCodeOnly=true;cfg.GenerateReport=false;
[ga,sc,mr,lo,up,ini,ao,ho,po,so]=pacTunableDefaults();ro=residualControllerOptions();
pdfb={zeros(6,1),zeros(2,1),true,0,[0;0;-9.81],true,0,true,true,false,po,so,ro};
pac={zeros(6,1),zeros(6,1),[0;0;-9.81],true,true,0,0,true,[28;28],[5;5],true,false, ...
    ga,sc,mr,lo,up,ini,ao,ho,po,so,false,true,ro};
results=table;
for j=1:2
    name='pdfbResidualInverseBlock';args=pdfb;if j==2,name='pacResidualAdaptiveTunableBlock';args=pac;end
    target=fullfile(folder,name);ok=false;message="";
    try
        codegen('-config',cfg,name,'-args',args,'-d',target);ok=true;
    catch err
        message=string(getReport(err,'extended','hyperlinks','off'));
    end
    files=dir(fullfile(target,'**','*.c'));
    results=[results;table(string(name),ok,numel(files),message,'VariableNames',{'entrypoint','generationReturnedSuccess','cFilesEmitted','message'})]; %#ok<AGROW>
    writetable(results,fullfile(folder,'status.csv'));fprintf('CODEGEN %s: returned%d, C files%d. %s\n',name,ok,numel(files),message);
end
save(fullfile(folder,'status.mat'),'results');
end
