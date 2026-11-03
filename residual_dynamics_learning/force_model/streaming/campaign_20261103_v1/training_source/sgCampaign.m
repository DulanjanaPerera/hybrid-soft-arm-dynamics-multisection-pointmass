function sgCampaign(root)
base=fullfile(root,'residual_dynamics_learning','force_model');src=fullfile(base,'streaming');addpath(base,src);
folder=fullfile(src,'campaign_20261103_v1');assert(~isfolder(folder));mkdir(folder);
cfg=sgConfig();data=fullfile(src,'data_20261103_v1');
protocol=struct('created',char(datetime('now')),'data',data,'oldOnly',true,'externalLoaded',false, ...
    'seeds',[113 211 307],'ranking','.7 one-step + .3 half-second macro vector RMSE on old holdouts', ...
    'historyRule','Test 1s burn-in on screen leader; replace .5s only if old score improves.', ...
    'finalistRule','Screen top two; leading history variant may replace leader. Three seeds each. No lucky-seed selection.', ...
    'representativeSeed',113,'externalPolicy','Freeze all finalist seeds before November2 evaluation.', ...
    'budget','1200 maximum updates, common patience4 after minimum600; checkpoint every150');
save(fullfile(folder,'PREDECLARED_PROTOCOL.mat'),'protocol','cfg');
mkdir(fullfile(folder,'training_source'));sources=[dir(fullfile(src,'*.m'));dir(fullfile(base,'*.m'))];hashes=table;
for j=1:numel(sources)
    p=fullfile(sources(j).folder,sources(j).name);copyfile(p,fullfile(folder,'training_source',sources(j).name));
    hashes=[hashes;table(string(p),string(sgHash(p)),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
end
p=fullfile(data,'data.mat');hashes=[hashes;table(string(p),string(sgHash(p)),'VariableNames',{'source','sha256'})];
writetable(hashes,fullfile(folder,'source_hashes.csv'));
widths={32,64,[32 32],32,64,[32 32]};horizons=[.5 .5 .5 1 1 1];screen=table;
for j=1:6
    c=cfg;c.widths=widths{j};c.horizon_s=horizons(j);id=char('A'+j-1);
    run=fullfile(folder,[id '_seed113']);fprintf('\nSTART SCREEN %s widths%s horizon%.1fs\n',id,mat2str(c.widths),c.horizon_s);
    model=sgTrain(data,run,c);p=model.params;names=fieldnames(p);count=0;
    for k=1:numel(names),count=count+numel(p.(names{k}));end
    screen=[screen;table(string(id),string(mat2str(c.widths)),c.horizon_s,c.burnIn_s,count, ...
        model.oldSelectionScore_mm,model.selectedUpdate,string(run), ...
        'VariableNames',{'id','widths','horizon_s','burnIn_s','parameters','oldScore_mm','selectedUpdate','folder'})]; %#ok<AGROW>
    writetable(screen,fullfile(folder,'screen.csv'));
end
[~,rank]=sort(screen.oldScore_mm);leaders=screen(rank(1:2),:);
M=load(fullfile(leaders.folder(1),'force_model.mat'));c=M.model.cfg;c.burnIn_s=1;
run=fullfile(folder,[char(leaders.id(1)) '_burn1_seed113']);model=sgTrain(data,run,c);
history=table(leaders.oldScore_mm(1),model.oldSelectionScore_mm,'VariableNames',{'burnHalfScore_mm','burnOneScore_mm'});
writetable(history,fullfile(folder,'history_comparison.csv'));
if model.oldSelectionScore_mm<leaders.oldScore_mm(1)
    leaders.folder(1)=string(run);leaders.burnIn_s(1)=1;leaders.oldScore_mm(1)=model.oldSelectionScore_mm;
    leaders.selectedUpdate(1)=model.selectedUpdate;
end
finalists=table;
for j=1:2
    M=load(fullfile(leaders.folder(j),'force_model.mat'));c=M.model.cfg;
    for seed=protocol.seeds
        c.seed=seed;
        if seed==113,run=char(leaders.folder(j));model=M.model;
        else,run=fullfile(folder,sprintf('%s_final_seed%d',leaders.id(j),seed));model=sgTrain(data,run,c);end
        finalists=[finalists;table(leaders.id(j),seed,string(run),model.oldSelectionScore_mm, ...
            string(sgHash(fullfile(run,'force_model.mat'))), ...
            'VariableNames',{'id','seed','folder','oldScore_mm','sha256'})]; %#ok<AGROW>
        writetable(finalists,fullfile(folder,'finalists.csv'));
    end
end
summary=groupsummary(finalists,'id',{'mean','std'},'oldScore_mm');writetable(summary,fullfile(folder,'finalist_development_summary.csv'));
for j=1:height(hashes),assert(strcmp(sgHash(hashes.source(j)),hashes.sha256(j)),'Training source changed during campaign.');end
D=load(p,'manifest');for j=1:height(D.manifest),assert(strcmp(sgHash(D.manifest.source(j)),D.manifest.sha256(j)));end
protocol.status='ALL_FINALISTS_FROZEN_BEFORE_EXTERNAL';protocol.frozenAt=char(datetime('now'));
save(fullfile(folder,'FINALISTS_FROZEN.mat'),'protocol','leaders','finalists','summary');
fid=fopen(fullfile(folder,'FINALISTS_FROZEN.md'),'w');fprintf(fid,'All six finalist weights frozen using October19/27 only.\nThree seeds per configuration; seed113 is the predeclared representative.\nNo November2 loaded by this campaign. See finalists.csv for hashes.\n');fclose(fid);
disp(screen);disp(summary);fprintf('CAMPAIGN COMPLETE; external evaluation can now begin.\n');
end
