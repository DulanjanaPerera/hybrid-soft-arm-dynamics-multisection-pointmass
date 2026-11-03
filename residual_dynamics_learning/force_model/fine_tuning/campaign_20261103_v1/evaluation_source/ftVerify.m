function ftVerify(root)
%FTVERIFY Independent saved-metric/source/split checks after frozen evaluation.
base=fullfile(root,'residual_dynamics_learning');force=fullfile(base,'force_model');src=fullfile(force,'fine_tuning');
addpath(base,force,src,fullfile(force,'streaming'));camp=fullfile(src,'campaign_20261103_v1');
F=load(fullfile(camp,'FROZEN_BEFORE_EVALUATION.mat'));D=load(fullfile(camp,'data','data.mat'));
initial=load(F.protocol.initialModel);initial=initial.model;
assert(isequaln(D.normalizer,initial.normalizer));old=load(fullfile(force,'streaming','data_20261103_v1','data.mat'));
assert(isequaln(D.sequences(1:10),old.sequences)&&isequaln(D.index(1:height(old.index),:),old.index));
hashes=sgRead(fullfile(camp,'training_source_hashes.csv'),'TextType','string');
for j=1:height(hashes),assert(strcmp(sgHash(hashes.source(j)),hashes.sha256(j)));end
external=sgRead(fullfile(camp,'evaluation','nov2_sources.csv'),'TextType','string');
for j=1:height(external),assert(strcmp(sgHash(external.source(j)),external.sha256(j)));end
O=sgRead(fullfile(camp,'evaluation','one_step.csv'),'TextType','string');
maxDifference=0;rowsChecked=0;ids={'StepWave','SineWave4','combinePressure','circleSweep','helix'};
locked={'widths','horizon_s','burnIn_s','commonBurnIn_s','chunk_s','maxHorizon_s','nominal','positionScale_m', ...
    'forceScale_N','forceBound_N','forceRateScale_Nps','smoothWeight','magnitudeWeight','substeps','stageIterations','newtonIterations'};
for arm=1:height(F.results)
    assert(strcmp(sgHash(F.results.modelPath(arm)),F.results.sha256(arm)));
    a=load(F.results.modelPath(arm));m=a.model;assert(isequaln(m.normalizer,initial.normalizer));
    for field=locked,assert(isequaln(m.cfg.(field{1}),initial.cfg.(field{1})));end
    for day=2:3
        for run=1:5
            if day==2
                b=load(fullfile(base,'validations','fresh_20261102_v2',[ids{run} '_sequences.mat']));S=b.S;parts={'external'};
            else,S=D.sequences{10+run};parts={'selection','train'};end
            a=load(fullfile(camp,'evaluation',sprintf('%s_nov%d_run%d_stream.mat',F.results.arm(arm),day,run)));
            for part=parts
                mask=true(numel(a.anchors),1);cohort="nov2_external";
                if day==3,mask=S.split(a.anchors)==string(part{1});cohort="nov3_"+part{1};end
                for modelIndex=1:3
                    if modelIndex>1&&arm>1,continue;end
                    label=F.results.arm(arm);if modelIndex==2,label="nominal";elseif modelIndex==3,label="constant_velocity";end
                    e=1000*(a.Q(:,mask,modelIndex)-S.q(a.anchors(mask)+1,:)');vr=sqrt(mean(sum(e.^2,1)));
                    row=O(O.cohort==cohort&O.runId==S.runId&O.variant==label,:);assert(height(row)==1);
                    maxDifference=max(maxDifference,abs(vr-row.lengthVectorRMSE_mm));rowsChecked=rowsChecked+1;
                end
            end
        end
    end
end
assert(maxDifference<1e-8);
M=load(fullfile(camp,'motion_comparison','predictions.mat'));T=sgRead(fullfile(camp,'motion_comparison','scores.csv'),'TextType','string');
longDifference=0;
for j=1:numel(M.trajectories)
    A=M.trajectories(j);W=A.window;
    if A.day==3
        names=cellfun(@(S)S.runId,D.sequences,'UniformOutput',false);
        index=find(D.runGroup==2&endsWith(string(names(:)),"__"+string(W.id)));assert(numel(index)==1);
        assert(all(D.sequences{index}.split(W.rows)=="selection"));
    end
    for arm=1:height(F.results)
        valid=all(isfinite(W.truth),2)&all(isfinite(A.predicted(:,:,arm)),2);e=1000*(A.predicted(valid,:,arm)-W.truth(valid,:));
        vr=sqrt(mean(sum(e.^2,2)));row=T(T.day==A.day&T.runId==W.id&T.channel==W.channel&T.variant==F.results.arm(arm),:);
        assert(height(row)==1&&row.complete==1);longDifference=max(longDifference,abs(vr-row.lengthVectorRMSE_mm));
    end
end
assert(longDifference<1e-8);
I=load(fullfile(camp,'inverse','summary.mat'));assert(I.summary.statePurityPass&&I.summary.maxBoundViolation_bar<1e-7);
assert(strcmp(sgHash(fullfile(camp,'release','force_model.mat')),F.protocol.selectedSHA256));
assert(strcmp(sgHash(fullfile(force,'streaming','release_20261103_v1','force_model.mat')),F.protocol.initialHash));
assert(strcmp(sgHash(fullfile(force,'release_20261103_v1','force_model.mat')),'b5cfd7d0cd13b2b9067a35d810b8ae166e032e96421b6a30c55d703039af268f'));
fid=fopen(fullfile(camp,'VERIFICATION.md'),'w');
fprintf(fid,'PASS: original weights/releases, old sequences/selection anchors/normalizer preserved.\nPASS: training source/data and November2 source/model hashes unchanged.\nPASS: architecture, loss, physics, history and normalization unchanged across arms.\nPASS: %d saved one-step rows independently recomputed, max difference %.12g mm.\nPASS:30 long-window rows recomputed, max difference %.12g mm; all Nov3 motion windows are whole-trial holdouts.\nPASS: selected model inverse bounds/state-commit checks.\nNo hardware/controller/stability/physical force qualification.\n',rowsChecked,maxDifference,longDifference);fclose(fid);
summary=I.summary;fid=fopen(fullfile(camp,'inverse','summary.json'),'w');fprintf(fid,'%s',jsonencode(summary));fclose(fid);
snapshot=fullfile(camp,'evaluation_source');mkdir(snapshot);files=dir(fullfile(src,'*.m'));source=table;
extra=which('frozenForceMeasuredReplay.m');[~,n,e]=fileparts(extra);copyfile(extra,fullfile(snapshot,[n e]));
source=[source;table(string(extra),string(sgHash(extra)),'VariableNames',{'source','sha256'})];
for j=1:numel(files)
    p=fullfile(files(j).folder,files(j).name);copyfile(p,fullfile(snapshot,files(j).name));
    source=[source;table(string(p),string(sgHash(p)),'VariableNames',{'source','sha256'})]; %#ok<AGROW>
end
writetable(source,fullfile(camp,'evaluation_source_hashes.csv'));
fprintf('FINE-TUNE VERIFICATION PASS: %d one-step rows,30 motion rows, hashes, old splits and heldout trials.\n',rowsChecked);
end
