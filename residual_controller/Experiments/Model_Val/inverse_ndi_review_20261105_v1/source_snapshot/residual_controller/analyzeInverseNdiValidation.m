function summary=analyzeInverseNdiValidation(root,folder)
%ANALYZEINVERSENDIVALIDATION Offline inverse/NDI signal audit, no hardware I/O.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
base=fullfile(root,'residual_controller','Experiments','Model_Val');
if nargin<2,folder=fullfile(base,'inverse_ndi_review_20261105_v1');end
assert(~isfolder(folder),'Use a fresh output folder.');mkdir(folder);
[~,~,~,po,io]=setupResidualShadow(root);network=residualControllerNetwork(io(1));
[fields,reasons]=ndiDiagnosticsSchema();writetable(fields,fullfile(folder,'ndi_fields.csv'));writetable(reasons,fullfile(folder,'ndi_reject_bits.csv'));
files={'circlesweep_FWD_INV_test.mat','helix_FWD_INV_test.mat'};summary=table;reasonCounts=table;
for f=1:2
    saved=load(fullfile(base,files{f}),'out');out=saved.out;L=out.logsout;
    [t,residualP]=sig(L,'residualPressurePrediction_bar');N=numel(t);
    [~,nominalP]=sig(L,'nominalPressurePrediction_bar');
    B=bus(L,'targetForce_N');I=struct();for name=fieldnames(B)',[tt,x]=matrix(B.(name{1}));assert(isequal(t,tt));I.(name{1})=x;end
    B=bus(L,'scoreValid');D=struct();for name=fieldnames(B)',[tt,x]=matrix(B.(name{1}));assert(isequal(t,tt));D.(name{1})=x;end
    B=bus(L,'historyPre');C=struct();for name=fieldnames(B)',[tt,x]=matrix(B.(name{1}));assert(isequal(t,tt));C.(name{1})=x;end
    S=struct();names={'sensorValid','orientationValid','baseOreintationValid','configurationValid', ...
        'phiOrientation','pollStart_s','hostTime_s','overrun','position','meas_position', ...
        'meas_positionTip_offsetCorrected','moduleTipXYZValid','raw5','raw6','raw7','dd_length', ...
        'P1_NI_DAQ','P2_NI_DAQ','P3_NI_DAQ'};
    for name=names,[tt,x]=sig(L,name{1});S.(name{1})=interp1(tt,x,t+1e-10,'previous','extrap');end
    [tt,x]=sig(L,'thetaObservable','/MATLAB System');S.thetaObservable=interp1(tt,x,t+1e-10,'previous','extrap');
    [tt,x]=sig(L,'ndiDiagnostics');assert(isequal(t,tt));M=reshape(x,N,32,2);
    p=[S.P1_NI_DAQ,S.P2_NI_DAQ,S.P3_NI_DAQ];
    tip=logical(S.sensorValid(:,2));baseValid=logical(S.sensorValid(:,1));ori=logical(S.orientationValid);
    thetaObs=logical(S.thetaObservable);config=logical(S.configurationValid);valid=logical(I.valid);
    finite=valid&all(isfinite([residualP,nominalP]),2);same=finite&max(abs(residualP-nominalP),[],2)<1e-8;
    positionFinite=all(isfinite(S.position),2);rawTipFinite=all(isfinite([S.raw5,S.raw6,S.raw7]),2);
    positionHeld=[false;max(abs(diff(S.position)),[],2)<1e-12];
    correctedHeld=[false;max(abs(diff(S.meas_positionTip_offsetCorrected)),[],2)<1e-12];
    r=struct();r.file=files{f};r.rows=N;r.endTime_s=t(end);r.stopEvent=out.SimulationMetadata.ExecutionInfo.StopEvent;
    r.hasSimulationError=~isempty(out.SimulationMetadata.ExecutionInfo.ErrorDiagnostic);
    r.baseInvalid=nnz(~baseValid);r.tipInvalid=nnz(~tip);r.orientationInvalid=nnz(~ori);
    r.orientationInvalidWithBothSamplesValid=nnz(~ori&baseValid&tip);
    r.thetaNotObservable=nnz(~thetaObs);r.thetaNotObservableDespiteOrientationValid=nnz(~thetaObs&ori);
    phi=S.phiOrientation(~thetaObs&ori);r.validButThetaUnobservablePhiRange_rad=[min(phi),max(phi)];
    r.configurationInvalidAtNDIRate=nnz(~config);r.configurationInvalidDespiteOrientationValid=nnz(~config&ori);
    r.configurationValidDespiteTipInvalid=nnz(config&~tip);
    r.tipInvalidRawDriverStatusMissing=nnz(~tip&M(:,5,2)==2);
    r.tipInvalidRawDriverUsable=nnz(~tip&M(:,14,2)==1);
    r.tipInvalidRawPositionFinite=nnz(~tip&M(:,11,2)==1);
    r.tipOutOfVolume=nnz(M(:,8,2)==1);r.tipPartialOutOfVolume=nnz(M(:,9,2)==1);
    r.tipRepeatedFrame=nnz(M(:,18,2)==0);r.unpairedFrame=nnz(M(:,19,2)==0);
    r.tipInvalidMaskedRawXYZFinite=nnz(~tip&rawTipFinite);
    r.tipInvalidInternalPositionFinite=nnz(~tip&positionFinite);
    r.tipInvalidInternalPositionExactlyHeld=nnz(~tip&positionHeld);
    publicHeld=[false;max(abs(diff(S.meas_position)),[],2)<1e-12];
    r.tipInvalidPublicXYZFinite=nnz(~tip&all(isfinite(S.meas_position),2));
    r.tipInvalidPublicXYZExactlyHeld=nnz(~tip&publicHeld);
    r.orientationInvalidCorrectedTipExactlyHeld=nnz(~ori&correctedHeld);
    r.orientationInvalidCorrectedTipFinite=nnz(~ori&all(isfinite(S.meas_positionTip_offsetCorrected),2));
    r.forwardReady=nnz(D.ready);r.forwardScored=nnz(D.scoreValid);r.inverseValid=nnz(valid);
    r.forwardRejectedDespiteValidOrientation=nnz(ori&~logical(D.accepted));
    r.forwardAcceptedDespiteInvalidOrientation=nnz(~ori&logical(D.accepted));
    r.forwardRejectedThetaUnobservableOnly=nnz(ori&~thetaObs&~logical(D.accepted));
    r.orientationInvalidAndOverrun=nnz(~ori&logical(S.overrun));
    r.inverseFinite=nnz(finite);r.inverseEqual=nnz(same);r.fractionInverseEqual=r.inverseEqual/max(1,r.inverseFinite);
    r.pressureDifferenceMax_bar=max(abs(residualP(finite,:)-nominalP(finite,:)),[],'all');
    r.inverseConverged=nnz(valid&I.converged);r.inverseLimited=nnz(valid&I.limited);
    r.nominalLimited=nnz(valid&I.nominalLimited);r.inverseIterations=unique(I.iterations(valid))';
    r.residualInverseErrorMedianP95Max_N=stats(vecnorm(I.inverseError_N(valid,:),2,2));
    r.nominalInverseErrorMedianP95Max_N=stats(vecnorm(I.nominalInverseError_N(valid,:),2,2));
    r.rawNeuralForceMedianP95Max_N=stats(vecnorm(D.rawForce_N(logical(C.valid),:),2,2));
    r.maxStepFromActual_bar=max(abs(residualP(finite,:)-p(finite,:)),[],'all');
    r.maximumPollInterval_s=max(diff(S.pollStart_s));r.overrunRows=nnz(S.overrun);
    % No assumption that equal curves mean no learned correction. Reproduce
    % logged inverse outputs using source-default options, both common modes.
    probe=find(finite);probe=probe(unique(round(linspace(1,numel(probe),min(120,numel(probe))))));
    differences=zeros(1,2);convergedMismatch=zeros(1,2);
    for common=0:1
        for k=probe'
            ctx=makeContext(C,k);[pr,pn,d]=rlResidualInverseShadowBlock(ctx,S.dd_length(k,:)',po,logical(common),io);
            differences(common+1)=max(differences(common+1),max(abs([pr'-residualP(k,:),pn'-nominalP(k,:)])));
            convergedMismatch(common+1)=convergedMismatch(common+1)+double(d.converged~=logical(I.converged(k)));
        end
    end
    r.defaultReplayMaxDifference_fixedVariable_bar=differences;
    r.defaultReplayConvergenceMismatches_fixedVariable=convergedMismatch;
    % Non-exclusive diagnostic bit counts; pair-invalid appears on BASE too,
    % and does not mean the base tracker itself failed.
    for sensor=1:2
        mask=uint32(M(:,23,sensor));
        for bit=1:height(reasons)
            present=bitand(mask,uint32(reasons.value(bit)))~=0;
            reasonCounts=[reasonCounts;table(string(files{f}),sensor,string(reasons.condition{bit}),nnz(present), ...
                nnz(present&~ori),'VariableNames',{'source','sensor','condition','rows','orientationInvalidRows'})]; %#ok<AGROW>
        end
    end
    % Selected equal-command cases, chosen by largest current residual force,
    % demonstrate the effect of the hypothetical slew cap OFFLINE only.
    equalRows=find(same);[~,order]=sort(vecnorm(D.rawForce_N(equalRows,:),2,2),'descend');
    selected=equalRows(order(1:min(8,numel(order))));counterfactual=table;
    [~,bestMode]=min(differences);common=bestMode==2;if differences(2)<1e-8,common=true;end
    for k=selected'
        ctx=makeContext(C,k);noSlew=po;noSlew(8)=1e6;
        [pr,pn,d]=rlResidualInverseShadowBlock(ctx,S.dd_length(k,:)',noSlew,common,io);
        lower=max(po(4),p(k,:)'-po(8)*ctx.horizon_s);upper=min(po(5),p(k,:)'+po(8)*ctx.horizon_s);
        atEdge=abs(residualP(k,:)'-lower)<1e-8|abs(residualP(k,:)'-upper)<1e-8;
        % Feasible grid checks whether the logged local projected solve can
        % miss a better pressure inside its SAME box/rails/correction caps.
        bestError=norm(I.inverseError_N(k,:));bestP=residualP(k,:)';
        if common
            cfg=residualPressureConfig(po,[0;0;1;.1;.1;1;Inf],true);
            for a1=0:4
                for a2=0:4
                    for a3=0:4
                        candidate=lower+(upper-lower).*[a1;a2;a3]/4;
                        s=(candidate(2)+candidate(3))/2-candidate(1);aa=(candidate(3)-candidate(2))/2;
                        if abs(s-po(2))>po(6)+1e-9||abs(aa-po(3))>po(7)+1e-9,continue;end
                        force=residualForceMap(candidate,ctx.measured,ctx.featureDt_s,ctx.historyPre,network,cfg,zeros(2,1));
                        err=norm(force-I.targetForce_N(k,:)');
                        if err<bestError,bestError=err;bestP=candidate;end
                    end
                end
            end
        end
        counterfactual=[counterfactual;table(t(k),sum(atEdge),residualP(k,:),nominalP(k,:),pr',pn', ...
            max(abs(pr-pn)),norm(I.inverseError_N(k,:)),norm(d.inverseError_N),bestError,bestP', ...
            'VariableNames',{'time_s','axesOnSlewOrRailEdge','loggedResidual_bar','loggedNominal_bar', ...
            'noSlewResidual_bar','noSlewNominal_bar','noSlewDifferenceMax_bar','loggedInverseError_N', ...
            'noSlewInverseError_N','sameBoundsGridBestError_N','sameBoundsGridBestPressure_bar'})]; %#ok<AGROW>
    end
    tag=erase(files{f},'.mat');caseFolder=fullfile(folder,tag);mkdir(caseFolder);
    event=~ori|~config|~thetaObs;
    events=table(t(event),baseValid(event),tip(event),ori(event),config(event),thetaObs(event), ...
        S.phiOrientation(event),M(event,5,2),M(event,8,2),M(event,9,2),M(event,18,2),M(event,19,2), ...
        M(event,23,2),rawTipFinite(event),positionFinite(event),positionHeld(event),logical(D.ready(event)),logical(I.valid(event)), ...
        'VariableNames',{'time_s','baseValid','tipValid','orientationValid','configurationValid','thetaObservable', ...
        'phiOrientation_rad','tipTrackerStatus','tipOutOfVolume','tipPartialOutOfVolume','tipFrameChanged','pairedFrames', ...
        'tipRejectMask','maskedRawTipXYZFinite','internalPositionFinite','internalPositionHeld','forwardReady','inverseValid'});
    writetable(events,fullfile(caseFolder,'validity_events.csv'));writetable(counterfactual,fullfile(caseFolder,'offline_constraint_counterfactuals.csv'));
    f0=fopen(fullfile(caseFolder,'summary.json'),'w');fwrite(f0,jsonencode(r,PrettyPrint=true));fclose(f0);
    save(fullfile(caseFolder,'audit.mat'),'r','t','S','I','D','C','M','p','residualP','nominalP','same','finite','probe','counterfactual');
    fig=figure('Visible','off','Color','w','Position',[20 20 1250 1000]);tiledlayout(4,1);
    nexttile;stairs(t,double(tip));hold on;stairs(t,double(ori),'--');stairs(t,double(config),':');ylabel('Flag');legend('Tip sample valid','Orientation valid','XYZ-fit configuration valid','Location','eastoutside');
    nexttile;stairs(t,double(thetaObs));hold on;stairs(t,double(D.ready),'--');ylabel('Flag');legend('Theta observable','Forward ready','Location','eastoutside');
    nexttile;plot(t,residualP,'LineWidth',1.2);hold on;plot(t,nominalP,'--');ylabel('Predicted pressure (bar)');
    nexttile;plot(t,vecnorm(I.inverseError_N,2,2));hold on;plot(t,vecnorm(I.nominalInverseError_N,2,2),'--');ylabel('Numerical inverse error (N)');xlabel('Simulation time (s)');legend('Residual model','Nominal model','Location','eastoutside');
    ax=findall(fig,'Type','axes');for j=1:numel(ax),grid(ax(j),'on');xlim(ax(j),[t(1) t(end)]);end
    sgtitle(strrep(tag,'_',' '));exportgraphics(fig,fullfile(caseFolder,'inverse_and_validity.png'),'Resolution',180);
    exportgraphics(fig,fullfile(caseFolder,'inverse_and_validity.pdf'),'ContentType','vector');savefig(fig,fullfile(caseFolder,'inverse_and_validity.fig'));close(fig);
    summary=[summary;table(string(files{f}),N,r.tipInvalid,r.orientationInvalid,r.thetaNotObservableDespiteOrientationValid, ...
        r.configurationInvalidAtNDIRate,r.inverseValid,r.inverseConverged,r.fractionInverseEqual,r.pressureDifferenceMax_bar, ...
        'VariableNames',{'source','rows','tipInvalid','orientationInvalid','thetaUnobservableWithValidOrientation', ...
        'configurationInvalid','inverseValid','inverseConverged','fractionInverseEqual','maxPressureDifference_bar'})]; %#ok<AGROW>
    disp(r);disp(counterfactual);
end
writetable(summary,fullfile(folder,'summary.csv'));writetable(reasonCounts,fullfile(folder,'ndi_reason_counts.csv'));disp(summary);
end

function c=makeContext(C,k)
c=residualShadowContext();names=fieldnames(c);
for j=1:numel(names),v=C.(names{j})(k,:)';if strcmp(names{j},'valid'),v=logical(v);end;c.(names{j})=v;end
end
function B=bus(L,field)
B=[];
for j=1:L.numElements,e=L.getElement(j);if isstruct(e.Values)&&isfield(e.Values,field),B=e.Values;return;end;end
assert(~isempty(B),'Missing bus with field %s',field);
end
function v=stats(x)
x=x(isfinite(x));if isempty(x),v=[NaN NaN NaN];else,v=[median(x) prctile(x,95) max(x)];end
end
function [t,x]=sig(L,name,path)
if nargin<3,path='';end
for j=1:L.numElements
    e=L.getElement(j);p=strjoin(string(convertToCell(e.BlockPath)),'|');
    if strcmp(regexprep(e.Name,'[<>]',''),name)&&contains(p,path),[t,x]=matrix(e.Values);return;end
end
error('Missing signal %s (%s).',name,path);
end
function [t,x]=matrix(v)
t=double(v.Time(:));raw=double(v.Data);
if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
