function [s,command]=pacTuneOnlineStep(s,u,c)
%PACTUNEONLINESTEP Event-driven HOST supervisor. Pure state in/state out.
% Commands are requests, never hardware actions. See ONLINE_WORKFLOW.md.
% u: t,healthy,holdReady,appliedId,appliedGamma,scoreReady,scoreId,
% metrics=[l2mm l3mm velocityRMS roughness],coverage,limitedFraction.
% scoreReady means the entire batch, including faults/limits, was audited.
assert(isfinite(u.t)&&isscalar(u.t),'pacTune:OnlineTime','Invalid clock.');
if isempty(s)
    assert(c.period_s>0&&c.circles>=c.scoreCircles&&c.scoreCircles>=1);
    assert(c.maxGammaFactor>=1&&c.confirmations>=1&&c.settle_s>=0);
    assert(all(c.targets>0)&&all(isfinite(c.targets)));
    assert(strcmp(c.statePolicy,'preserve'),'Only preserve-state protocol is supported.');
    % Validate optimizer configuration before producing any command.
    pacTuneOptimizerSuggest(c.optimizer,pacTuneOptimizerObservations());
    s=struct('phase','hold','since',u.t,'lastTime',u.t,'id',0, ...
        'gamma',c.optimizer.baselineGamma(:),'observations',pacTuneOptimizerObservations(), ...
        'consecutivePasses',0,'reason','','events',{{}});
end
command=struct('mode','hold','apply',false,'id',s.id,'gamma',s.gamma, ...
    'preserveAdaptiveState',true,'phase',s.phase,'reason',s.reason);
if strcmp(s.phase,'stopped'),command.mode='stop';return;end
if ~u.healthy||u.t<s.lastTime
    s=stop(s,'unhealthy_or_clock_reversal');
end
s.lastTime=u.t;
switch s.phase
    case 'hold'
        if u.holdReady
            if height(s.observations)>=c.optimizer.maxObservations
                s=stop(s,'budget_exhausted');
            else
                s.id=s.id+1;s.phase='apply';s.since=u.t;
            end
        elseif u.t-s.since>c.holdTimeout_s
            s=stop(s,'hold_timeout');
        end
    case 'apply'
        if ~u.holdReady
            s=stop(s,'lost_hold_during_apply');
        elseif u.appliedId==s.id && isequal(u.appliedGamma(:),s.gamma)
            s.phase='settle';s.since=u.t;
        elseif u.t-s.since>c.applyTimeout_s
            s=stop(s,'apply_timeout');
        end
    case 'settle'
        if ~u.holdReady,s=stop(s,'lost_hold_during_settle');
        elseif u.t-s.since>=c.settle_s
            s.phase='run';s.since=u.t;
        end
    case 'run'
        if u.t-s.since>=c.period_s*c.circles
            s.phase='score';s.since=u.t;
        end
    case 'score'
        if u.scoreReady&&u.scoreId==s.id&&u.holdReady
            m=u.metrics(:)';
            if numel(m)~=4||any(~isfinite(m))||any(m<0)|| ...
                    ~isfinite(u.coverage)||u.coverage<c.minimumCoverage||u.coverage>1|| ...
                    ~isfinite(u.limitedFraction)||u.limitedFraction<0||u.limitedFraction>1
                s=stop(s,'invalid_score');
            elseif u.limitedFraction>0
                s=stop(s,'pressure_limiting');
            else
                pass=all(m<=c.targets);
                s.observations(end+1,:)={string(s.id),s.gamma(1),s.gamma(2), ...
                    s.gamma(3),m(4)/c.targets(4),true,true,pass,c.period_s*c.circles};
                s.events{end+1}=struct('id',s.id,'gamma',s.gamma,'metrics',m, ...
                    'time',u.t,'statePolicy',c.statePolicy);
                if pass,s.consecutivePasses=s.consecutivePasses+1;
                else,s.consecutivePasses=0;end
                if s.consecutivePasses>=c.confirmations
                    s=stop(s,'targets_repeated_operationally');
                elseif height(s.observations)>=c.optimizer.maxObservations
                    s=stop(s,'budget_exhausted');
                elseif pass
                    % Repeat this exact candidate before exploring further.
                    s.phase='hold';s.since=u.t;
                else
                    % Local bounded exploration; preserve actual tested gains.
                    oc=c.optimizer;
                    oc.gammaBounds(:,1)=max(oc.gammaBounds(:,1),s.gamma/c.maxGammaFactor);
                    oc.gammaBounds(:,2)=min(oc.gammaBounds(:,2),s.gamma*c.maxGammaFactor);
                    oc.baselineGamma=s.gamma;
                    % Existing BO normalizes history against its bounds: train
                    % only on observations inside this local trust region.
                    o=s.observations;g=o{:,{'gammaK','gammaD','gammaH'}};
                    inside=all(g>=oc.gammaBounds(:,1)'&g<=oc.gammaBounds(:,2)',2);
                    ticId=tic;p=pacTuneOptimizerSuggest(oc,o(inside,:));elapsed=toc(ticId);
                    if elapsed>c.optimizerTimeout_s
                        s=stop(s,'optimizer_timeout');
                    elseif ~strcmp(p.action,'propose')
                        s=stop(s,p.reason);
                    else
                        s.gamma=p.gamma;s.phase='hold';s.since=u.t;
                    end
                end
            end
        elseif u.t-s.since>c.holdTimeout_s
            s=stop(s,'score_or_hold_timeout');
        end
end
command.id=s.id;command.gamma=s.gamma;command.phase=s.phase;command.reason=s.reason;
command.apply=strcmp(s.phase,'apply');
if strcmp(s.phase,'run'),command.mode='run';
elseif strcmp(s.phase,'stopped'),command.mode='stop';
else,command.mode='hold';end
end

function s=stop(s,reason)
s.phase='stopped';s.reason=reason;
end
