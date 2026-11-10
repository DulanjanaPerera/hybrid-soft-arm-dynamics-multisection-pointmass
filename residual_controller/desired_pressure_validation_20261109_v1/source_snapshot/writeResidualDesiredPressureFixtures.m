function writeResidualDesiredPressureFixtures(root,folder)
%WRITERESIDUALDESIREDPRESSUREFIXTURES Offline native reference data only.
assert(~isfolder(folder),'Use a fresh folder.');mkdir(folder);
[~,~,~,po]=setupResidualShadow(root);[so,~]=setupResidualSinglePassShadow(root);
inputs=[];expected=[];
for units=[64 128]
    so(1)=units;
    for k=0:100
        t=.05*k;c=residualShadowContext();c.units=units;
        c.valid=~ismember(k,[0 1 2 40 41]);c.measured=[.002*sin(t);.001*cos(t);.002*cos(t);-.001*sin(t)];
        c.historyPre(1:units)=.1*sin((1:units)'*.2+t);c.sampleTime_host_s=100+t;
        c.pressure_bar=[1.5+1.4*sin(t);1.5+1.3*cos(t);1.5+.8*sin(2*t)];
        motion=[.012*sin(t);.01*cos(t);.012*cos(t);-.01*sin(t);-.012*sin(t);-.01*cos(t)];
        for common=[false true]
            [p,pn,d]=rlResidualDesiredPressureBlock(c,motion,po,common,so);
            inputs(end+1,:)=[flat(c);motion;po;common;so]'; %#ok<AGROW>
            expected(end+1,:)=[p;pn;flat(d)]'; %#ok<AGROW>
        end
    end
end
a=load(fullfile(root,'residual_controller','Experiments','Model_Val', ...
    'target_mismatch_review_20261109_v1','analysis.mat'),'t','C','qd','vd','ad');so(1)=128;
for k=unique(round(linspace(1,numel(a.t),360)))
    c=residualShadowContext();for f=fieldnames(c)',v=a.C.(f{1})(k,:)';if strcmp(f{1},'valid'),v=logical(v);end;c.(f{1})=v;end
    motion=[a.qd(k,:)';a.vd(k,:)';a.ad(k,:)'];
    [p,pn,d]=rlResidualDesiredPressureBlock(c,motion,po,true,so);
    inputs(end+1,:)=[flat(c);motion;po;true;so]'; %#ok<AGROW>
    expected(end+1,:)=[p;pn;flat(d)]'; %#ok<AGROW>
end
writematrix(inputs,fullfile(folder,'inputs.csv'));writematrix(expected,fullfile(folder,'expected.csv'));
fprintf('Wrote %d desired-pressure native-reference rows.\n',size(inputs,1));
end
function x=flat(s)
x=[];for f=fieldnames(s)',v=double(s.(f{1}));x=[x;v(:)];end %#ok<AGROW>
end
