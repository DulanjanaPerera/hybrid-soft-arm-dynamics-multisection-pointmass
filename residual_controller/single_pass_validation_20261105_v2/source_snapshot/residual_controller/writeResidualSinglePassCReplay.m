function writeResidualSinglePassCReplay(root,folder)
%WRITERESIDUALSINGLEPASSCREPLAY Reference fixtures only, no SLX/device I/O.
assert(~isfolder(folder),'Choose a fresh fixture folder.');mkdir(folder);
[~,~,~,po]=setupResidualShadow(root);[so,~]=setupResidualSinglePassShadow(root);
inputs=[];expected=[];
for units=[64 128]
    so(1)=units;
    for k=0:100
        t=.05*k;c=residualShadowContext();c.units=units;
        c.valid=~ismember(k,[0 1 2 40 41]);c.measured=[.012*sin(t);.01*cos(t);.012*cos(t);-.01*sin(t)];
        c.historyPre(1:units)=.1*sin((1:units)'*.2+t);c.sampleTime_host_s=100+t;
        c.pressure_bar=[1.5+1.4*sin(t);1.5+1.3*cos(t);1.5+.8*sin(2*t)];
        acceleration=.02*[cos(t);sin(t)];
        for common=[false true]
            [p,pn,d]=rlResidualSinglePassShadowBlock(c,acceleration,po,common,so);
            inputs(end+1,:)=[flat(c);acceleration;po;common;so]'; %#ok<AGROW>
            expected(end+1,:)=[p;pn;flat(d)]'; %#ok<AGROW>
        end
    end
end
% Real saved histories and pressures, including dropout/cold contexts.
audit=fullfile(root,'residual_controller','Experiments','Model_Val','inverse_ndi_review_20261105_v1');
for name={'circlesweep_FWD_INV_test','helix_FWD_INV_test'}
    a=load(fullfile(audit,name{1},'audit.mat'));so(1)=128;
    for k=unique(round(linspace(1,numel(a.t),180)))
        c=residualShadowContext();
        for field=fieldnames(c)',v=a.C.(field{1})(k,:)';if strcmp(field{1},'valid'),v=logical(v);end;c.(field{1})=v;end
        acceleration=a.S.dd_length(k,:)';
        [p,pn,d]=rlResidualSinglePassShadowBlock(c,acceleration,po,true,so);
        inputs(end+1,:)=[flat(c);acceleration;po;true;so]'; %#ok<AGROW>
        expected(end+1,:)=[p;pn;flat(d)]'; %#ok<AGROW>
    end
end
writematrix(inputs,fullfile(folder,'inputs.csv'));writematrix(expected,fullfile(folder,'expected.csv'));
fprintf('Wrote %d single-pass native-reference rows.\n',size(inputs,1));
end
function x=flat(s)
x=[];for name=fieldnames(s)',v=double(s.(name{1}));x=[x;v(:)];end %#ok<AGROW>
end
