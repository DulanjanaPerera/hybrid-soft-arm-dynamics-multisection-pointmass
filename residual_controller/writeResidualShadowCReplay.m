function writeResidualShadowCReplay(root,folder)
%WRITERESIDUALSHADOWCREPLAY MATLAB-reference fixtures; no device/model I/O.
assert(~isfolder(folder),'Choose a fresh fixture folder.');mkdir(folder);
[o,~,~,po,io]=setupResidualShadow(root);
cases=table;
for units=[64 128]
    clear rlResidualShadowBlock
    o(1)=units;io(1)=units;inputs=[];outputs=[];inverseInputs=[];inverseOutputs=[];
    offset=0;
    for k=0:220
        t=.05*k;if k==150,offset=offset+.2;end
        host=2000+.0508*k+offset;frame=k;
        if k==200,frame=k-1;end
        q=.003*[sin(.4*t);cos(.4*t)];lengths=[-sum(q);q];
        p=[1.4+.1*sin(t);1.6;1.8+.1*sin(t)];
        valid=~ismember(k,[90 91 130 131 132]);stamp=t;if k==180,stamp=t-.2;end
        reset=k==205;
        [a,b,c,d,context]=rlResidualShadowBlock(lengths,valid,true,host,frame,p,stamp,t,reset,o);
        inputs(end+1,:)=[lengths;valid;true;host;frame;p;stamp;t;reset;o]'; %#ok<AGROW>
        outputs(end+1,:)=[a(:);b(:);c(:);flatten(d);flatten(context)]'; %#ok<AGROW>
        acceleration=.001*[sin(t);cos(t)];
        for common=[false true]
            [pi,pn,di]=rlResidualInverseShadowBlock(context,acceleration,po,common,io);
            inverseInputs(end+1,:)=[flatten(context);acceleration;po;common;io]'; %#ok<AGROW>
            inverseOutputs(end+1,:)=[pi;pn;flatten(di)]'; %#ok<AGROW>
        end
    end
    id=sprintf('shadow%d',units);writematrix(inputs,fullfile(folder,[id '_inputs.csv']));writematrix(outputs,fullfile(folder,[id '_expected.csv']));
    iid=sprintf('inverse%d',units);writematrix(inverseInputs,fullfile(folder,[iid '_inputs.csv']));writematrix(inverseOutputs,fullfile(folder,[iid '_expected.csv']));
    cases=[cases;table(string(id),string('forward'),size(inputs,1),'VariableNames',{'id','kind','rows'}); ...
        table(string(iid),string('inverse'),size(inverseInputs,1),'VariableNames',{'id','kind','rows'})]; %#ok<AGROW>
end
writetable(cases,fullfile(folder,'cases.csv'));disp(cases);
end
function out=flatten(s)
out=[];names=fieldnames(s);for j=1:numel(names),v=double(s.(names{j}));out=[out;v(:)];end %#ok<AGROW>
end
