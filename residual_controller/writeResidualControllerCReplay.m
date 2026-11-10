function writeResidualControllerCReplay(root,folder)
%WRITERESIDUALCONTROLLERCREPLAY Generate pure component fixtures for native C.
assert(~isfolder(folder));mkdir(folder);addpath(root,fullfile(root,'residual_controller'),fullfile(root,'passivity_adaptive_control'));
[ga,sc,mr,lo,up,ini,ao,ho,~,~]=pacTunableDefaults();manifest=table;
for kind=1:2
    for units=[64 128]
        for common=[false true]
            for enabled=[false true]
                clear pdfbResidualInverseBlock pacResidualAdaptiveTunableBlock
                [po,so,ro]=pdfbResidualSettings();so(1:3)=[1;3;1];ro(1)=units;
                if kind==2,po(9)=.8;so(6)=.0005;end
                X=[];Y=[];
                for k=0:160
                    t=k*.05;q=zeros(2,1);v=q;
                    if t>4,q=.004*[sin(t);cos(t)];v=.004*[cos(t);-sin(t)];end
                    m=[q;v;zeros(2,1)];des=[q+.001;v;zeros(2,1)];valid=~ismember(k,[95 131 132 133]);
                    % Small base-gravity changes exercise exact re-expression.
                    gravity=[.01*sin(t);-.02*sin(t);-9.81];
                    if kind==1
                        args={m,[.01;-.01],valid,t,gravity,true,t,true,enabled,common,po,so,ro};out=cell(1,10);
                        [out{:}]=pdfbResidualInverseBlock(args{:});
                    else
                        args={m,des,gravity,valid,true,t,t,true,[28;28],[5;5],true,common,ga,sc,mr,lo,up,ini,ao,ho,po,so,false,enabled,ro};out=cell(1,29);
                        [out{:}]=pacResidualAdaptiveTunableBlock(args{:});
                    end
                    X=[X;flat(args)];d=out{end};out=out(1:end-1);
                    for field=fieldnames(d)',out{end+1}=d.(field{1});end %#ok<AGROW>
                    Y=[Y;flat(out)]; %#ok<AGROW>
                end
                id=sprintf('kind%d_units%d_common%d_enabled%d',kind,units,common,enabled);
                writematrix(X,fullfile(folder,[id '_input.csv']));writematrix(Y,fullfile(folder,[id '_expected.csv']));
                manifest=[manifest;table(string(id),kind,units,common,enabled,size(X,1),size(X,2),size(Y,2), ...
                    'VariableNames',{'id','controller','units','variableCommon','enabled','ticks','inputs','outputs'})]; %#ok<AGROW>
            end
        end
    end
end
writetable(manifest,fullfile(folder,'cases.csv'));fprintf('Wrote%d frozen component replay cases. No hardware/model run.\n',height(manifest));
end
function x=flat(values)
x=[];for j=1:numel(values),x=[x reshape(double(values{j}),1,[])];end %#ok<AGROW>
end
