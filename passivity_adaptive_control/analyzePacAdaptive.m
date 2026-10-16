function report=analyzePacAdaptive(filePath)
%ANALYZEPACADAPTIVE Print adaptive diagnostics from an experiment recording.
% analyzePacAdaptive() uses base-workspace out; explicit out also accepted.
% MAT-path input saves its report; workspace input writes no files.
% Offline only. Also use PAC_circle/analyzePacStaticSweep for tracking metrics.
if nargin<1 || isempty(filePath)
    assert(evalin('base','exist(''out'',''var'')')==1, ...
        'Define out in the base workspace or pass out/a MAT file path.');
    filePath=evalin('base','out');
end
fromFile=ischar(filePath) || (isstring(filePath) && isscalar(filePath));
if fromFile
    S=load(filePath); assert(isfield(S,'out'),'Recording must contain out.');
    sourceOut=S.out;
else
    sourceOut=filePath;
    filePath='workspace out';
end
ds=sourceOut.logsout;v=read('parameterEstimate');t=v.Time(:);
th=matrix(v);active=matrix(read('adaptationActive'));bound=matrix(read('parameterAtBound'));
status=matrix(read('status'));limited=matrix(read('limited'));
assert(size(th,2)==10 && numel(status)==numel(t),'Unexpected log dimensions.');
report=sprintf('ADAPTIVE PAC REPORT: %s\nSamples %d; time %.2f..%.2f s\n',filePath,numel(t),t(1),t(end));
report=[report sprintf('Parameter order: k1 k2 k3 d2 d3 h1 h2 h3 b2 b3\nInitial %s\nFinal %s\nMin %s\nMax %s\n', ...
 mat2str(th(1,:),5),mat2str(th(end,:),5),mat2str(min(th),5),mat2str(max(th),5))];
report=[report sprintf('Adaptive updates %d; bound flag %d; limited %d; fault %d samples\n', ...
 nnz(active),nnz(bound),nnz(limited),nnz(status==5))];
for name={'hysteresisForce_N','adaptiveBias_N','adaptiveForce_N','allocationResidual_N'}
 x=matrix(read(name{1}));ix=status==3;
 report=[report sprintf('%s control mean %s; peak abs %s N\n',name{1}, ...
 mat2str(mean(x(ix,:),1,'omitnan'),5),mat2str(max(abs(x(ix,:)),[],1,'omitnan'),5))];
end
fprintf('%s',report);
if fromFile
[folder,base]=fileparts(filePath);fid=fopen(fullfile(folder,[base '_adaptive_report.txt']),'w');
assert(fid>=0);fprintf(fid,'%s',report);fclose(fid);
end
figure('Name','Adaptive parameter history');tiledlayout(4,1);
groups={1:3,4:5,6:8,9:10};labels={'spring coefficients (N/m)','damping (N s/m)','history weights (N)','bias (N)'};
for i=1:4,nexttile;plot(t,th(:,groups{i}));ylabel(labels{i});end
xlabel('time (s)');
 function v=read(name)
  ix=find(strcmp(ds.getElementNames,name));assert(numel(ix)==1,'Need one logged %s.',name);v=ds{ix}.Values;
 end
end
function x=matrix(v)
raw=double(v.Data);n=numel(v.Time);
if v.IsTimeFirst,x=reshape(raw,n,[]);else,x=reshape(raw,[],n)';end
end
