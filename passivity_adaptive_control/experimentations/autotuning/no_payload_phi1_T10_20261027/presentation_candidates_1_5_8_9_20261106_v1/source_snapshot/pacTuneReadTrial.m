function trial=pacTuneReadTrial(file,cfg)
%PACTUNEREADTRIAL Read native logged signals with explicit provenance.
% Does not run models. Ambiguous names require cfg.selectors.<field>.
if nargin<2,cfg=pacTuneCampaignConfig();end
S=load(file,'out');assert(isfield(S,'out'),'pacTune:MissingOut','MAT must contain out.');
ds=S.out.logsout;n=ds.numElements;
index=(1:n)';names=strings(n,1);paths=names;units=names;shapes=names;
first=NaN(n,1);last=first;count=first;
for k=1:n
    e=ds{k};names(k)=string(e.Name);
    paths(k)=strjoin(string(convertToCell(e.BlockPath)),' -> ');
    if isa(e.Values,'timeseries')
        v=e.Values;units(k)=string(v.DataInfo.Units);
        shapes(k)=string(mat2str(size(v.Data)));count(k)=numel(v.Time);
        if count(k)>0,first(k)=v.Time(1);last(k)=v.Time(end);end
    end
end
trial.inventory=table(index,names,paths,units,shapes,count,first,last);
trial.file=char(java.io.File(file).getCanonicalPath());
trial.sha256=hashFile(file);trial.schemaVersion=cfg.schemaVersion;
required={'parameterEstimate','status','modelValid','orientationValid', ...
 'configurationValid','baseOreintationValid','thetaObservable','measurementTime_s', ...
 'supervisorOptions','sampled des_length','meas_length','d_length','filtered_d_length', ...
 'theta','phi','des_theta','des_phi','limited','adaptationActive', ...
 'pressureCommand_bar','allocationResidual_N','adaptationGamma','parameterLower', ...
 'parameterUpper','parameterInitial','parameterMaxRate','parameterScale', ...
 'adaptationOptions','hystoryOptions','pressureOptions','lambda_gains','damping_gains', ...
 'enableAdaptation','enableCommonPressureAdaptation','freezeHysteresisWeights'};
optional={'s_mPerS','historyState','tauBlended_N','forceBias_N', ...
 'dampingForce_N','hysteresisForce_N','feedbackForce_N','inertialForce_N','staticForce_N'};
trial.requiredFields=cellfun(@matlab.lang.makeValidName,required,'UniformOutput',false);
trial.signals=struct;trial.missing=strings(0,1);
for name=[required optional]
    field=matlab.lang.makeValidName(name{1});ix=find(names==string(name{1}));
    if isfield(cfg.selectors,field)
        choice=cfg.selectors.(field);
        if isnumeric(choice),ix=ix(ix==choice);else,ix=ix(paths(ix)==string(choice));end
    end
    if isempty(ix)
        if ismember(name,required),trial.missing(end+1,1)=string(name{1});end
        continue
    end
    assert(isscalar(ix),'pacTune:AmbiguousSignal', ...
        'Signal %s has multiple matches; inspect inventory and set a selector.',name{1});
    v=ds{ix}.Values;assert(isa(v,'timeseries'),'pacTune:SignalType','Expected timeseries.');
    t=double(v.Time(:));raw=double(v.Data);
    if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
    assert(all(isfinite(t))&&all(diff(t)>0),'pacTune:Time','Nonmonotonic time: %s',name{1});
    trial.signals.(field)=struct('t',t,'x',x,'index',ix,'name',name{1}, ...
        'blockPath',char(paths(ix)),'units',char(units(ix)),'rawShape',size(raw));
end
end

function h=hashFile(file)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(file,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid),md.update(fread(fid,1048576,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
