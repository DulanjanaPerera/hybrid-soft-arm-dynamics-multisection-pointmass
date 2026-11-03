function auditResidualInterfaces(root,outputFolder)
%AUDITRESIDUALINTERFACES Read both acquisition cohorts; no training/model run.
assert(~isfolder(outputFolder),'Choose a new output directory.');
cohorts={'day1_root','day2_20261027'};
locations={fullfile(root,'RL_data_collection'),fullfile(root,'RL_data_collection','20261027')};
manifest=struct([]);inventory=struct([]);
for day=1:2
    files=dir(fullfile(locations{day},'*.mat'));
    for k=1:numel(files)
        path=fullfile(files(k).folder,files(k).name);d=load(path,'out');
        assert(isfield(d,'out'),'Missing out: %s',path);ds=d.out.logsout;
        runId=[cohorts{day} '__' erase(files(k).name,'.mat')];
        meta=d.out.SimulationMetadata;
        runRow=struct('cohort',cohorts{day},'runId',runId, ...
            'source',path,'bytes',files(k).bytes,'fileModified',files(k).date, ...
            'model',char(meta.ModelInfo.ModelName),'modelVersion',char(meta.ModelInfo.ModelVersion), ...
            'timingMetadata',jsonencode(meta.TimingInfo),'signals',ds.numElements); %#ok<AGROW>
        if isempty(manifest),manifest=runRow;else,manifest(end+1)=runRow;end %#ok<AGROW>
        for j=1:ds.numElements
            e=ds{j};p=e.BlockPath;parts=cell(1,p.getLength);
            for m=1:p.getLength,parts{m}=p.getBlock(m);end
            block=strjoin(parts,' | ');name=char(e.Name);v=e.Values;
            samples=0;channels=0;t0=NaN;t1=NaN;dt=NaN;lo='';hi='';units='';
            if isa(v,'timeseries')&&(isnumeric(v.Data)||islogical(v.Data))
                t=double(v.Time(:));samples=numel(t);
                if samples>0
                    x=double(v.Data);
                    if v.IsTimeFirst,x=reshape(x,samples,[]);else,x=reshape(x,[],samples)';end
                    channels=size(x,2);t0=t(1);t1=t(end);dt=median(diff(t));
                    lo=mat2str(min(x,[],1,'omitnan'),8);hi=mat2str(max(x,[],1,'omitnan'),8);
                    units=char(v.DataInfo.Units);
                end
            end
            role='unclassified';
            if ~isempty(regexp(name,'^P[123#]_NI_DAQ$','once'))&&contains(block,'Pneumatic Valves')
                role='final_command_candidate';
            elseif contains(lower(name),'pressure')&&contains(lower(block),'inverse')
                role='inverse_prediction_not_physical_input';
            elseif contains(block,'MATLAB System')&&ismember(name, ...
                    {'thetaOrientation','phiOrientation','lengthChnage_m','lengthChange_m','meas_length'})
                role='physical_response_orientation_derived';
            elseif contains(lower(name),'des_')||contains(lower(name),'reference')
                role='reference_or_prediction_check_provenance';
            elseif contains(lower(name),'valid')||contains(lower(name),'frame')|| ...
                    contains(lower(name),'host')||contains(lower(name),'poll')
                role='quality_or_time';
            end
            row=struct('cohort',cohorts{day},'runId',runId,'index',j, ...
                'name',name,'blockPath',block,'role',role,'samples',samples, ...
                'channels',channels,'start_s',t0,'end_s',t1,'medianDt_sim_s',dt, ...
                'units',units,'minimum',lo,'maximum',hi); %#ok<AGROW>
            if isempty(inventory),inventory=row;else,inventory(end+1)=row;end %#ok<AGROW>
        end
        fprintf('Audited %s (%d signals)\n',runId,ds.numElements);
    end
end
mkdir(outputFolder);runs=struct2table(manifest);signals=struct2table(inventory);
writetable(runs,fullfile(outputFolder,'runs.csv'));
writetable(signals,fullfile(outputFolder,'signals.csv'));
writetable(signals(strcmp(signals.role,'final_command_candidate'),:), ...
    fullfile(outputFolder,'physical_command_candidates.csv'));
save(fullfile(outputFolder,'interface_audit.mat'),'runs','signals');
assert(height(runs)==10,'Expected five files per cohort; inspect manifest if changed.');
end
