function [sensors,conditions]=decodeNdiDiagnostics(diagnostics)
%DECODENDIDIAGNOSTICS Read one32x2 diagnostic sample without device access.
assert(isequal(size(diagnostics),[32 2])&&all(diagnostics(1,:)==1), ...
    'Expected one schema1 diagnostic sample,32x2(base,tip).');
[fields,reasons]=ndiDiagnosticsSchema();
sensors=array2table(diagnostics','VariableNames',fields.field');
sensors.sensor=["base";"tip"];sensors=movevars(sensors,'sensor','Before',1);
conditions=table(strings(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{'sensor','value','condition'});
for j=1:2
    mask=diagnostics(23,j);
    if ~isfinite(mask)
        conditions=[conditions;{sensors.sensor(j),NaN,"reader_stage_unavailable"}]; %#ok<AGROW>
        continue;
    end
    assert(mask>=0&&mask==floor(mask)&&mask<=double(intmax('uint32')));
    for k=1:height(reasons)
        if bitand(uint32(mask),uint32(reasons.value(k)))~=0
            conditions=[conditions;{sensors.sensor(j),reasons.value(k),string(reasons.condition{k})}]; %#ok<AGROW>
        end
    end
end
end
