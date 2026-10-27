function manifest=pacTuneTrialManifest()
%PACTUNETRIALMANIFEST Operator-completed metadata for one manually run trial.
manifest.schemaVersion='pac-manifest-1';
manifest.kind='screening'; % screening or confirmation
manifest.conditions=pacTuneManualConfig().conditions;
manifest.resetVerified=false;
manifest.operator='';
manifest.notes='';
end
