root=fileparts(fileparts(mfilename('fullpath')));
S=load(fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_tr1.mat')); ds=S.out.logsout;
[t,th]=read(ds,'parameterEstimate');
names={'modelValid','status','limited','adaptationActive','dampingEstimate_NsPerM','commonPressure_bar','pressureCommand_bar','allocationResidual_N','lambda_gains','damping_gains','theta','phi','des_theta','des_phi','orientationValid','thetaObservable','sampled des_length','meas_length'};
D=struct;
for k=1:numel(names)
 [tt,x]=read(ds,names{k}); D.(matlab.lang.makeValidName(names{k}))=interp1(tt,x,t,'previous','extrap');
end
z=all(D.dampingEstimate_NsPerM==0,2);
fprintf('ZERO %d INVALID %d HOLD %d FAULT %d\n',nnz(z),nnz(z & ~D.modelValid),nnz(z & D.status==4),nnz(D.status==5));
fprintf('RATE %s\n',mat2str(max(abs(diff(th)./diff(t))),6));
for name={'lambda_gains','damping_gains'}
 x=D.(name{1}); ix=[1;find(any(diff(x)~=0,2))+1]; fprintf('%s CHANGES [time values]\n',name{1}); disp([t(ix) x(ix,:)]);
end
e=atan2(sin(D.des_theta-D.theta),cos(D.des_theta-D.theta)); ep=D.des_phi-D.phi;
el=D.sampledDes_length-D.meas_length(:,2:3);
fprintf('LAP thetaMAE thetaMean phiMAE lengthRMSEmm validPct adaptPct limitPct commonMinMax residualPeak\n');
for lap=1:30
 ix=t>=20+(lap-1)*10 & t<20+lap*10; good=ix & D.modelValid & D.orientationValid & D.thetaObservable;
 fprintf('%d %.5f %.5f %.5f %s %.2f %.2f %.2f %s %.5f\n',lap,mean(abs(e(good))),mean(e(good)),mean(abs(ep(good))),mat2str(1000*sqrt(mean(el(good,:).^2)),5),100*mean(D.modelValid(ix)),100*mean(D.adaptationActive(ix)),100*mean(D.limited(ix)),mat2str([min(D.commonPressure_bar(ix)) max(D.commonPressure_bar(ix))],5),max(vecnorm(D.allocationResidual_N(ix & D.modelValid,:),2,2)));
end
ix=t>=20 & t<320; fprintf('PRESSURE %s %s\n',mat2str(min(D.pressureCommand_bar(ix,:)),5),mat2str(max(D.pressureCommand_bar(ix,:)),5));
function [t,x]=read(ds,name)
 names=ds.getElementNames; v=ds{find(strcmp(names,name),1)}.Values;t=v.Time(:);raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
