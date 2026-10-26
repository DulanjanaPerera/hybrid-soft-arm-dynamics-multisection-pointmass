root=fileparts(fileparts(mfilename('fullpath')));folder=fullfile(root,'passivity_adaptive_control/experimentations/adaptCommonP/1rad_10s_gainTesting_v8_analysis');
S=load(fullfile(fileparts(folder),'1rad_10s_gainTesting_v8.mat'));ds=S.out.logsout;
[t,th]=read(ds,'parameterEstimate');D=struct;
for name={'parameterLower','parameterUpper','parameterMaxRate','adaptationActive','modelValid','pressureCommand_bar','pressureOptions','parameterAtBound'}
 [tt,x]=read(ds,name{1});D.(name{1})=interp1(tt,x,t,'previous','extrap');
end
fid=fopen(fullfile(folder,'bound_and_rate_audit.txt'),'w');fclose(fid);diary(fullfile(folder,'bound_and_rate_audit.txt'));
labels={'k1','k2','k3','d2','d3','h1','h2','h3','b2','b3'};
circle=t>=20&t<320;late=t>=290&t<320;
fprintf('Stored PREUPDATE estimates: circleN%d lateN%d\n',nnz(circle),nnz(late));
for j=1:10
 lo=th(:,j)<=D.parameterLower(:,j)+1e-9;hi=th(:,j)>=D.parameterUpper(:,j)-1e-9;
 fprintf('%s lower%g upper%g circleMin%g circleMax%g lowerCount%d upperCount%d lateLowerCount%d lateUpperCount%d returnLowerCount%d\n',labels{j},D.parameterLower(1,j),D.parameterUpper(1,j),min(th(circle,j)),max(th(circle,j)),nnz(circle&lo),nnz(circle&hi),nnz(late&lo),nnz(late&hi),nnz(t>=320&lo));
end
rates=diff(th)./diff(t);active=D.adaptationActive(1:end-1)>0;
for win=[20 221.65 271.45 290;221.65 271.45 320 320]
 ix=t(1:end-1)>=win(1)&t(1:end-1)<win(2)&active;
 fprintf('RATEWIN%.2f-%.2f Nactive%d capCounts%s peak%s\n',win(1),win(2),nnz(ix),mat2str(sum(abs(rates(ix,:))>=D.parameterMaxRate(ix,:)-1e-6)),mat2str(max(abs(rates(ix,:))),6));
end
for mask={circle,late}
 ix=mask{1};p=D.pressureCommand_bar;steps=diff(p)./diff(t);ig=ix(1:end-1)&ix(2:end);
 fprintf('PRESSUREWIN%g-%g ranges%s %s peakRate%s fraction5xCompressedExceeds2barPerSec%s\n',min(t(ix)),max(t(ix)),mat2str(min(p(ix,:)),6),mat2str(max(p(ix,:)),6),mat2str(max(abs(steps(ig,:))),6),mat2str(100*mean(5*abs(steps(ig,:))>2+1e-8),6));
end
fprintf('2s reference:40 samples/circle at .05s; q amplitude13mm -> velocityPeak%g mm/s accelerationPeak%g mm/s^2\n',13*pi,13*pi^2);
diary off
function [t,x]=read(ds,name)
ix=find(strcmp(ds.getElementNames,name),1);v=ds{ix}.Values;t=double(v.Time(:));raw=double(v.Data);if v.IsTimeFirst,x=reshape(raw,numel(t),[]);else,x=reshape(raw,[],numel(t))';end
end
