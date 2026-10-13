function results=analyzePacStaticSweep(filePath,makePlots)
%ANALYZEPACSTATICS WEEP Detect all static holds from logged references.
% Offline only. Uses logsout, not the assumed reference-block schedule.
% Saves *_sweep_report.txt, *_sweep_analysis.mat and *_sweep_summary.csv.
if nargin<2,makePlots=true;end
S=load(filePath);ds=S.out.logsout;
v=signal(ds,'des_phi');t=v.Time(:);phi=double(v.Data(:));
v=signal(ds,'des_theta');theta=interp1(v.Time,double(v.Data(:)),t,'previous',NaN);
constant=[false;abs(diff(phi))<1e-7 & abs(atan2(sin(diff(theta)),cos(diff(theta))))<1e-7];
mask=constant & phi>1e-4;
beg=find(diff([false;mask])==1);fin=find(diff([mask;false])==-1);
keep=t(fin)-t(beg)>=10;beg=beg(keep);fin=fin(keep);
assert(~isempty(beg),'No >=10 s constant nonzero holds found.');
n=numel(beg);rows=nan(n,21);runs=cell(n,1);
report=sprintf('PAC STATIC SWEEP REPORT: %s\nDetected %d holds from logged reference.\n',filePath,n);
for i=1:n
    window=[t(beg(i)),t(fin(i))];
    try
        captured=evalc('r=analyzePacStatic(filePath,window,false,false);'); %#ok<NASGU>
        runs{i}=r;a=r.signals;
        ix=r.time>=r.settledStart_s & r.time<=window(2);
        good=ix & a.modelValid==1 & a.orientationValid==1 & a.status==3;
        assert(any(good),'No valid fully engaged settled samples.');
        q=mean(a.fitered_length(good,:),1);e=mean(r.lengthError_m(good,:),1);
        f=mean(a.feedbackForce_N(good,:),1);
        rows(i,:)=[i mean(a.des_theta(good)) mean(a.des_phi(good)) ...
            mean(a.phiOrientation(good)) mean(abs(r.phiError_rad(good))) ...
            mean(abs(r.thetaError_rad(good & a.phiOrientation>.05))) ...
            q*1000 e*1000 f mean(a.forceBias_N(good,:),1) ...
            mean(a.limited(ix)) mean(a.modelValid(ix)) r.inferredLambda r.inferredKs ...
            mean(a.allAboveActivation(ix))];
        report=[report sprintf('\nHold %d %.2f..%.2f s: theta/phiDesired %.4f/%.4f; phiActual/MAE %.4f/%.4f; thetaMAE %.4f rad\n  length error %s mm; feedback %s N; bias %s N\n  lambda %s; Ks %s; valid %.3f; limited %.3f; aboveActivation %.3f\n  pressures %s bar; engaged %d/%d; max |allocationResidual| %s N\n',i,window,rows(i,2:6),mat2str(e*1000,5),mat2str(f,5),mat2str(rows(i,13:14),5),mat2str(r.inferredLambda,5),mat2str(r.inferredKs,5),rows(i,16),rows(i,15),rows(i,21),mat2str(mean(r.pressure_bar(good,:),1),5),nnz(good),nnz(ix),mat2str(max(abs(a.allocationResidual_N(good,:)),[],1),5))];
    catch err
        report=[report sprintf('\nHOLD %d %.2f..%.2f FAILED: %s\n',i,window,err.message)];
    end
end
headers={'holdId','theta_rad','phiDesired_rad','phiActual_rad','phiMAE_rad', ...
 'thetaMAE_rad','q2_mm','q3_mm','e2_mm','e3_mm','feedback2_N','feedback3_N', ...
 'bias2_N','bias3_N','limitedFraction','validFraction','lambda2','lambda3', ...
 'Ks2','Ks3','aboveActivationFraction'};
summary=array2table(rows,'VariableNames',headers);
% Exploratory static residual: feedback ~= deltaK*q+b. Scale predictors to mm.
usable=all(isfinite(rows(:,[7 8 11 12])),2) & rows(:,15)==0 & rows(:,16)>.95;
fit=struct('available',false);
X=[rows(usable,7:8),ones(nnz(usable),1)];Y=rows(usable,11:12);
if size(X,1)>=4 && rank(X)==3
    B=X\Y;pred=X*B;fit.available=true;
    fit.deltaK_NPerM=B(1:2,:)'*1000;fit.offset_N=B(3,:)';
    fit.rmse_N=sqrt(mean((pred-Y).^2,1));fit.predictorCondition=cond(X);
    report=[report sprintf('\nEXPLORATORY feedback ~= deltaK*q+b\nDeltaK (N/m): %s\nb (N): %s\nTraining RMSE (N): %s; scaled predictor condition %.4g\n', ...
        mat2str(fit.deltaK_NPerM,5),mat2str(fit.offset_N',5),mat2str(fit.rmse_N,5),fit.predictorCondition)];
    report=[report sprintf('This fit combines mechanics, pressure mapping and trim errors; it is NOT physical stiffness identification.\nNo symmetry/positive-definiteness constraint imposed. Do not insert directly into the controller.\n')];
end
report=[report sprintf('\nUse repeated theta=0/phi=.2 holds to assess drift. Compare their bias, errors and feedback.\nPlateau is defined by the reference; physical settling is not guaranteed. Inspect plots for drift/oscillation.\n')];
fprintf('%s',report);disp(summary);
results=struct('file',filePath,'summary',summary,'holds',{runs},'fit',fit,'report',report);
[folder,base]=fileparts(filePath);
fid=fopen(fullfile(folder,[base '_sweep_report.txt']),'w');assert(fid>=0);
fprintf(fid,'%s\n%s',report,evalc('disp(summary)'));fclose(fid);
save(fullfile(folder,[base '_sweep_analysis.mat']),'results');
writetable(summary,fullfile(folder,[base '_sweep_summary.csv']));
if makePlots
    figure('Name','PAC static sweep summary');tiledlayout(2,2);
    nexttile;plot(rows(:,1),rows(:,3:4),'o-');legend('desired phi','actual phi');ylabel('rad');
    nexttile;plot(rows(:,1),rows(:,9:10),'o-');legend('l2 error','l3 error');ylabel('mm');
    nexttile;plot(rows(:,1),rows(:,11:12),'o-');legend('feedback 2','feedback 3');ylabel('N');
    nexttile;plot(rows(:,1),rows(:,6),'o-');ylabel('theta MAE (rad)');xlabel('hold');
end
end
function v=signal(ds,name)
ix=find(strcmp(ds.getElementNames,name));assert(numel(ix)==1,'Need one signal %s.',name);v=ds{ix}.Values;
end
