function Data = f20220211_1_batchGenDynaResponses_varParams_1(T, A)
% Dr. Isuru's original dynamic model. THis function sweeps two variables and generate data from the simulink


% clc
% close all

% tau=[.25 .5 1 2 3 4];
tau=T;
% amp=[1 2 3 4 5 6];
amp=A;
trial=10;

n=2000;
sRate=10;

Data = zeros(n*sRate,11,trial);

% for i=1:length(amp)
%     for j=1:length(tau)
        tic
        for k=1:trial
            fprintf("Aplitude=%0.4f,  Tau=%0.4f,  Trial=%d \n",amp, tau, k)
            %             f20181022_3_genDynaResponses_5(n,tau(j),amp(i),sRate,k)
            Data(:,:,k)=f20220211_2_genDynaResponses_varParams_1(n,tau,amp,sRate,k);
        end
        toc
%     end
% end

end

