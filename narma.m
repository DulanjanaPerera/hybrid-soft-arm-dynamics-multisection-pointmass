function NMSE = narma(In)
S = 90;
Trial = 10;
Initial = 400;
Training = 1000;
Evaluation = 600;
T = Initial + Training + Evaluation;

u = zeros(T + 2,1);
Input = zeros(T + 2,1);
x = zeros(S + 1,T + 2);
states = zeros(Training,S + 1);
states_tr = zeros(S + 1,Training);

Output_network_10order = 0;
Output_10order = zeros(T + 2,1);
desired_output_10order = zeros(Training,1);
w_out_10order = zeros(S + 1,1);
w_out_10order_final = zeros(S + 1,1);
NMSE_10order = zeros(Trial,1);

try

    Average_NMSE_10order = 0;
    for trial = 1:Trial
        u(1:T,1) = rand(T,1);
        x(S,1:T) = 1.0;
        try

%             in = strcat("Input_csv\trial_",num2str(trial),".csv");
%             input = importfile(in, [1, Inf]);
            input = In(:,:,trial);
            [r,~] = size(input);
            it = floor(r/10);
            for i=1:it
                u(i,1) = input((i-1)*10+1,2);
                temp = reshape(input((i-1)*10+1:i*10,3:11)',[],1);
                x(1:S,i)=temp;
            end
        catch ME
            disp("File reading error");
            rethrow(ME);
        end
        desired_output_10order(:,1) = 0;
        w_out_10order(:,1) = 0;
        w_out_10order_final(:,1) = 0;
        tr = 1;
        noise_level = 0.00001;
        for t = Initial+1:Training + Initial
            for m = 1:S + 1
                if ( rand()> 0.5)
                    states(tr,m) = x(m,t) + noise_level * rand();
                else
                    states(tr,m) = x(m,t) - noise_level * rand();
                end
            end
            tr = tr + 1;
        end
        
        state2 = states'*states;
        mat = (state2+state2')*0.5;
        eigenvalues = real(eig(mat));

        Output_10order(1:T,1) = 0.2;
        Input(2:T,1) = u(1:T-1,1) * 0.2;
        tr = 1;
        for t = 1:T
            if (t > 10 && t < T)
                Output_10order(t + 1,1) = 0.3 * Output_10order(t,1) + 0.05 * Output_10order(t,1) * (Output_10order(t,1) + Output_10order(t - 1,1) + Output_10order(t - 2,1) + Output_10order(t - 3,1) + Output_10order(t - 4,1) + Output_10order(t - 5,1) + Output_10order(t - 6,1) + Output_10order(t - 7,1) + Output_10order(t - 8,1) + Output_10order(t - 9,1))...
                    + 1.5 * Input(t - 9,1) * Input(t,1) + 0.1;
            end
            if (t >= Initial+1 && t < Training + Initial+1)
                desired_output_10order(tr,1) = Output_10order(t + 1,1);
                tr = tr + 1;
            end
        end
        R_10order = desired_output_10order;
        P_10order = states'*R_10order;

        min_AIC_10order = 0;
        for ll = 1:S + 1
            df = ll;
            lambda = S + 1;
%             lambda_next = S + 1;
            for iteration = 0:100
                A = 0;
                B = 0;
                for n = 1:S + 1
                    A = A + eigenvalues(n) / (eigenvalues(n) + lambda);
                    B = B + eigenvalues(n) / ((eigenvalues(n) + lambda) * (eigenvalues(n) + lambda));
                end
                lambda_next = lambda + (A - df) / B;
                lambda = lambda_next;
            end
            if (lambda < 0)
                lambda = 0;
            end
            dia_vector = ones(1,S + 1)*lambda;
            I = diag(dia_vector);
            W_10order = (state2 + I)\P_10order;
            w_out_10order(1:S+1,1) = W_10order(1:S+1, 1); %probably (:,1) is fine.
            RSS_10order = 0;
            for t = Initial+1:Initial + Training
                Output_network_10order = 0;
                for j = 1:S + 1
                    Output_network_10order = Output_network_10order + w_out_10order(j,1) * x(j,t);
                end
%                 Output_network_10order = sum(w_out_10order(1:S + 1,1).*x(1:S + 1,t));
                RSS_10order = RSS_10order + (Output_10order(t + 1,1) - Output_network_10order) * (Output_10order(t + 1,1) - Output_network_10order);
            end

            AIC_10order = Training * log(RSS_10order) + 2.0 * df;
            if ll == 1
                min_AIC_10order = AIC_10order;
            end

            if (min_AIC_10order > AIC_10order)
                min_AIC_10order = AIC_10order;
                w_out_10order_final(:,1) = w_out_10order(:,1);
            end
        end

        NMSE_10order(trial,1) = 0;
        mean_output_10order = 0;
        MSE_10order = 0;

        for t = T - Evaluation+1:T - 1
            Output_network_10order = 0;
            for j = 1:S + 1
                Output_network_10order = Output_network_10order + w_out_10order_final(j,1) * x(j,t);
            end
%             Output_network_10order = (w_out_10order_final(1:S + 1,1)') * (x(1:S + 1,t));
            mean_output_10order = mean_output_10order + (Output_10order(t + 1,1)) *  (Output_10order(t + 1,1));
            MSE_10order = MSE_10order +  (Output_10order(t + 1,1) - Output_network_10order) * (Output_10order(t + 1,1) - Output_network_10order);
        end%T

        NMSE_10order(trial,1) = MSE_10order / mean_output_10order;
        Average_NMSE_10order = Average_NMSE_10order + NMSE_10order(trial,1) / Trial;
    end
    NMSE = Average_NMSE_10order;
%     disp(NMSE);
catch ME
    disp("Program doesn't start");
    rethrow(ME);
end
end