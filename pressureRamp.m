function tau = pressureRamp(tspan, maxP, start, endt)

tau = zeros(size(tspan));
minP = 0.0;
minP_idx = tspan<=start;
tau(minP_idx) = minP;
m = (maxP-minP) / (endt - start);
maxP_idx = star<tspan && tspan <= endt;
tau(maxP_idx) = m.*1

    
end