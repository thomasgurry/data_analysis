function [ CI_vals ] = calc_weight_CI( alphas )

% Calculates the 95% confidence intervals for each weight

covariances = zeros(length(alphas));
alpha0 = sum(alphas);
denom_factor = (alpha0^2)*(alpha0 + 1);

for i = 1:length(alphas)
    for j = 1:length(alphas)
        
        if(i==j)
            kronDel = 1;
        else
            kronDel = 0;
        end
        
        covariances(i,j) = (alphas(i)*alpha0*kronDel - alphas(i)*alphas(j))/denom_factor;
        
    end    
end

CI_vals = zeros(1,length(alphas));

for i = 1:length(alphas)
    
    CI_vals(i) = 1.54*1.96*sqrt(covariances(i,i));

end

end
