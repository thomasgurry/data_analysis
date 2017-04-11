function [rate_const] = calc_rate(ddG,T)

% Calculates rate constant for a given activation energy and 
% temperature using the Eyring equation

kb = 1.38065e-23;
R = 8.31446;
h = 6.62607e-27;

rate_const = (kb*T/h)*exp(-ddG/(R*T));

end