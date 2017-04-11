% Function for calculating the fraction of beta-content from a set of
% phi/psi angles and a value for the tau parameter.
% For further details about the tau parameter, see 
% Vitalis A, et al.  Biophys. J. 97(1):303-311 (2009).

function [f_beta] = calc_fbeta(phi_vals,psi_vals,tau,r_beta)

    f_beta_indiv = zeros(length(phi_vals),1); % f_beta for each residue
    phi_b = -152;
    psi_b = 142;
    r_b = r_beta;
    
    phi_vals = radtodeg(phi_vals);
    psi_vals = radtodeg(psi_vals); 
    
    for i = 1:length(phi_vals)
        
        % Take minimum angle difference        
        phi_diff = phi_vals(i)-phi_b;
        if(phi_diff < 0)
            shifted_phi_diff = phi_diff + 360;
            if(abs(shifted_phi_diff) < abs(phi_diff))
                phi_diff = shifted_phi_diff;
            end
        else
            shifted_phi_diff = phi_diff - 360;
            if(abs(shifted_phi_diff) < abs(phi_diff))
                phi_diff = shifted_phi_diff;
            end
        end
        
        % Take minimum angle difference
        psi_diff = psi_vals(i)-psi_b;
        if(psi_diff < 0)
            shifted_psi_diff = psi_diff + 360;
            if(abs(shifted_psi_diff) < abs(psi_diff))
                psi_diff = shifted_psi_diff;
            end
        else
            shifted_psi_diff = psi_diff - 360;
            if(abs(shifted_psi_diff) < abs(psi_diff))
                psi_diff = shifted_psi_diff;
            end
        end
        
        d=(phi_diff)^2 + (psi_diff)^2;
        d=sqrt(d);
        d = d-r_b;
        
        if(d <= 0)
            f_beta_indiv(i) = 1;
        else
            d = d^2;
            f_beta_indiv(i) = exp(-tau*d);
        end
        
    end
    
    f_beta = sum(f_beta_indiv)/length(phi_vals);
    
end
        
        
        