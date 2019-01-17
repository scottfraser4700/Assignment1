function [Phi,dPhidr] = LJ(r, Epsilon, Sigma)
%LJ returns leonard jones potential for a given atomic configuration

Phi = 4*Epsilon*(Sigma^12*r^-12 - Sigma^6*r^-6);
dPhidr = 4*Epsilon*(Sigma^12*(-12)*r^-13);
end

