function [I] = calcI_alt(COV)

% ------------------------------------------------------------------------------
% COV:      covariance matrix of X
% N:        size of X
% I:        integration
% ------------------------------------------------------------------------------
% Returns the integration for a system X.
% System dynamics is characterised by its covariance matrix COV.
% Utilizes....
% Olaf Sporns, Indiana University, 2003
% ------------------------------------------------------------------------------

N = size(COV,1);
pie1 = 2*pi*exp(1);

U = chol(COV);
d = 2*sum(log(diag(U)));

I = (sum(log(pie1*diag(COV)))/2) - d/2 - log(pie1^N)/2;
