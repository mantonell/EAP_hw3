%% Numerical solution to Problem 1
beta = 0.89;
mu_C = 0.0189;
mu_D = 0.0189;
gamma = 2;
sigma_C = 0.015;
sigma_D = 0.112;
rho = 0.2;

%% Question 1:
% The term that appears in the infinite sum.
A = beta*exp(mu_D-gamma*mu_C)*...
    exp(.5*(sigma_D^2+gamma^2*sigma_C^2)-gamma*rho*sigma_C*sigma_D);

% The Price-Dividend ratio.
H = A/(1-A);

%% Question 2:
% Expected value of the risky return.
m = (H+1)/H*exp(mu_D+.5*sigma_D^2);
% Standard error of the risky return.
sigma = (H+1)/H*exp(mu_D+.5*sigma_D^2)*sqrt(exp(sigma_D^2)-1);