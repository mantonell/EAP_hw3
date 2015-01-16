% Illustration how to discretize an AR(1)
% x_t = b + a x_{t-1} + eps_t
% eps_t N(0,sigeps^2)
clc; clear all;

%% Calibration of the AR(1) model.

stdy = 0.02;
xm = 0.016;
a = 0.5;

b = xm*(1-a);
sigeps = stdy*sqrt(1-a^2);


%% Determine the boundaries of the domain.

alpha   = 0.00001;
mu      = 0;
sigma   = 1;
k_alpha = norminv(1-alpha/2,mu,sigma);
yL      = k_alpha*stdy;

% Number of discretization points.
N = 60;
% Grid of boundaries.
Bk = linspace(-yL,yL,N+1);
Bk = xm+Bk;
% This will be grid with boundaries but where extremes are +/- infty.
Bkub = Bk;
Bkub(1) = -Inf;
Bkub(end)=Inf;
xj    = (Bk(2:N+1)+Bk(1:N))/2;

%% Matrix of transition probabilities and steady state distribution.

P=zeros(N,N); % matrix of transition probabilities
for i=1:N
    P(i,:)=normcdf( (Bkub(2:N+1) - b - a*xj(i))./sigeps,0,1)- ...
           normcdf( (Bkub(1:N)   - b - a*xj(i))./sigeps,0,1);
end

% now hunt steady state distribution
AM=P'-eye(N);
AM(end,:)=ones(1,N);
bv=zeros(N,1);
bv(end)=1;
s= AM\bv; % same as inv(AM)*bv


% and verify again
% figure()
% hold on;
% plot(xj,s,'ok');
% I=normpdf(xj,b/(1-a),stdy);
% plot(xj,I/sum(I),'+m');
% hold off;

%% Value function iteration.

gamma = 2;
beta = 0.89;

T = @(v) repmat(beta*exp(xj*(1-gamma)),[N 1]).*P * (1+v);

% Initialization of the value function algorithm.
H = zeros(N,1);
Hnext = T(H);
% Tolerance.
tol = 0.0001;

while norm(Hnext-H,inf)>tol
    H = Hnext;
    Hnext = T(H);
end

%% Solution as a linear system.

A = eye(N) - repmat(beta*exp(xj*(1-gamma)),[N 1]).*P;
B = repmat(beta*exp(xj*(1-gamma)),[N 1]).*P*ones(N,1);
f = linsolve(A,B);

% Plot.
figure(1);
hold on;
plot(xj,H,'r'); plot(xj,f,'g');
legend('value function iteration','linear system');
title('Price-dividend ratio');
xlabel('c'); ylabel('P/D');
hold off;

%% Risk-free rate and expected value of the risky return.

h = exp((1-gamma)*xm)*exp(.5*(1-gamma+gamma^2)*stdy^2);
H = beta*h/(1-beta*h);

% Expected value of the risky return.
ER = (H+1)/H*exp(xm+.5*stdy^2);
% Risk-free rate.
Rf = exp(gamma*xm-.5*(gamma*stdy)^2)/beta;

%% Simulation part.
% simulate from this transition probability
% initialize random number generators to something truly normal
randn('state',sum(100*clock));
rand('twister',sum(100*clock));

NSim  = 10000;
xtsim = zeros(NSim,1);
Rsim = ones(NSim,1);
stsim = zeros(NSim,1); % keeps track of state in which one is
Pc    = cumsum(P,2); % sums the various columns of trans prob
idx   = round(N/2); % current state, just any state
for SimCtr=1:NSim
    u=rand(1,1);
    Rsim(SimCtr,1) = Rsim(SimCtr,1)/f(idx);
    for j=1:N % find the position of the next state
        if Pc(idx,j)>u
            idx=j;
            break
        end
    end
    stsim(SimCtr,1) = idx;
    xtsim(SimCtr,1) = xj(idx);
    Rsim(SimCtr,1) = Rsim(SimCtr,1)*(1+f(idx))*exp(xj(idx));
end

figure(2);
hold on;
plot(xtsim); xlabel('t'); ylabel('c_t');
title('Simulation of consumption rate');
hold off;

figure(3);
hold on;
plot(xtsim); xlabel('t'); ylabel('R_t');
title('Simulation of stock return');
hold off;