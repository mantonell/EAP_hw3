% Illustration how to discretize an AR(1)
% x_t = b + a x_{t-1} + eps_t
% eps_t N(0,sigeps^2)
close all, clc, clear all

% determine the boundaries of the domain
alpha   = 0.00001;
mu      = 0;
sigma   = 1;
k_alpha = norminv(1-alpha/2,mu,sigma)
a       = 0.8;
b       = 2;
sigeps  = 3;
stdy    = sigeps/sqrt(1-a^2);
yL      = k_alpha*stdy;

xm      = b/(1-a);

N=60;
Bk=linspace(-yL,yL,N+1) % grid of boundaries
Bk=xm+Bk;
Bkub=Bk; % will be grid with boundaries but where extremes are +/- infty
Bkub(1)=-Inf
Bkub(end)=Inf
xj    = (Bk(2:N+1)+Bk(1:N))/2;

% trace a little plot to show that everything is nice and cool
hold on
plot([Bk(1)-1,Bk(end)+1],[0,0],'-k')
for j=1:N+1
    plot([Bk(j);Bk(j)],[-0.5;0.5],'-r')
end
for j=1:N
    plot(xj(j),0,'ok')
end
axis([Bk(1)-1,Bk(end)+1,-1,1])
hold off

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
s= AM\bv % same as inv(AM)*bv


% and verify again
figure()
hold on
plot(xj,s,'ok')
I=normpdf(xj,b/(1-a),stdy);
plot(xj,I/sum(I),'+m')
hold off

% simulate from this transition probability
% initialize random number generators to something truly normal
randn('state',sum(100*clock))
rand('twister',sum(100*clock))

NSim  = 10000;
xtsim = zeros(NSim,1);
stsim = zeros(NSim,1); % keeps track of state in which one is
Pc    = cumsum(P,2); % sums the various columns of trans prob
idx   = round(N/2); % current state, just any state
for SimCtr=1:NSim
    u=rand(1,1);
    for j=1:N % find the position of the next state
        if Pc(idx,j)>u
            idx=j;
            break
        end
    end
    stsim(SimCtr,1) = idx;
    xtsim(SimCtr,1) = xj(idx);
end

plot(xtsim)
m=mean(xtsim);
s=std(xtsim);

fprintf('Mean, theoretical and empirical %12.6f %12.6f\n',b/(1-a),m)
fprintf('Std , theoretical and empirical %12.6f %12.6f\n',stdy,s)


