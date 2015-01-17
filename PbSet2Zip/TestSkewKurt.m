function [equity_premium]=TestSkewKurt(SS,KK)
%% Problem 3.
%clear; clc;

% Parameters of the model.
a     = .1; 
B     = .7; 
sig   = .1;
gamma = 2;
beta  = .9;

abfreud=zeros(100,2);
load -ascii coefffreud4.txt;
abfreud(:,2)=coefffreud4;
Nj=50;
xw=gauss(Nj,abfreud);

uj = xw(:,1); % abscissas
wj = xw(:,2); % weights with respect to w(x)=exp(-x^4)


%% Definition of the f function.
fnam = 'SKKUBench.txt';
[ag, bg, sg, kg] = textread(fnam,'%f %f %f %f');
InitStruct.ag=ag;
InitStruct.bg=bg;
InitStruct.sg=sg;
InitStruct.kg=kg;
% Skewness and kurtosis.
SkC=SS;
KuC=KK;
% Lambda-parameters.
[L0C,L1C,L2C,L3C,L4C] = Canon7_Gautschi(0,1,SkC,KuC,xw,InitStruct);

%clear('xw','abfreud','coefffreud4');

f = @(x) exp(L0C+L1C*x+L2C*x^2+L3C*x^3+L4C*x^4);

%% Steady-state distribution.

% Determine h, the steady-state distribution.
LDA = zeros(Nj,Nj);
U = (-L4C)^(.25)*B;
for j=1:Nj
    for k=1:Nj
        X=f(uj(j)/U-B*uj(k)/U)*exp(uj(k)^4);
        LDA(j,k)=X*wj(k)/U;
    end
end

ImII=eye(Nj)-LDA;
II1=1/U*wj.*exp(uj.^4);

A=[II1'; ImII(2:Nj,:)];
bb = zeros(Nj,1); 
bb(1)=1;

h = A\bb; % The ergodic distribution

%% Construction of the Markov chain.

% Transition probability matrix.
II=zeros(Nj,Nj);
for j=1:Nj
    for k=1:Nj
        II(k,j)=f(uj(j)/(-L4C)^.25-B/(-L4C)^.25*uj(k))*...
            exp(uj(j)^4)*wj(j)/(-L4C)^.25;
    end
end
IIs  = sum(II,2);
TrPr = II./repmat(IIs,1,Nj); % transition probability

%% Simulation.

TrPrs=cumsum(TrPr,2);
xj=uj/(-L4C)^.25;

NSim=1000;
xt=zeros(NSim,1); % corresponds to the variable from which one simulates
idx=zeros(NSim,1); % corresponds to the index of the state 
idx(1)=floor(Nj/2); % starts in the middle of nowhere
xt(1) = xj(idx(1));
for j=2:NSim
    u=rand(1,1);
    z = TrPrs(idx(j-1),:); % a line in the transition prob
    for L=Nj-1:(-1):1
        if u>=z(L)
            jj=L+1;
            break
        elseif (u<z(1))&&(L==1)
            jj=1;
        end
    end
   idx(j) = jj;
   xt(j)  = xj(jj);
end

% Verification of dynamics.
T=NSim;
y=xt(2:T);
xx=[ones(T-1,1) xt(1:T-1)];
res = ols(y,xx);
disp('regression xt on xt(-1)')
prt(res)

%% Steady-state price-dividend process.

% Consumption growth rate.
m=a/(1-B); 
cj=m+sig*xj;

% Finds P/D process.
psi=zeros(Nj,Nj);
for j=1:Nj
    for k=1:Nj
        psi(j,k) = beta*exp((1-gamma)*cj(k));
    end
end
b=psi.*TrPr*ones(Nj,1);
H=(eye(Nj)-psi.*TrPr)\b;

%% Simulation of returns.
R=(H(idx(2:end))+1)./H(idx(1:end-1)).*exp(cj(idx(2:end)));

% Summary statistics
mean(R)
std(R)
min(R)
max(R)
skewness(R)
kurtosis(R)

% Equity premium
Rf=1/(beta*mean(exp(-gamma*cj(idx))));
equity_premium=mean(R)-Rf;
end