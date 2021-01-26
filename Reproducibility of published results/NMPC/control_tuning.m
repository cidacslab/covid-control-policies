function [uc,Psic,S,E,Ia,Is,H,U,R,D,Nw]=control_tuning(Q,Usig)

close all
clc

%% ============ Load data ============
load 'input_data/covidbahia_CIDACS'
load 'input_data/Psi_Ba13set'
load 'input_data/parameters'
load 'input_data/string_ba'
parameter=x;
coviddata=covidbahia;
data_in = datetime(2020,03,06);
strig=strig/100;

%% ============ Filter Psi ============
global psir costd 
costd = inf;

npd = length(Psi);
Nfo = 7;            % Filter order

for k = 1:npd
    if (k-Nfo) < 0
        Psif(k) = 1/k*sum(Psi(1:k));
    else
        Psif(k) = 1/Nfo*sum(Psi(k-Nfo+1:k));
    end
end

%% ============ Initial Conditions ============

Psi_dado    = Psif(1:size(coviddata,1));
D0          = 0;
N0          = 14930634;
R0          = 0;
H0          = 0;
U0          = 0;
Is0         = 2.015439771376298e-06;
Ia0         = 1.8028646508967777e-06;
E0          = 1.7639153732952095e-06;
S0          = (1-Is0-Ia0-E0);
Nw0         = 0;
Model_0     = [Psi_dado(1),N0*[S0,E0,Ia0,Is0,H0,U0,R0,D0,Nw0]];


Psic = Model_0(1);    
S = Model_0(2);
E = Model_0(3);
Ia= Model_0(4);
Is= Model_0(5);
H = Model_0(6);
U = Model_0(7);
R = Model_0(8);
D = Model_0(9);
Nw = Model_0(10);


%% ============ Control and Simulation Settings ============
Hc=21;
tsim=33*7;
Ts=7;

Gain = [34.4/33.7 46.0/40.9 42.2/49.2 41.2/49.2 39.1/40.9 39.0/39.1];
N_obd = max(Gain);

%% ============ Simulation High Adherence for differents Q ============
for k=1:Ts:tsim

    Th = linspace(k-1,k+Hc-1,Hc+1);

    l = fix(k/13)+1;

    if l>size(parameter,2)
        l=size(parameter,2);
    end

    par = parameter(:,l);

    cost=covid_opt(x,Th,Model_0,N0,Usig,par,N_obd,Q);

    xv = psir*ones(1,length(Usig));
    [minValue,closestIndex] = min(abs(Usig-xv));

    x0 = Usig(closestIndex);
    uc(k:k+Ts) = Usig(closestIndex);

    t = linspace(k-1,k+Ts-1,Ts+1);
    [t,y] = ode45(@(t,y)covid_odes_2(t,y,uc(k),N0,N_obd),t,Model_0);

Psic=[Psic;y(2:end,1)];
S=[S;y(2:end,2)];
E=[E;y(2:end,3)];
Ia=[Ia;y(2:end,4)];
Is=[Is;y(2:end,5)];
H=[H;y(2:end,6)];
U=[U;y(2:end,7)];
R=[R;y(2:end,8)];
D=[D;y(2:end,9)];
Nw=[Nw;y(2:end,10)];

N0=N0-D(end);
Model_0 = [Psic(end),S(end),E(end),Ia(end),Is(end),H(end),U(end),R(end),D(end),Nw(end)]';
costd=inf;
end

close all

end