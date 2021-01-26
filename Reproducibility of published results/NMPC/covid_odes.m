function dydt = covid_odes(t,y,x,N0,par,N_obd)
uc = x;

%% Parameters
nominal = 0;

if nominal==0 
beta=par(1);
h=par(2);
qsi=par(3);
mi_u=par(4);
gamma_u=par(5);
gamma_h=par(6);
p=par(7);
ome_u=par(8);
ome_h=par(9);
mi_h=par(10);
delta=par(11);
gamma_a = 1/3.5;
gamma_s = 1/4;
kappa = 1/4;
  
else

kappa = 1/4;
gamma_a = 1/3.5;
gamma_s = 1/4;
gamma_h = 0.1776863042741119;
gamma_u = 0.13342706158133355;
mi_u = 0.4;
qsi = 0.53;
h = 0.06287612499693644;
mi_h = 0.15;
ome_h = 0.14;
ome_u = 0.29;
delta = 0.30906304338495505;
p = 0.2;

% Beta
if t<20.178
beta=2.1317;    %beta=1.3987731952032998;
elseif (t>=28.178-8)&&(t< 96.94-24)
beta=1.7645;    %0.9614724422279308;    
elseif (t>=72.94)&&(t< 148)
beta=1.1281;    %0.6657552424857321; 
else 
    beta=1;
end

end
%% ODE's

Psi=y(1);
S=y(2);
E=y(3);
Ia=y(4);
Is=y(5);
H=y(6);
U=y(7);
R=y(8);
D=y(9);
Nw=y(10);
N=N0-D;

if t < 59
   nlhos = 466;
   nluti = 422;
else
    nlhos = 1610;
    nluti = 1210;
end

tau = 0.4;
K = N_obd;

dPsidt = tau*K*uc - tau*Psi;
dSdt = -(1-Psi)*beta*S*(Is+delta*Ia)/N;
dEdt = (1-Psi)*beta*S*(Is+delta*Ia)/N - kappa*E;
dIadt = (1-p)*kappa*E - gamma_a*Ia;
dIsdt = p*kappa*E - gamma_s*Is;
dHdt = h*qsi*gamma_s*Is + (1-mi_u+ome_u*mi_u)*gamma_u*U - gamma_h*H;
dUdt = h*(1-qsi)*gamma_s*Is + ome_h*gamma_h*H - gamma_u*U;
dRdt = gamma_a*Ia + (1-h)*gamma_s*Is + (1-mi_h)*(1-ome_h)*gamma_h*H;
dDdt = (1-ome_h)*mi_h*gamma_h*H + (1-ome_u)*mi_u*gamma_u*U;
dNwdt = p * kappa * E ;

dydt = [dPsidt; dSdt; dEdt; dIadt; dIsdt; dHdt; dUdt; dRdt; dDdt; dNwdt];
end