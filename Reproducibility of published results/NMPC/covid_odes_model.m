function dydt = covid_odes_model(t,y,Psif,N0,x)
N=1;
Psi=Psif(floor(t+1));

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

if t<20.178
beta=2.1317;    %beta=1.3987731952032998;
elseif (t>=28.178-8)&&(t< 72.94)
beta=1.7645;    %0.9614724422279308;    
elseif (t>=72.94)&&(t< 148)
beta=1.1281;    %0.6657552424857321; 
else 
    beta=1;
end

S=y(1);
E=y(2);
Ia=y(3);
Is=y(4);
H=y(5);
U=y(6);
R=y(7);
D=y(8);
Nw=y(9);

dSdt = -(1-Psi)*beta*S*(Is+delta*Ia)/N;
dEdt = (1-Psi)*beta*S*(Is+delta*Ia)/N - kappa*E;
dIadt = (1-p)*kappa*E - gamma_a*Ia;
dIsdt = p*kappa*E - gamma_s*Is;
dHdt = h*qsi*gamma_s*Is + (1-mi_u+ome_u*mi_u)*gamma_u*U - gamma_h*H;
dUdt = h*(1-qsi)*gamma_s*Is + ome_h*gamma_h*H - gamma_u*U;
dRdt = gamma_a*Ia + (1-h)*gamma_s*Is + (1-mi_h)*(1-ome_h)*gamma_h*H;
dDdt = (1-ome_h)*mi_h*gamma_h*H + (1-ome_u)*mi_u*gamma_u*U;
dNwdt = p * kappa * E ;

dydt = [dSdt; dEdt; dIadt; dIsdt; dHdt; dUdt; dRdt; dDdt; dNwdt];

end