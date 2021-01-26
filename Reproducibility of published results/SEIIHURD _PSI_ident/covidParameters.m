function parameters=covidParameters(Dados,Model_0,N0,t,Psi)

%x0=[1 ; 1/4; 1/3.5; 1/4; 0.18; 0.13; 0.4; 0.53; 0.06];
%lb=[0 ; 1/6;  1/3.7; 1/5; 1/12; 1/12; 0.4; 0.5 ; 0.05];
%ub=[3 ; 1/3; 1/3.24; 1/3;  1/4;  1/3; 0.5; 0.99; 0.25];    
    
%[ Beta, h, xi, mi_u, gamma_u, gamma_h, p, omega_u, omega_h, mi_h, delta]
x0=[1 ; 0.06; 0.53; 0.4; 0.13; 0.18; 0.2; 0.29; 0.14; 0.15; 0.31 ];
lb=[0 ; 0.05; 0.5; 0.4; 1/12; 1/12; 0.13; 0.1; 0.1; 0.1; 0.01 ];
ub=[3 ; 0.25; 0.99; 0.5; 1/3; 1/4; 0.5; 0.3; 0.3; 0.2; 0.75 ];

    % Otimização com gradiente MSE
    options = optimoptions('fmincon','Algorithm','interior-point');
    options = optimoptions(options,'Display', 'off');
    options = optimoptions(options,'PlotFcn', { @optimplotfval });
    [x,fval,exitflag] = ...
    fmincon(@(x)covid_opt(x,t,Model_0,Dados,N0,Psi),x0,[],[],[],[],lb,ub,[],options);
    
    parameters = x;
    
end


%% =====================================================================
%  Functions para calcular o custo e simular o modelo
%  =====================================================================

function mse_total=covid_opt(x,t,Model_0,Dados,N0,Psi)

% Calcula o valor 
 
[t,S,E,Ia,Is,H,U,R,D,Nw]=solver_model(x,t,Model_0,N0,Psi);


% Inicia valores dos erros
mse_h=0;
mse_u=0;
mse_d=0;
mse_c=0;

for i=1:size(Dados,1)
    
     mse_c=mse_c+(Dados(i,1)-Nw(i))^2;
     mse_h=mse_h+(Dados(i,2)-H(i))^2;
     mse_u=mse_u+(Dados(i,3)-U(i))^2;
     mse_d=mse_d+(Dados(i,4)-D(i))^2;
     
end

mse_total=1*mse_c+5*mse_d+35*mse_h+35*mse_u;

end

%% Solver Model SIRD
function [t,S,E,Ia,Is,H,U,R,D,Nw]=solver_model(x,t,Model_0,N,Psi)

% Resolve o modelo de no tempo t
options=odeset('NonNegative',(1:9));
[t,y] = ode45(@(t,y)covid_odes(t,y,x,N,Psi),t,Model_0,options);

S=y(:,1);
E=y(:,2);
Ia=y(:,3);
Is=y(:,4);
H=y(:,5);
U=y(:,6);
R=y(:,7);
D=y(:,8);
Nw=y(:,9);

end

%% COVID ODEs
function dydt = covid_odes(t,y,x,N,Psi)

Psi=Psi(floor(t));

% Variáveis de Decisão %[ Beta, h, xi, mi_u, gamma_u, gamma_h, p, omega_u, omega_h, mi_h, delta
beta=x(1);
h=x(2);
qsi=x(3);
mi_u=x(4);
gamma_u=x(5);
gamma_h=x(6);
p=x(7);
ome_u=x(8);
ome_h=x(9);
mi_h=x(10);
delta=x(11);


% Parâmetros fixos
   %qsi = 0.53;
   gamma_a = 1/3.5;
   gamma_s = 1/4;
   k = 1/4;
   %mi_h = 0.15;
   %ome_h = 0.14;
   %ome_u = 0.29;
   %delta = 0.31;
   %h = 0.28;
   %p = 0.2;
   %gamma_h = 0.14;
   %gamma_u = 0.14;
   %mi_u = 0.4;

% Variáveis do modelo
S=y(1);
E=y(2);
Ia=y(3);
Is=y(4);
H=y(5);
U=y(6);
R=y(7);
D=y(8);
Nw=y(9);

% Resolve as EDO`s 
dSdt = -(1-Psi)*beta*S*(Is+delta*Ia)/N;
dEdt = (1-Psi)*beta*S*(Is+delta*Ia)/N - k*E;
dIadt = (1-p)*k*E - gamma_a*Ia;
dIsdt = p*k*E - gamma_s*Is;
dHdt = h*qsi*gamma_s*Is + (1-mi_u+ome_u*mi_u)*gamma_u*U - gamma_h*H;
dUdt = h*(1-qsi)*gamma_s*Is + ome_h*gamma_h*H - gamma_u*U;
dRdt = gamma_a*Ia + (1-h)*gamma_s*Is + (1-mi_h)*(1-ome_h)*gamma_h*H;
dDdt = (1-ome_h)*mi_h*gamma_h*H + (1-ome_u)*mi_u*gamma_u*U;
dNwdt = p*k*E ;

dydt = [dSdt; dEdt; dIadt; dIsdt; dHdt; dUdt; dRdt; dDdt; dNwdt];

end