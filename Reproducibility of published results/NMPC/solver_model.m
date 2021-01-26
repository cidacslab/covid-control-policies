%% Solver Model SIRD
function [t,Psi,S,E,Ia,Is,H,U,R,D,Nw]=solver_model(x,t,Model_0,N0,par,N_obd)


% Resolve o modelo de no tempo t
%options=odeset('NonNegative',(1:8));
[t,y] = ode45(@(t,y)covid_odes(t,y,x,N0,par,N_obd),t,Model_0);

Psi=y(:,1);
S=y(:,2);
E=y(:,3);
Ia=y(:,4);
Is=y(:,5);
H=y(:,6);
U=y(:,7);
R=y(:,8);
D=y(:,9);
Nw=y(:,10);

end
