clc
clear all
close all

%% ============ Load Data ============
load 'input_data/covidbahia_CIDACS'
load 'input_data/Psi_Ba13set'
coviddata=covidbahia;

Cases_dado      = coviddata(:,1);
R_dado          = coviddata(:,7);
D_dado          = coviddata(:,5);
H_dado          = coviddata(:,2);
U_dado          = coviddata(:,3);
data_in         = datetime(2020,03,06);
Nw_cases        = diff(Cases_dado);
Nw_deaths       = diff(D_dado);

%% ============ Initial Condition ============
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
Model_0     = N0*[S0,E0,Ia0,Is0,H0,U0,R0,D0,Nw0];

%% ============ Identification Settings ============

t_pred  = 25;
t_opt   = 13;                                   
tf      = size(coviddata,1)-t_pred;           

npd = length(Psi);
Nfo = 8;            % Filter order
for k = 1:npd
    if (k-Nfo) < 0
        Psif(k) = 1/k*sum(Psi(1:k));
    else
        Psif(k) = 1/Nfo*sum(Psi(k-Nfo+1:k));
    end
end

%% ============ Initiate Vectors ============
S       = Model_0(1,1);
E       = Model_0(1,2);
Ia      = Model_0(1,3);
Is      = Model_0(1,4);
H       = Model_0(1,5);
U       = Model_0(1,6);
R       = Model_0(1,7);
D       = Model_0(1,8);
Nw      = Model_0(1,9);
t_p     = 1; 
l = 1;

%% ============ Star Identification ============
for i = 1:t_opt:tf
    
    if i+t_opt+1>tf
        dif=tf-i;
        t_opt=dif;
    end

    t        = linspace(i,i+t_opt,t_opt+1);
    Dados    = [Cases_dado(i:i+t_opt),H_dado(i:i+t_opt),U_dado(i:i+t_opt),D_dado(i:i+t_opt)];
    
    parameters=covidParameters(Dados,Model_0,N0,t,Psif);
      
    x(:,l) = parameters;
    
    [t_m,S_m,E_m,Ia_m,Is_m,H_m,U_m,R_m,D_m,Nw_m]=solver_model(x(:,l),t,Model_0,N0,Psif);
    
    t_p = [t_p;t(2:end)'];
    S = [S;S_m(2:end)]; 
    E = [E;E_m(2:end)]; 
    Ia = [Ia;Ia_m(2:end)]; 
    Is = [Is;Is_m(2:end)]; 
    H = [H;H_m(2:end)]; 
    U = [U;U_m(2:end)]; 
    R = [R;R_m(2:end)]; 
    D = [D;D_m(2:end)];
    Nw = [Nw;Nw_m(2:end)];  
    Model = [S,E,Ia,Is,H,U,R,D,Nw];
    
    beta_t(i:i+t_opt)=x(1,l);
    h_t(i:i+t_opt)=x(2,l);
    qsi_t(i:i+t_opt)=x(3,l);
    mi_u_t(i:i+t_opt)=x(4,l);
    gamma_u_t(i:i+t_opt)=x(5,l);
    gamma_h_t(i:i+t_opt)=x(6,l);
    p_t(i:i+t_opt)=x(7,l);
    ome_u_t(i:i+t_opt)=x(8,l);
    ome_h_t(i:i+t_opt)=x(9,l);
    mi_h_t(i:i+t_opt)=x(10,l);
    delta_t(i:i+t_opt)=x(11,l);
    
    l=l+1;
    
    Model_0         = [S(end),E(end),Ia(end),Is(end),H(end),U(end),R(end),D(end),Nw(end)];
end

close all

%% ============ Start Prediction ============

[t_m,S_m,E_m,Ia_m,Is_m,H_m,U_m,R_m,D_m,Nw_m]=solver_model(x(:,end),tf:1:tf+t_pred,Model_0,N0,Psif);

t_p(tf+1:tf+t_pred) = t_m(2:end);
S(tf+1:tf+t_pred) = S_m(2:end); 
E(tf+1:tf+t_pred) = E_m(2:end); 
Ia(tf+1:tf+t_pred) = Ia_m(2:end); 
Is(tf+1:tf+t_pred) = Is_m(2:end); 
H(tf+1:tf+t_pred) = H_m(2:end); 
U(tf+1:tf+t_pred) = U_m(2:end); 
R(tf+1:tf+t_pred) = R_m(2:end); 
D(tf+1:tf+t_pred) = D_m(2:end);
Nw(tf+1:tf+t_pred) = Nw_m(2:end);
    
beta_t(tf+1:tf+t_pred)=x(1,l-1);
h_t(tf+1:tf+t_pred)=x(2,l-1);
qsi_t(tf+1:tf+t_pred)=x(3,l-1);
mi_u_t(tf+1:tf+t_pred)=x(4,l-1);
gamma_u_t(tf+1:tf+t_pred)=x(5,l-1);
gamma_h_t(tf+1:tf+t_pred)=x(6,l-1);
p_t(tf+1:tf+t_pred)=x(7,l-1);
ome_u_t(tf+1:tf+t_pred)=x(8,l-1);
ome_h_t(tf+1:tf+t_pred)=x(9,l-1);
mi_h_t(tf+1:tf+t_pred)=x(10,l-1);
delta_t(tf+1:tf+t_pred)=x(11,l-1);

%% ============ Data and parameters ============

t_dado=1:1:length(Cases_dado);
t_dado = data_in + caldays(0:(length(t_dado)-1));

t_p = data_in + caldays(0:(length(t_p)-1));

% [ Beta, h, xi, mi_u, gamma_u, gamma_h, p, omega_u, omega_h, mi_h, delta]
beta_ori=[2.1317*ones(1,20),1.7645*ones(1,73-20),1.1281*ones(1,148-73),ones(1,(size(beta_t,2))-148)];
param_ori = [beta_ori;ones(1,(size(beta_t,2))).*[0.06; 0.53; 0.4; 0.13; 0.18; 0.2; 0.29; 0.14; 0.15; 0.31]];
param_opt = [beta_t;h_t;qsi_t;mi_u_t;gamma_u_t;gamma_h_t;p_t;ome_u_t;ome_h_t;mi_h_t;delta_t];
G_t = param_opt./param_ori;
G_w = x(2:end,:)./[0.06; 0.53; 0.4; 0.13; 0.18; 0.2; 0.29; 0.14; 0.15; 0.31];

Nw_d=diff(D);
Nw_i=diff(Nw);

%% ============ Model Erros ============

r_c = Cases_dado-Nw(1:length(Cases_dado));
r_nwi = (Nw_cases)-Nw_i(1:length(Nw_cases));
r_nwd = (Nw_deaths)-Nw_d(1:length(Nw_deaths));
r_d = D_dado-D(1:length(D_dado));
r_h = H_dado-H(1:length(H_dado));
r_u = U_dado-U(1:length(U_dado));

SSE_c = norm(r_c).^2;
SSE_nwi = norm(r_nwi).^2;
SSE_nwd = norm(r_nwd).^2;
SSE_d = norm(r_d).^2;
SSE_h = norm(r_h).^2;
SSE_u = norm(r_u).^2;

SST_c = norm(Cases_dado-mean(Nw(1:length(Cases_dado))))^2;
SST_nwi = norm(Nw_cases-mean(Nw_i(1:length(Nw_cases))))^2;
SST_nwd = norm(Nw_deaths-mean(Nw_d(1:length(Nw_deaths))))^2;
SST_d = norm(D_dado-mean(D(1:length(D_dado))))^2;
SST_h = norm(H_dado-mean(H(1:length(H_dado))))^2;
SST_u = norm(U_dado-mean(U(1:length(U_dado))))^2;

% R2 - correlation factor
R2_c = 1 - SSE_c/SST_c;
R2_nwi = 1 - SSE_nwi/SST_nwi;
R2_nwd = 1 - SSE_nwd/SST_nwd;
R2_d = 1 - SSE_d/SST_d;
R2_h = 1 - SSE_h/SST_h;
R2_u = 1 - SSE_u/SST_u;

% Error in %
r_c = (abs(r_c)./Cases_dado)*100;
r_nwi = (abs(r_nwi)./(Nw_cases))*100;
r_nwd = (abs(r_nwd)./Nw_deaths)*100;
r_d = (abs(r_d)./D_dado)*100;
r_h = (abs(r_h)./H_dado)*100;
r_u = (abs(r_u)./U_dado)*100;
    
%% ============ Plot Curves ============

% figure(1)
% plot(t_p,G_t(1,:),'m','LineWidth',2)
% hold on
% plot(t_p,G_t(2,:),'r','LineWidth',2)
% plot(t_p,G_t(3,:),'k','LineWidth',2)
% plot(t_p,G_t(4,:),'b','LineWidth',2)
% plot(t_p,G_t(5,:),'g','LineWidth',2)
% plot(t_p,G_t(6,:),'y','LineWidth',2)
% plot(t_p,G_t(7,:),'r--','LineWidth',2)
% plot(t_p,G_t(8,:),'k--','LineWidth',2)
% plot(t_p,G_t(9,:),'b--','LineWidth',2)
% plot(t_p,G_t(10,:),'g--','LineWidth',2)
% plot(t_p,G_t(11,:),'y--','LineWidth',2)
% xlabel('Time [days]')
% ylabel('Gains')
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% grid on
% legend('g_1(t)','g_2(t)','g_3(t)','g_4(t)','g_5(t)','g_6(t)','g_7(t)','g_8(t)','g_9(t)','g_{10}(t)','g_{11}(t)')

figure(2)
subplot(2,2,3)
H_min = H*0.95;
H_max = H*1.05;
shadedplot(t_p,H_min',H_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
plot(t_p,H,'b','LineWidth',3);
plot(t_dado(1:tf),H_dado(1:tf),'k.','MarkerSize',10);plot(t_dado(tf+1:end),H_dado(tf+1:end),'r.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(end)])
ylabel('Clinical beds occupancy')

figure(2)
subplot(2,2,4)
U_min = U*0.95;
U_max = U*1.05;
shadedplot(t_p,U_min',U_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
plot(t_p,U,'b','LineWidth',3);
plot(t_dado(1:tf),U_dado(1:tf),'k.','MarkerSize',10);plot(t_dado(tf+1:end),U_dado(tf+1:end),'r.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(end)])
ylabel('ICU beds occupancy')

figure(2)
subplot(2,2,1)
Ic_min = (Nw)*0.95;
Ic_max = (Nw)*1.05;
shadedplot(t_p,Ic_min',Ic_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
plot(t_p,Nw,'b','LineWidth',3);
plot(t_dado(1:tf),Cases_dado(1:tf),'k.','MarkerSize',10);plot(t_dado(tf+1:end),Cases_dado(tf+1:end),'r.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(end)])
ylabel('Cumulative infected cases')

figure(2)
subplot(2,2,2)
D_min = D*0.95;
D_max = D*1.05;
shadedplot(t_p,D_min',D_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
plot(t_p,D,'b','LineWidth',3);
plot(t_dado(1:tf),D_dado(1:tf),'k.','MarkerSize',10);plot(t_dado(tf+1:end),D_dado(tf+1:end),'r.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(end)])
ylabel('Fatal cases')

% figure(3)
% Nw_i_min = Nw_i*0.8;
% Nw_i_max = Nw_i*1.2;
% shadedplot(t_p(2:end),Nw_i_min',Nw_i_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p(2:end),Nw_i,'b','LineWidth',3);
% p4=plot(t_dado(2:end),Nw_cases,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(2) t_dado(end)])
% xlabel('Time [days]')
% ylabel('New Cases')

% figure(4)
% Nw_d_min = Nw_d*0.8;
% Nw_d_max = Nw_d*1.2;
% shadedplot(t_p(2:end),Nw_d_min',Nw_d_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p(2:end),Nw_d,'b','LineWidth',3);
% p4=plot(t_dado(2:end),Nw_deaths,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(2) t_dado(end)])
% xlabel('Time [days]')
% ylabel('New Fatal Cases')

figure(5)
plot(t_p,Psi,'k--','LineWidth',1);
hold on
plot(t_p,Psif,'k','LineWidth',3);
title('Isolation Index')
set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
grid on
xlim([t_dado(2) t_dado(end)])
ylabel('\psi (t)')

% figure(6)
% gamma_a = 1/3.5;
% gamma_s = 1/4;
% k = 1/4;
% Ro_t = (S/N0)'.*(beta_t.*(1-Psif).*p_t./gamma_s + beta_t.*(1-Psif).*(1-p_t).*delta_t/gamma_a);
% plot(t_p,Ro_t,'LineWidth',3)
% hold on
% grid on
% xlabel('Time [days]')
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% ylabel('R_t(t)')
% xlim([t_p(1) t_p(end)])
% ylim([0 5])
% title('Effective reproduction number')

lim10 = 10*ones(length(t_dado),1);
figure(7)
plot(t_dado(2:end),r_nwi,'g*',t_dado,r_d,'b',t_dado(2:end),r_nwd,'k+')
hold on
plot(t_dado,r_h,'k-',t_dado,r_u,'b--',t_dado,r_c,'r+')
plot(lim10,'--r')
xlabel('dias')
ylabel('Erro Data - Modelo [%]')
xlim([t_dado(1) t_dado(end)])
ylim([0 50])
set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
legend(sprintf('New Cases:R^{2} = %0.4f',R2_nwi),sprintf('D:R^{2} = %0.4f',R2_d),sprintf('New Deaths:R^{2} = %0.4f',R2_nwd),sprintf('H:R^{2} = %0.4f',R2_h),sprintf('U:R^{2} = %0.4f',R2_u),sprintf('Cases:R^{2} = %0.4f',R2_c)) 
title('SEIIHURD+\Psi and R^2')

%% ============ Save data ============

    saveas(figure(7),[pwd '/Figures/seiihurd_opt_error.fig']);
    saveas(figure(2),[pwd '/Figures/seiihurd_opt_painel.fig']);
    saveas(figure(5),[pwd '/Figures/psif.fig']);
    save('output_data/seiihurd_ident_data.mat')
    save('output_data/parameters.mat','x')

%% ============ Solver Model SEIIHURD ============
function [t,S,E,Ia,Is,H,U,R,D,Nw]=solver_model(x,t,Model_0,N,Psi)

% Resolve o modelo de no tempo t
options=odeset('NonNegative',(1:8));
[t,y] = ode45(@(t,y)covid_odes(t,y,x,N,Psi),t,Model_0',options);

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

%% ============ COVID ODE ============
function dydt = covid_odes(t,y,x,N,Psi)

Psi=Psi(floor(t));

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

gamma_a = 1/3.5;
gamma_s = 1/4;
k = 1/4;

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
