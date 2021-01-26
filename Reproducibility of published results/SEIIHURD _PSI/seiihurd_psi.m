close all
clear all
clc


%% ============ Load Data ============
load 'input_data/covidbahia_CIDACS'
load 'input_data/Psi_Ba13set'

npd = length(Psi);
Nfo = 8;            % Filter order
for k = 1:npd
    if (k-Nfo) < 0
        Psif(k) = 1/k*sum(Psi(1:k));
    else
        Psif(k) = 1/Nfo*sum(Psi(k-Nfo+1:k));
    end
end

coviddata=covidbahia;

Cases_dado      = coviddata(:,1)';
R_dado          = coviddata(:,7)';
D_dado          = coviddata(:,5)';
H_dado          = coviddata(:,2)';
U_dado          = coviddata(:,3)';
data_in = datetime(2020,03,06);
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
Model_0 =  [S0,E0,Ia0,Is0,H0,U0,R0,D0,Nw0];

%% ============ Model ============
    
tsim = length(U_dado)-1; 
    
[t,y] = ode45(@(t,y)covid_odes(t,y,Psif),0:1:tsim,Model_0');

Sm=y(:,1);
Em=y(:,2);
Iam=N0.*y(:,3);
Ism=N0.*y(:,4);
Hm=N0.*y(:,5);
Um=N0.*y(:,6);
Rm=N0.*y(:,7);
Dm=N0.*y(:,8);
Nwm=N0.*y(:,9);    

%% ============ DATA ============

t_p = data_in + caldays(0:(length(U_dado)-1));
t_dado=1:1:length(H_dado);
t_dado = data_in + caldays(0:(length(t_dado)-1));

Nw_i=diff(Nwm);
Nw_d=diff(Dm);

r_nwi = (Nw_cases)-Nw_i(1:length(Nw_cases))';
r_nwd = (Nw_deaths)-Nw_d(1:length(Nw_deaths))';
r_h = H_dado-Hm(1:length(H_dado))';
r_u = U_dado-Um(1:length(U_dado))';

SSE_nwi = norm(r_nwi).^2;
SSE_nwd = norm(r_nwd).^2;
SSE_h = norm(r_h).^2;
SSE_u = norm(r_u).^2;

SST_nwi = norm(Nw_cases-mean(Nw_i(1:length(Nw_cases))))^2;
SST_nwd = norm(Nw_deaths-mean(Nw_d(1:length(Nw_deaths))))^2;
SST_h = norm(H_dado-mean(Hm(1:length(H_dado))))^2;
SST_u = norm(U_dado-mean(Um(1:length(U_dado))))^2;

% R2 - correlation factor
R2_nwi = 1 - SSE_nwi/SST_nwi;
R2_nwd = 1 - SSE_nwd/SST_nwd;
R2_h = 1 - SSE_h/SST_h;
R2_u = 1 - SSE_u/SST_u;

% Error in %
r_nwi = (abs(r_nwi)./(Nw_cases))*100;
r_nwd = (abs(r_nwd)./Nw_deaths)*100;
r_h = (abs(r_h)./H_dado)*100;
r_u = (abs(r_u)./U_dado)*100;

%% ============ Plot ============

figure(1)
subplot(2,2,3)
H_min = Hm*0.95;
H_max = Hm*1.05;
shadedplot(t_p,H_min',H_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p,Hm,'b','LineWidth',3);
p4=plot(t_dado,H_dado,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(tsim)])
ylabel('Clinical beds occupancy')
      
figure(1)
subplot(2,2,4)
U_min = Um*0.95;
U_max = Um*1.05;
shadedplot(t_p,U_min',U_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p,Um,'b','LineWidth',3);
p4=plot(t_dado,U_dado,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(1) t_dado(tsim)])
ylabel('ICU beds occupancy')
    
% figure(3)
% Ic_min = (Nwm)*0.8;
% Ic_max = (Nwm)*1.2;
% shadedplot(t_p,Ic_min',Ic_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p,Nwm,'b','LineWidth',3);
% p4=plot(t_dado,Cases_dado,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(1) t_dado(tsim)])
% xlabel('Time [days]')
% ylabel('Cumulative Infected Cases')

% figure(4)
% D_min = Dm*0.8;
% D_max = Dm*1.2;
% shadedplot(t_p,D_min',D_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p,Dm,'b','LineWidth',3);
% p4=plot(t_dado,D_dado,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(1) t_dado(tsim)])
% xlabel('Time [days]')
% ylabel('Fatal Cases')

figure(1)
subplot(2,2,1)
Nw_i_min = Nw_i*0.95;
Nw_i_max = Nw_i*1.05;
shadedplot(t_p(2:end),Nw_i_min',Nw_i_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p(2:end),Nw_i,'b','LineWidth',3);
p4=plot(t_dado(2:end),Nw_cases,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(2) t_dado(tsim)])
ylabel('New cases')

figure(1)
subplot(2,2,2)
Nw_d_min = Nw_d*0.95;
Nw_d_max = Nw_d*1.05;
shadedplot(t_p(2:end),Nw_d_min',Nw_d_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
hold on
p3=plot(t_p(2:end),Nw_d,'b','LineWidth',3);
p4=plot(t_dado(2:end),Nw_deaths,'k.','MarkerSize',10);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
xlim([t_dado(2) t_dado(tsim)])
ylabel('New fatal cases')

% figure(6)
% R_min = Rm*0.8;
% R_max = Rm*1.2;
% shadedplot(t_p,R_min',R_max',[0.9 0.9 0.95],[0.9 0.9 0.95]);
% hold on
% p3=plot(t_p,Rm,'b','LineWidth',3);
% p4=plot(t_dado,R_dado,'k*','LineWidth',2);
% legend([p3 p4],{'SEIIHURD+\Psi Model','Official Data - SESAB'})
% set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
% xlim([t_dado(1) t_dado(tsim)])
% xlabel('Time [days]')
% ylabel('Recovered Cases')

lim10 = 10*ones(length(t_dado),1);
figure(10)
plot(t_dado(2:end),r_nwi,'b',t_dado(2:end),r_nwd,'k+')
hold on
plot(t_dado,r_h,'k-',t_dado,r_u,'b--')
plot(lim10,'--r')
xlabel('dias')
ylabel('Erro Data - Modelo [%]')
xlim([t_dado(2) t_dado(end)])
ylim([0 50])
set(gca,'FontSize',30,'LineWidth',2,'FontWeight','Bold')
legend(sprintf('New Cases:R^{2} = %0.4f',R2_nwi),sprintf('New Deaths:R^{2} = %0.4f',R2_nwd),sprintf('H:R^{2} = %0.4f',R2_h),sprintf('U:R^{2} = %0.4f',R2_u)) 
title('SEIIHURD+\Psi and R^2')

%% ============ Save Data ============

saveas(figure(1),[pwd '/Figures/seiihurd_painel.fig']);
saveas(figure(10),[pwd '/Figures/output_scen1.fig']);
save('output_data/seiihurd_psi_data.mat')

%% ============ Model simulation ============
function dydt = covid_odes(t,y,Psif)
N=1;
psi=Psif(floor(t+1));

k = 1/4;
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

dSdt = -(1-psi)*beta*S*(Is+delta*Ia)/N;
dEdt = (1-psi)*beta*S*(Is+delta*Ia)/N - k*E;
dIadt = (1-p)*k*E - gamma_a*Ia;
dIsdt = p*k*E - gamma_s*Is;
dHdt = h*qsi*gamma_s*Is + (1-mi_u+ome_u*mi_u)*gamma_u*U - gamma_h*H;
dUdt = h*(1-qsi)*gamma_s*Is + ome_h*gamma_h*H - gamma_u*U;
dRdt = gamma_a*Ia + (1-h)*gamma_s*Is + (1-mi_h)*(1-ome_h)*gamma_h*H;
dDdt = (1-ome_h)*mi_h*gamma_h*H + (1-ome_u)*mi_u*gamma_u*U;
dNwdt = p* k * E ;

dydt = [dSdt; dEdt; dIadt; dIsdt; dHdt; dUdt; dRdt; dDdt; dNwdt];

end