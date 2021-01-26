clear
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
Nfo = 8;            % Filter order

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

Model_0_ho = Model_0;
Model_0_lo = Model_0;

Psic_ho = Model_0(1);    
S_ho = Model_0(2);
E_ho = Model_0(3);
Ia_ho= Model_0(4);
Is_ho= Model_0(5);
H_ho = Model_0(6);
U_ho = Model_0(7);
R_ho = Model_0(8);
D_ho = Model_0(9);
Nw_ho = Model_0(10);
N0_ho = N0;

Psic_lo = Model_0(1);    
S_lo = Model_0(2);
E_lo = Model_0(3);
Ia_lo= Model_0(4);
Is_lo= Model_0(5);
H_lo = Model_0(6);
U_lo = Model_0(7);
R_lo = Model_0(8);
D_lo = Model_0(9);
Nw_lo = Model_0(10);
N0_lo = N0;
%% ============ Define Scenario ============

% USER CHOICE
SCENARIO = 1;

switch SCENARIO
    case 1 
    Usig = 0.01*[21.62 36.71 38.10 40.87 49.21];
    Q = [8e4 3e4 1e4];
    case 2
    Usig = 0.01*[21.62 36.71 38.10 39.08 40.87 41.50 42.23 43.25 44.84 45.90 47.05 49.21 49.80 51.59 54.56 54.69 55.61 57.94 60.12 62.70];
    Q = [8e4 3e4 1e4];
    case 3
    Usig = 0.01*[21.62 36.71 38.10 40.87 49.21];
    Q = [8e4 3e4 1e4
        5e6 1e6 2e5
        1e4 5e3 1e3
        5e3 1e3 1e3];
        for i=1:size(Q,1)
         [uc(:,i),Psic(:,i),S(:,i),E(:,i),Ia(:,i),Is(:,i),H(:,i),U(:,i),R(:,i),D(:,i),Nw(:,i)]=control_tuning(Q(i,:),Usig);
        end
    figure_Q(uc,Psic,S,E,Ia,Is,H,U,R,D,Nw,Q);
    return
end

%% ============ Control and Simulation Settings ============
Hc=21;
tsim=32*7;
Ts=7;

Gain = [34.4/33.7 46.0/40.9 42.2/49.2 41.2/49.2 39.1/40.9 39.0/39.1];
N_obd_lo = min(Gain);
N_obd_ho = max(Gain);

%% ============ Simulation High Adherence ============

for k=1:Ts:tsim
         
    Th = linspace(k-1,k+Hc-1,Hc+1);
    
    l = fix(k/13)+1;
    
    if l>size(parameter,2)
        l=size(parameter,2);
    end
    
    par = parameter(:,l);
     
    cost=covid_opt(x,Th,Model_0_ho,N0_ho,Usig,par,N_obd_ho,Q);

    xv = psir*ones(1,length(Usig));
    [minValue,closestIndex] = min(abs(Usig-xv));
    
    x0 = Usig(closestIndex);
    uc_ho(k:k+Ts) = Usig(closestIndex);
       
    t = linspace(k-1,k+Ts-1,Ts+1);
    [t,y] = ode45(@(t,y)covid_odes_2(t,y,uc_ho(k),N0_ho,N_obd_ho),t,Model_0_ho');

Psic_ho=[Psic_ho;y(2:end,1)];
S_ho=[S_ho;y(2:end,2)];
E_ho=[E_ho;y(2:end,3)];
Ia_ho=[Ia_ho;y(2:end,4)];
Is_ho=[Is_ho;y(2:end,5)];
H_ho=[H_ho;y(2:end,6)];
U_ho=[U_ho;y(2:end,7)];
R_ho=[R_ho;y(2:end,8)];
D_ho=[D_ho;y(2:end,9)];
Nw_ho=[Nw_ho;y(2:end,10)];

N0_ho=N0_ho-D_ho(end);
Model_0_ho = [Psic_ho(end),S_ho(end),E_ho(end),Ia_ho(end),Is_ho(end),H_ho(end),U_ho(end),R_ho(end),D_ho(end),Nw_ho(end)];
costd=inf;
end


%% ============ Simulation Low Adherence ============
for k=1:Ts:tsim
         
    Th = linspace(k-1,k+Hc-1,Hc+1);
    
    l = fix(k/13)+1;
    
    if l>size(parameter,2)
        l=size(parameter,2);
    end
    
    par = parameter(:,l);
     
    cost=covid_opt(x,Th,Model_0_lo,N0_lo,Usig,par,N_obd_lo,Q);

    xv = psir*ones(1,length(Usig));
    [minValue,closestIndex] = min(abs(Usig-xv));
    
    x0 = Usig(closestIndex);
    uc_lo(k:k+Ts) = Usig(closestIndex);
       
    t = linspace(k-1,k+Ts-1,Ts+1);
    [t,y] = ode45(@(t,y)covid_odes_2(t,y,uc_lo(k),N0_lo,N_obd_lo),t,Model_0_lo');

Psic_lo=[Psic_lo;y(2:end,1)];
S_lo=[S_lo;y(2:end,2)];
E_lo=[E_lo;y(2:end,3)];
Ia_lo=[Ia_lo;y(2:end,4)];
Is_lo=[Is_lo;y(2:end,5)];
H_lo=[H_lo;y(2:end,6)];
U_lo=[U_lo;y(2:end,7)];
R_lo=[R_lo;y(2:end,8)];
D_lo=[D_lo;y(2:end,9)];
Nw_lo=[Nw_lo;y(2:end,10)];

N0_lo=N0_lo-D_lo(end);
Model_0_lo = [Psic_lo(end),S_lo(end),E_lo(end),Ia_lo(end),Is_lo(end),H_lo(end),U_lo(end),R_lo(end),D_lo(end),Nw_lo(end)];
costd=inf;
end

close all



%% ============ Model Simulation =============

Cases_dado      = coviddata(:,1)';
R_dado          = coviddata(:,7)';
D_dado          = coviddata(:,5)';
H_dado          = coviddata(:,2)';
U_dado          = coviddata(:,3)';
Nw_cases        = diff(Cases_dado);
Nw_deaths       = diff(D_dado);

Model_0 =  [S0,E0,Ia0,Is0,H0,U0,R0,D0,Nw0];

for j = size(Psif,2):(tsim)
    if mod((j-173),21) == 0
        Psif(j+1)=max(min(Psif),Psif(end)*0.98);
    else
        Psif(j+1)=Psif(end);
    end
end

[t,y] = ode45(@(t,y)covid_odes_model(t,y,Psif,N0),(0:1:tsim),Model_0');
  
Sm=y(:,1);
Em=y(:,2);
Iam=N0.*y(:,3);
Ism=N0.*y(:,4);
Hm=N0.*y(:,5);
Um=N0.*y(:,6);
Rm=N0.*y(:,7);
Dm=N0.*y(:,8);
Nwm=N0.*y(:,9);    

npts=length(t);

%% ============ DATA TO PLOT AND INDICES ============

nlhos1 = ones(59,1)*466;
nluti1 = ones(59,1)*422;
nlhos2 = ones(npts-59,1)*1610;
nluti2 = ones(npts-59,1)*1210;

M_Psi = mean(Psi);
Var_Psi = sum(Psi-min(Psi(1:size(Psi,1))))/size(Psi,1);

M_Psif = mean(Psif);
Var_Psif = sum((Psif)-min(Psif))/size(Psif,2);
M_Psic_ho = mean(Psic_ho);
Var_Psic_ho = sum(Psic_ho-min(Psic_ho))/size(Psic_ho,1);
M_uc_ho = mean(uc_ho);
Var_uc_ho = sum(uc_ho-min(uc_ho))/size(uc_ho,2);
M_Psic_lo = mean(Psic_lo);
Var_Psic_lo = sum(Psic_lo-min(Psic_lo))/size(Psic_lo,1);
M_uc_lo = mean(uc_lo);
Var_uc_lo = sum(uc_lo-min(uc_lo))/size(uc_lo,2);

statistc_DATA = [M_Psif Var_Psif 
             M_Psic_ho Var_Psic_ho 
             M_uc_ho Var_uc_ho 
             M_Psic_lo Var_Psic_lo 
             M_uc_lo Var_uc_lo];

K_psi =  Psi(end)/strig(end);

strig = [strig;ones(tsim-length(strig)+1,1)*K_psi*strig(end)];

tplotsDATA   = datetime(2020,3,06) + caldays(0:(length(Cases_dado)-1));
tplotsbeforemodel = datetime(2020,3,06) + caldays(0:(npts-1));

Per_cases_ho= 1-Nw_ho(end)/Nwm(end)
Per_D_ho= 1-D_ho(end)/Dm(end)
Per_H_ho= 1-max(H_ho)/max(Hm)
Per_U_ho= 1-max(U_ho)/max(Um)

Per_cases_lo= Nw_lo(end)/Nwm(end)
Per_D_lo= D_lo(end)/Dm(end)

Psi_fev = 0.2845;
%% ============ Plots ============

figure(4)
subplot(2,1,1)
plot(tplotsbeforemodel,(Nw_ho),'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); grid on; hold on;
plot(tplotsbeforemodel,(Nw_lo),'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3); grid on; hold on;
plot(tplotsbeforemodel,(Nwm),'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); grid on; hold on;
plot(tplotsDATA,(Cases_dado),'k.','MarkerSize',10); grid on; hold on;
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('Cumulative cases')
title('(A)', 'FontSize', 20);
legend('High adherence','Low adherence','Validated model','Real data','Orientation','horizontal')
legend boxoff
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(4)
subplot(2,1,2)
plot(tplotsbeforemodel,Dm,'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); 
grid on; hold on;
plot(tplotsbeforemodel,D_ho,'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); 
plot(tplotsbeforemodel,D_lo,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3); 
plot(tplotsDATA,D_dado,'k.','MarkerSize',10); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('Fatal cases')
title('(B)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(2)
subplot(2,2,1)
Nwm_i=diff(Nwm);
Nw_i_ho=diff(Nw_ho);
Nw_i_lo=diff(Nw_lo);
Nw_dado=diff(Cases_dado);
fill([tplotsbeforemodel(npd) tplotsbeforemodel(npts) tplotsbeforemodel(npts) tplotsbeforemodel(npd)], [0 0 1e4 1e4], [0.87, 0.87, 0.87],'FaceAlpha', 0.5, 'EdgeColor','w');
grid off; hold on;
p(1)=plot(tplotsbeforemodel(2:end),(Nw_i_ho),'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); 
p(2)=plot(tplotsbeforemodel(2:end),(Nw_i_lo),'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3); 
p(3)=plot(tplotsbeforemodel(2:end),(Nwm_i),'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); 
p(4)=plot(tplotsDATA(2:end),(Nw_dado),'k.','MarkerSize',10); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('New cases')
title('(A)', 'FontSize', 20);
legend([p(1) p(2) p(3) p(4)],'High adherence','Low adherence','Validated model','Real data','Orientation','vertical','Location', 'best')
legend boxoff
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(2)
subplot(2,2,2)
NwDm_i=diff(Dm);
NwD_i_ho=diff(D_ho);
NwD_i_lo=diff(D_lo);
NwD_dado=diff(D_dado);
fill([tplotsbeforemodel(npd) tplotsbeforemodel(npts) tplotsbeforemodel(npts) tplotsbeforemodel(npd)], [0 0 1.5e2 1.5e2], [0.87, 0.87, 0.87],'FaceAlpha', 0.5, 'EdgeColor','w');
grid off; hold on;
plot(tplotsbeforemodel(2:end),(NwDm_i),'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); 
plot(tplotsbeforemodel(2:end),(NwD_i_ho),'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); 
plot(tplotsbeforemodel(2:end),(NwD_i_lo),'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3);
plot(tplotsDATA(2:end),(NwD_dado),'k.','MarkerSize',10); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('New deaths')
title('(B)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(2)
subplot(2,2,3)
fill([tplotsbeforemodel(npd) tplotsbeforemodel(npts) tplotsbeforemodel(npts) tplotsbeforemodel(npd)], [0 0 2.5e3 2.5e3], [0.87, 0.87, 0.87],'FaceAlpha', 0.5, 'EdgeColor','w');
grid off; hold on;
plot(tplotsbeforemodel,Hm,'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); 
plot(tplotsbeforemodel,H_ho,'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); 
plot(tplotsbeforemodel,H_lo,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3); 
plot(tplotsDATA,H_dado,'k.','MarkerSize',10); 
plot(tplotsbeforemodel,[nlhos1; nlhos2],'k-.','LineWidth',1); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('Clinical beds occupancy')
title('(C)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(2)
subplot(2,2,4)
fill([tplotsbeforemodel(npd) tplotsbeforemodel(npts) tplotsbeforemodel(npts) tplotsbeforemodel(npd)], [0 0 2e3 2e3], [0.87, 0.87, 0.87],'FaceAlpha', 0.5, 'EdgeColor','w');
grid off; hold on;
plot(tplotsbeforemodel,Um,'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); 
plot(tplotsbeforemodel,U_ho,'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); 
plot(tplotsbeforemodel,U_lo,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3); 
plot(tplotsDATA,U_dado,'k.','MarkerSize',10); 
plot(tplotsbeforemodel,[nluti1; nluti2],'k-.','LineWidth',1); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('ICU beds occupancy')
title('(D)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(1)
subplot(2,1,1)
plot(tplotsbeforemodel,uc_ho,'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); 
hold on; grid on;
plot(tplotsbeforemodel,uc_lo,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3); 
plot(tplotsbeforemodel(1:npd),strig(1:npd),'Color',[0, 0.4470, 0.7410],'LineStyle','-','LineWidth',3);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('Stringency index')
title('(A)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])
ylim([0 1])


figure(1)
subplot(2,1,2)
plot(tplotsbeforemodel,Psic_ho,'Color',[0.9300, 0.4350, 0.0980],'LineWidth',3); 
grid on; hold on;
plot(tplotsbeforemodel,Psic_lo,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',3); 
plot(tplotsbeforemodel(1:npd),Psif(1:npd),'Color',[0, 0.4470, 0.7410],'LineStyle','-','LineWidth',3); plot(tplotsbeforemodel(npd+1:end),Psif(npd+1:end),'Color',[0, 0.4470, 0.7410],'LineStyle',':','LineWidth',3);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('SMRI')
title('(B)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])
ylim([0 1])
legend('High adherence','Low adherence','Observed index','Assumed index','Orientation','horizontal','Location','best')
legend boxoff

figure(3)
plot(tplotsbeforemodel(1:length(Psi)),strig(1:length(Psi)),'Color',[0, 0.4470, 0.7410],'LineStyle','-','LineWidth',3); grid on; hold on;
plot(tplotsbeforemodel(1:length(Psi)),Psi,'LineWidth',3); 
plot(tplotsbeforemodel(1:length(Psi)),Psi_fev*ones(1,size(Psi,1)),'k--','LineWidth',3);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('Stringency and SMRI')
xlim([tplotsbeforemodel(1) tplotsbeforemodel(length(Psi))])
legend('Stringency index','Observed SMRI','SMRI baseline','Orientation','horizontal')
legend boxoff

%% ============ Save data ============

saveas(figure(3),[pwd '/Figures/fig1.fig']);

switch SCENARIO
    case 1
    saveas(figure(1),[pwd '/Figures/input_scen1.fig']);
    saveas(figure(2),[pwd '/Figures/output_scen1.fig']);
    save('output_data/scenario_1.mat')
    case 2
    saveas(figure(1),[pwd '/Figures/input_scen2.fig']);
    saveas(figure(2),[pwd '/Figures/output_scen2.fig']);
    save('output_data/scenario_2.mat')
end
