function figure_Q(uc,Psic,S,E,Ia,Is,H,U,R,D,Nw,Q)

%% ============ Load data ============
load 'input_data/covidbahia_CIDACS'
load 'input_data/Psi_Ba13set'
load 'input_data/parameters'
load 'input_data/string_ba'
parameter=x;
coviddata=covidbahia;
data_in = datetime(2020,03,06);
strig=strig/100;

%% ============ Initial Conditions ============

npd = length(Psi);
Nfo = 8;            % Filter order

for k = 1:npd
    if (k-Nfo) < 0
        Psif(k) = 1/k*sum(Psi(1:k));
    else
        Psif(k) = 1/Nfo*sum(Psi(k-Nfo+1:k));
    end
end

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


tsim=33*7;
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

K_psi =  Psi(end)/strig(end);

strig = [strig;ones(tsim-length(strig)+1,1)*K_psi*strig(end)];

tplotsDATA   = datetime(2020,3,06) + caldays(0:(length(Cases_dado)-1));
tplotsbeforemodel = datetime(2020,3,06) + caldays(0:(npts-1));
         
%% ============ Plots ============

% figure(4)
% subplot(2,1,1)
% plot(tplotsbeforemodel,(Nwm),'g','LineWidth',3); grid on; hold on;
% plot(tplotsbeforemodel,(Nw),'b','LineWidth',3); grid on; hold on;
% plot(tplotsDATA,(Cases_dado),'k.','MarkerSize',10); grid on; hold on;
% set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
% xlabel('Time [days]')
% ylabel('Cumulative Cases')
% xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])
% 
% figure(4)
% subplot(2,1,2)
% plot(tplotsbeforemodel,Dm,'g','LineWidth',3); grid on; hold on;
% plot(tplotsbeforemodel,D,'b','LineWidth',3); grid on; hold on;
% plot(tplotsDATA,D_dado,'k.','MarkerSize',10); grid on; hold on;
% set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
% xlabel('Time [days]')
% ylabel('Fatal Cases')
% xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(2)
subplot(2,2,3)
fill([tplotsbeforemodel(npd) tplotsbeforemodel(npts) tplotsbeforemodel(npts) tplotsbeforemodel(npd)], [0 0 2.5e3 2.5e3], [0.87, 0.87, 0.87],'FaceAlpha', 0.5, 'EdgeColor','w');
grid off; hold on;
set(gca,'ColorOrderIndex',1)
plot(tplotsbeforemodel,Hm,'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); 
for i=1:size(Q,1)
plot(tplotsbeforemodel,H(:,i),'LineWidth',3); 
end
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
set(gca,'ColorOrderIndex',1)
plot(tplotsbeforemodel,Um,'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3);
for i=1:size(Q,1)
plot(tplotsbeforemodel,U(:,i),'LineWidth',3); 
end
plot(tplotsDATA,U_dado,'k.','MarkerSize',10); 
plot(tplotsbeforemodel,[nluti1; nluti2],'k-.','LineWidth',1); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('ICU beds occupancy')
title('(D)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(2)
subplot(2,2,1)
Nwm_i=diff(Nwm);
Nw_i=diff(Nw);
Nw_dado=diff(Cases_dado);
fill([tplotsbeforemodel(npd) tplotsbeforemodel(npts) tplotsbeforemodel(npts) tplotsbeforemodel(npd)], [0 0 1e4 1e4], [0.87, 0.87, 0.87],'FaceAlpha', 0.5, 'EdgeColor','w');
grid off; hold on;
set(gca,'ColorOrderIndex',1)
plot(tplotsbeforemodel(2:end),(Nwm_i),'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3); 
for i=1:size(Q,1)
plot(tplotsbeforemodel(2:end),Nw_i(:,i),'LineWidth',3);
end
plot(tplotsDATA(2:end),(Nw_dado),'k.','MarkerSize',10); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('New cases')
title('(A)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(2)
subplot(2,2,2)
NwDm_i=diff(Dm);
NwD_i=diff(D);
NwD_dado=diff(D_dado);
fill([tplotsbeforemodel(npd) tplotsbeforemodel(npts) tplotsbeforemodel(npts) tplotsbeforemodel(npd)], [0 0 1.5e2 1.5e2], [0.87, 0.87, 0.87],'FaceAlpha', 0.5, 'EdgeColor','w');
grid off; hold on;
set(gca,'ColorOrderIndex',1)
plot(tplotsbeforemodel(2:end),(NwDm_i),'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',3);
for i=1:size(Q,1)
plot(tplotsbeforemodel(2:end),(NwD_i(:,i)),'LineWidth',3);
end
plot(tplotsDATA(2:end),(NwD_dado),'k.','MarkerSize',10); 
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('New deaths')
title('(B)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])

figure(1)
subplot(2,1,2)
plot(tplotsbeforemodel(1:npd),Psif(1:npd),'Color',[0, 0.4470, 0.7410],'LineStyle','-','LineWidth',3);grid on; hold on; 
for i=1:size(Q,1)
plot(tplotsbeforemodel,Psic(:,i),'LineWidth',3);
end
plot(tplotsbeforemodel(npd+1:end),Psif(npd+1:end),'Color',[0, 0.4470, 0.7410],'LineStyle',':','LineWidth',3);
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('SMRI')
title('(B)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])
ylim([0 1])

figure(1)
subplot(2,1,1)
plot(tplotsbeforemodel(1:length(Psi)),strig(1:length(Psi)),'Color',[0, 0.4470, 0.7410],'LineStyle','-','LineWidth',3); grid on; hold on;
for i=1:size(Q,1)
plot(tplotsbeforemodel,uc(:,i),'LineWidth',3); 
end
set(gca,'FontSize',20,'LineWidth',2,'FontWeight','Bold')
ylabel('Stringency index')
title('(A)', 'FontSize', 20);
xlim([tplotsbeforemodel(1) tplotsbeforemodel(end)])
ylim([0 1])


%% ============ Save data ============


     %saveas(figure(1),[pwd '/Figures/input_scen3.fig']);
     %saveas(figure(2),[pwd '/Figures/output_scen3.fig']);
     save('output_data/scenario_3.mat')



end