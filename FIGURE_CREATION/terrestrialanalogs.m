%% TERRESTRIAL ANALOGS & PANCAKE FACTORS
 
% Developed by: Marissa Dudek (1)
% (1) University of North Carolina at Chapel Hill, Geological Sciences
 
% Clean workspace
clc; format long; 
set(0,'DefaultFigureWindowStyle','docked')
counter = 0;

% ----------------------------------------------------------------------------
%% Run 1st time to population ablation only atmospheric travel first
% ... First run the model with the second section commented out

% load('Workspace_AEM_TerrestrialAnalogs_AblationOnly.mat'); % Terrestrial Analog Pop
% FLAVRS = FLAVRS;
% 
% k = 1;
% subplot(4,2,1)
% plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'-','LineWidth',2); hold on; 
% xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
% legend('Ablation Only','Ablation & Breakup','location','northwest')
% title('Change in Velocity');
% subplot(4,2,2)
% plot(Emassr(k,:).',Ealt(k,:).'/1000.','-','LineWidth',2); hold on; 
% xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');
% title('Change in Mass');
% k = 2;
% subplot(4,2,3)
% plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'-','LineWidth',2); hold on; 
% xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
% subplot(4,2,4)
% plot(Emassr(k,:).',Ealt(k,:).'/1000.','-','LineWidth',2); hold on; 
% xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');
% k = 3;
% subplot(4,2,5)
% plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'-','LineWidth',2); hold on; 
% xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
% subplot(4,2,6)
% plot(Emassr(k,:).',Ealt(k,:).'/1000.','-','LineWidth',2); hold on; 
% xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');
% k = 4;
% subplot(4,2,7)
% plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'-','LineWidth',2); hold on; 
% xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
% subplot(4,2,8)
% plot(Emassr(k,:).',Ealt(k,:).'/1000.','-','LineWidth',2); hold on; 
% xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');

% ----------------------------------------------------------------------------
%% Run 2nd time to population ablation and fragmentation atmospheric travel first
% ... Second run the model with the first section commented out

load('Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.mat'); % Terrestrial Analog Pop
FLAVRS = FLAVRS;

k = 1;
subplot(4,2,1)
plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'--','LineWidth',2); hold on; 
xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
legend('Ablation','Ablation & Breakup','location','northwest')
title('Change in Velocity');
subplot(4,2,2)
plot(Emassr(k,:).',Ealt(k,:).'/1000.','--','LineWidth',2); hold on; 
xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');
title('Change in Mass');
k = 2;
subplot(4,2,3)
plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'--','LineWidth',2); hold on; 
xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
subplot(4,2,4)
plot(Emassr(k,:).',Ealt(k,:).'/1000.','--','LineWidth',2); hold on; 
xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');
k = 3;
subplot(4,2,5)
plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'--','LineWidth',2); hold on; 
xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
subplot(4,2,6)
plot(Emassr(k,:).',Ealt(k,:).'/1000.','--','LineWidth',2); hold on; 
xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');
k = 4;
subplot(4,2,7)
plot(Evel(k,:).'/1000,Ealt(k,:).'/1000,'--','LineWidth',2); hold on; 
xlim([0 22.5]); ylim([0 50]); xlabel('Velocity (km/s)'); ylabel('Altitude (km)');
subplot(4,2,8)
plot(Emassr(k,:).',Ealt(k,:).'/1000.','--','LineWidth',2); hold on; 
xlim([0 105]); ylim([0 50]); xlabel('Mass Remaining (%)');
 
% ----------------------------------------------------------------------------

plot(Epanfact.',Ealt.'/1000,'LineWidth',2); hold on; xlabel('Pancake Factor (Lz/L0)'); xlim([0 25]); ylim([0 200]); ylim([0 100]);



