%% PROJEVOL
 
% Developed by: Marissa Dudek (1)
% (1) University of North Carolina at Chapel Hill, Geological Sciences
 
% Clean workspace
clc; format long; 
set(0,'DefaultFigureWindowStyle','docked')
counter = 0;

% ----------------------------------------------------------------------------

% Build color code for projectile characteristics
color_rhop = zeros(n.projectiles,3);
width_Lp = zeros(n.projectiles,1);
for k = 1:n.projectiles
    if FLAVRS(k,5) == rhop.ch(1) % chondritic / blue
        color_rhop(k,:) = [0 0.4470 0.7410];
    elseif FLAVRS(k,5) == rhop.ir(1) % iron / red
        color_rhop(k,:) = [0.6350 0.0780 0.1840];
    elseif FLAVRS(k,5) == rhop.si(1) % stony iron / orange
        color_rhop(k,:) = [0.8500 0.3250 0.0980];
    elseif FLAVRS(k,5) == rhop.st(1) % stony / green
        color_rhop(k,:) = [0.4660 0.6740 0.1880];
    end
    if FLAVRS(k,2) >= 10000 % greater than 10 km in diameter
        width_Lp(k,:) = 3;
    elseif FLAVRS(k,2) >= 1000 % greater than 1 km in diameter
        width_Lp(k,:) = 2;
    elseif FLAVRS(k,2) >= 100 % greater than 100 m in diameter
        width_Lp(k,:) = 1;
    elseif FLAVRS(k,2) >= 10 % greater than 10 m in diameter
        width_Lp(k,:) = 0.5;
    elseif FLAVRS(k,2) >= 0 % greater than 0 m in diameter
        width_Lp(k,:) = 0.25;
    end
    if FLAVRS(k,5) == rhop.ch(1) % chondritic / blue
        color_rhop(k,:) = [0 0.4470 0.7410];
    elseif FLAVRS(k,5) == rhop.ir(1) % iron / red
        color_rhop(k,:) = [0.6350 0.0780 0.1840];
    elseif FLAVRS(k,5) == rhop.si(1) % stony iron / orange
        color_rhop(k,:) = [0.8500 0.3250 0.0980];
    elseif FLAVRS(k,5) == rhop.st(1) % stony / green
        color_rhop(k,:) = [0.4660 0.6740 0.1880];
    end
end
 
% ----------------------------------------------------------------------------
% AEM Manuscript Figure 8

% Import MATALB workspace
load('Workspace_AEM_TheoreticalPopulation_AblationFragmentation.mat'); % Theoretical projectile popultion (est.2024)
for k = 1:n.projectiles
%     if FLAVRS(k,11) == 45
%         if FLAVRS(k,4) == 20000
            loglog(FLAVRS(k,7),craterM(k,3),'o','Color',color_rhop(k,:)); hold on;
            loglog(FLAVRS(k,7),craterE(k,3),'x','Color',color_rhop(k,:)); hold on;
            loglog(FLAVRS(k,7),craterV(k,3),'*','Color',color_rhop(k,:));
            legend({'Mars','Earth','Venus'},'Location','northwest'); grid on;
            ylabel('Crater Diameter (km)'); xlabel('Initial Projectile Mass (kg)'); ylim([.1 1100])
            title('Projectile Diameter vs. Crater Diameter')
%         end
%     end
end

% % Import MATALB workspace
% load('Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.mat'); % Terrestrial analog popultion (est.2024)
% for k = 1:n.projectiles
%     loglog(FLAVRS(k,7),craterM(k,3),'o','Color','k'); hold on;
%     loglog(FLAVRS(k,7),craterE(k,3),'x','Color','k'); hold on;
%     loglog(FLAVRS(k,7),craterV(k,3),'*','Color','k');
%     legend({'Mars','Earth','Venus'},'Location','northwest'); grid on;
%     ylabel('Crater Diameter (km)'); xlabel('Initial Projectile Mass (kg)'); ylim([.1 1100])
%     title('Projectile Diameter vs. Crater Diameter')
% end

% ----------------------------------------------------------------------------
% AEM Manuscript Figure 9

% Import MATALB workspace
% load('Workspace_AEM_TheoreticalPopulation_AblationFragmentation.mat'); % Terrestrial analog popultion (est.2024)

% for k = 1:n.projectiles  
% subplot(3,2,1);plot(Mvel(k,:)/1000,Malt(k,:)/1000,'Color',color_rhop(k,:),'LineWidth',width_Lp(k,:)); hold on; ylabel('Altitude (km)'); xlabel('Velocity (km/s)'); xlim([0 40]);  ylim([0 150]);  
% subplot(3,2,2);plot(Mmassr(k,:),Malt(k,:)/1000,'Color',color_rhop(k,:),'LineWidth',width_Lp(k,:)); hold on; xlabel('Mass Remaining (%)'); xlim([0 100]);  ylim([0 150]); 
% subplot(3,2,3);plot(Evel(k,:)/1000,Ealt(k,:)/1000,'Color',color_rhop(k,:),'LineWidth',width_Lp(k,:)); hold on; ylabel('Altitude (km)'); xlabel('Velocity (km/s)'); xlim([0 40]); ylim([0 150]); 
% subplot(3,2,4);plot(Emassr(k,:),Ealt(k,:)/1000,'Color',color_rhop(k,:),'LineWidth',width_Lp(k,:)); hold on; xlabel('Mass Remaining (%)'); xlim([0 100]);  ylim([0 150]); 
% subplot(3,2,5);plot(Vvel(k,:)/1000,Valt(k,:)/1000,'Color',color_rhop(k,:),'LineWidth',width_Lp(k,:)); hold on; ylabel('Altitude (km)'); xlabel('Velocity (km/s)'); xlim([0 40]); ylim([0 150]); 
% subplot(3,2,6);plot(Vmassr(k,:),Valt(k,:)/1000,'Color',color_rhop(k,:),'LineWidth',width_Lp(k,:)); hold on; xlabel('Mass Remaining (%)'); xlim([0 100]); ylim([0 150]); 
% end

