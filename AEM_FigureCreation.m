%% Title:Atmospheric Effects Model (Figure Creation)

% Last update: 2024-01-05

% Developed by: Marissa Dudek (1)
% (1) University of North Carolina at Chapel Hill, Geological Sciences

% Clean workspace
close all; clc; format long; 
counter = 0;

%% Plot 

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
end
% ,'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)
% if FLAVRS(k,5) == rhop.ir(1); end

%% Atmospheric profile plots

% Crater stats
% ... histogram
    % figure(1); hist(craterL(:,3), 10); xlim([0 (5*10^6)]); ylim([0 (4*10^4)]);
    % figure(2); hist(craterM(:,3), 10); xlim([0 (5*10^6)]); ylim([0 (4*10^4)]);
    % figure(3); hist(craterE(:,3), 10); xlim([0 (5*10^6)]); ylim([0 (4*10^4)]);
    % figure(4); hist(craterV(:,3), 10); xlim([0 (5*10^6)]); ylim([0 (4*10^4)]);
% ... line plots
% figure(1)
% subplot(1,4,1); h1 = histogram(craterL(:,3)/1000, 5); xlabel('Crater Diameter (km)'); ylabel('Frequency'); xlim([0 250])
% subplot(1,4,2); h2 = histogram(craterM(:,3)/1000, 5); xlabel('Crater Diameter (km)'); ylabel('Frequency'); xlim([0 250])
% subplot(1,4,3); h3 = histogram(craterE(:,3)/1000, 5); xlabel('Crater Diameter (km)'); ylabel('Frequency'); xlim([0 250]) 
% subplot(1,4,4); h4 = histogram(craterV(:,3)/1000, 5); xlabel('Crater Diameter (km)'); ylabel('Frequency'); xlim([0 250])
% xlabel('Crater Diameter (km)'); ylabel('Frequency')

% Crater diameter and depth
% figure(1);
% for k = 1:n.projectiles
% loglog(craterL(:,3),craterL(:,4),'*','Color','k'); hold on;
% loglog(craterM(:,3),craterM(:,4),'*','Color',[0.6350 0.0780 0.1840]); hold on;
% loglog(craterE(:,3),craterE(:,4),'*','Color',[0 0.4470 0.7410]); hold on;
% loglog(craterV(:,3),craterV(:,4),'*','Color',[0.8500 0.3250 0.0980]); hold on;
% legend({'Lunar','Mars','Earth','Venus'},'Location','northwest'); grid on;
% xlabel('Final Crater Diameter (km)'); ylabel('Final Crater Depth (km)');
% title('Lunar Crater Diameter vs. Depth')
% end

% for k = 1:n.projectiles
% figure(2); 
% loglog(FLAVRS(k,2)/1000,craterL(k,3),'square','LineWidth',4,'Color',color_rhop(k,:),'MarkerSize',3); hold on; 
% loglog(FLAVRS(k,2)/1000,craterM(k,3),'o','Color',color_rhop(k,:)); hold on; 
% loglog(FLAVRS(k,2)/1000,craterE(k,3),'x','Color',color_rhop(k,:)); hold on; 
% loglog(FLAVRS(k,2)/1000,craterV(k,3),'*','Color',color_rhop(k,:)); 
% legend({'Lunar','Mars','Earth','Venus'},'Location','northwest'); grid on;
% ylabel('Crater Diameter (km)'); xlabel('Iniital Projectile Diameter (km)'); ylim([.1 1100])
% title('Projectile Diameter vs. Crater Diameter')


% figure(3); 
% loglog(FLAVRS(:,2)/1000,craterL(:,3),'Color','k'); hold on; 
% loglog(FLAVRS(:,2)/1000,craterM(:,3),'Color',[0.6350 0.0780 0.1840]); hold on; 
% loglog(FLAVRS(:,2)/1000,craterE(:,3),'Color',[0 0.4470 0.7410]); hold on; 
% loglog(FLAVRS(:,2)/1000,craterV(:,3),'Color',[0.8500 0.3250 0.0980]); 
% legend({'Lunar','Mars','Earth','Venus'},'Location','northwest'); grid on;
% ylabel('Crater Diameter (km)'); xlabel('Iniital Projectile Diameter (km)'); ylim([.1 1100])

%%
% figure(1)
% Venus
    msx = V_atmos.temp; msy = V_atmos.alt/1000;
    lx = (800/16000)*(V_atmos.press/1000); ly = V_atmos.alt/1000;
    evx = (800/150)*V_atmos.rho; evy = V_atmos.alt/1000;
    figure(1)
    a = axes('units','normalized','position',[.1 .35 .7 .6],'xtick',0:50:800);
    plot(msx,msy,'color','k','LineWidth',2); hold on; xlim([0 800]); 
    plot(lx,ly,'color','#0072b2','LineWidth',2); hold on
    plot(evx,evy,'color','#d55e00','LineWidth',2); grid on; grid minor;
    title('Venus'); ylabel('Altitude [km]'); 
    yline(50,'--'); yline(75,'--'); yline(110,'--');
    text(610,25,'Troposphere'); text(500,62.5,'Lower Stratosphere'); text(500,92.5,'Upper Statosphere'); text(500,145,'Thermosphere');
    xlabel(a,'Temperature [K]'); hold on;
    b = axes('units','normalized','position',[.1 .25 .7 0.000001],'xlim',[0 16000],'color','#0072b2','xtick',0:2000:16000,'XColor','#0072b2');
    xlabel(b,'Pressure [kPa]'); b.FontSize = 10;
    c = axes('units','normalized','position',[.1 .15 .7 0.000001],'xlim',[0 150],'color','#d55e00','xtick',0:10:150,'XColor','#d55e00');
    xlabel(c,'Density [kg/m^3]'); c.FontSize = 10;
    
% Earth
    msx = E_atmos.temp; msy = E_atmos.alt/1000;
    lx = (800/160)*(E_atmos.press/1000); ly = E_atmos.alt/1000;
    evx = (800/2.5)*E_atmos.rho; evy = E_atmos.alt/1000;
    figure(2)
    a = axes('units','normalized','position',[.1 .35 .7 .6],'xtick',0:50:800);
    plot(msx,msy,'color','k','LineWidth',2); hold on; xlim([0 800]);
    plot(lx,ly,'color','#0072b2','LineWidth',2); hold on
    plot(evx,evy,'color','#d55e00','LineWidth',2); grid on; grid minor;
    title('Earth'); ylabel('Altitude [km]'); 
    yline(12,'--'); yline(50,'--'); yline(85,'--');
    text(610,6,'Troposphere'); text(550,30,'Stratosphere'); text(550,67.5,'Mesosphere'); text(550,145,'Thermosphere');
    xlabel(a,'Temperature [K]'); hold on;
    b = axes('units','normalized','position',[.1 .25 .7 0.000001],'xlim',[0 160],'color','#0072b2','xtick',0:20:160,'XColor','#0072b2');
    xlabel(b,'Pressure [kPa]'); b.FontSize = 10;
    c = axes('units','normalized','position',[.1 .15 .7 0.000001],'xlim',[0 2.5],'color','#d55e00','xtick',0:0.2:2.5,'XColor','#d55e00');
    xlabel(c,'Density [kg/m^3]'); c.FontSize = 10;
    
% Mars
    msx = M_atmos.temp; msy = M_atmos.alt/1000;
    lx = (800/2)*(M_atmos.press/1000); ly = M_atmos.alt/1000;
    evx = (800/0.10)*M_atmos.rho; evy = M_atmos.alt/1000;
    figure(3)
    a = axes('units','normalized','position',[.1 .35 .7 .6],'xtick',0:50:800);
    plot(msx,msy,'color','k','LineWidth',2); hold on; xlim([0 800]);
    plot(lx,ly,'color','#0072b2','LineWidth',2); hold on
    plot(evx,evy,'color','#d55e00','LineWidth',2); grid on; grid minor;
    title('Mars'); ylabel('Altitude [km]'); 
    yline(40,'--'); yline(100,'--');
    text(500,25,'Troposphere'); text(500,70,'Stratosphere'); text(500,145,'Thermosphere');
    xlabel(a,'Temperature [K]'); hold on;
    b = axes('units','normalized','position',[.1 .25 .7 0.000001],'xlim',[0 2],'color','#0072b2','xtick',0:0.2:2,'XColor','#0072b2');
    xlabel(b,'Pressure [kPa]'); b.FontSize = 10;
    c = axes('units','normalized','position',[.1 .15 .7 0.000001],'xlim',[0 0.10],'color','#d55e00','xtick',0:0.02:0.10,'XColor','#d55e00');
    xlabel(c,'Density [kg/m^3]'); c.FontSize = 10;

%     subplot(1,2,1)
%     plot(Vvel(k,1:end-1)/1000,Valt(k,1:end-1)/1000, 'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on; 
%     ylabel('Altitude [km]'); xlabel('Velocity [km/s]');
%     subplot(1,2,2)
%     plot(Vmassr(k,1:end-1),Valt(k,1:end-1)/1000, 'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on; 
%     ylabel('Altitude [km]'); xlabel('Mass Remaining [km/s]');

% figure(1)
% % Venus
%     subplot(1,3,1)
%     semilogx(V_atmos.temp,V_atmos.alt/1000,'color','k','LineWidth',2); hold on;
%     semilogx(V_atmos.rho,V_atmos.alt/1000,'color','#0072b2','LineWidth',2); hold on;
%     semilogx(V_atmos.press,V_atmos.alt/1000,'color','#d55e00','LineWidth',2); hold on;
%     xline(0.1,'color',[0.5, 0.5, 0.5],'LineWidth',1.5);
%     yline(50,'--'); yline(75,'--'); yline(110,'--');
%     text(10e-14,25,'Troposphere'); text(10e-14,62.5,'Lower Stratosphere'); text(10e-14,92.5,'Upper Statosphere'); text(10e-14,145,'Thermosphere');
%     title({'Venus'}); ylabel('Altitude [km]'); axis([10e-15 10e8 0 200]);
%     legend({'Temperature (K)','Density (kg/m3)','Pressure (Pa)'},'Location','southoutside'); grid on; grid minor;
% % Earth
%     subplot(1,3,2)
%     semilogx(E_atmos.temp,V_atmos.alt/1000,'color','k','LineWidth',2); hold on;
%     semilogx(E_atmos.rho,V_atmos.alt/1000,'color','#0072b2','LineWidth',2); hold on;
%     semilogx(E_atmos.press,V_atmos.alt/1000,'color','#d55e00','LineWidth',2); hold on;
%     xline(0.1,'color',[0.5, 0.5, 0.5],'LineWidth',1.5);
%     yline(12,'--'); yline(50,'--'); yline(85,'--');
%     text(10e-14,6,'Troposphere'); text(10e-14,25,'Stratosphere'); text(10e-14,67.5,'Mesosphere'); text(10e-14,145,'Thermosphere');
%     title({'Earth'}); ylabel('Altitude [km]'); axis([10e-15 10e8 0 200]); grid on; grid minor;
%     legend({'Temperature (K)','Density (kg/m3)','Pressure (Pa)'},'Location','southoutside');
% % Mars
%     subplot(1,3,3)
%     semilogx(M_atmos.temp,V_atmos.alt/1000,'color','k','LineWidth',2); hold on;
%     semilogx(M_atmos.rho,V_atmos.alt/1000,'color','#0072b2','LineWidth',2); hold on;
%     semilogx(M_atmos.press,V_atmos.alt/1000,'color','#d55e00','LineWidth',2); hold on;
%     xline(0.1,'color',[0.5, 0.5, 0.5],'LineWidth',1.5);
%     yline(40,'--'); yline(100,'--');
%     text(10e-14,25,'Troposphere'); text(10e-14,70,'Stratosphere'); text(10e-14,145,'Thermosphere');
%     title({'Mars'}); ylabel('Altitude [km]'); axis([10e-15 10e8 0 200]); grid on; grid minor;
%     legend({'Temperature (K)','Density (kg/m3)','Pressure (Pa)'},'Location','southoutside'); 
% 
% figure (1)
% for k = 1:n.projectiles 
%     if impactM(k,1) ~= input.z_steps(1)
%         %plot(Vvel(k,1:end-1)/1000,Valt(k,1:end-1)/1000, 'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on; 
%         plot(Evel(k,1:end-1)/1000,Ealt(k,1:end-1)/1000,'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on;  
%         %plot(Mvel(k,1:end-1)/1000,Malt(k,1:end-1)/1000,'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on;
%         ylabel('Altitude (km)'); xlabel('Velocity (km/s)'); xlim([0 20]); ylim([0 200]);
%         %title("Vel:"+flag.vel(1)+"   Ang:"+flag.ang(1)+"   Mass:"+flag.mass(1)+"   Sf:"+flag.shape(1)+"   C_D:"+flag.drag(1)+"   C_L:"+flag.lift(1)+"   C_h_t:"+flag.htc(1)+"   C_h_a:"+flag.hoa(1)+"   Atmos_m_a_x:"+atmos.max(1));
%         %title("Vel:"+flag.vel(1)+"   Ang:"+flag.ang(1)+"   Mass:"+flag.mass(1)+"   Sf:"+flag.shape(1)+"   C_D:"+flag.drag(1)+"   C_L:"+flag.lift(1)+"   C_h_t:"+flag.htc(1)+"   C_h_a:"+flag.hoa(1));
%         title('Velocity')
%         %legend({'Venus','Earth','Mars'},'Location','northwest'); 
%         grid on; grid minor; 
%     end
% end
% 
% figure (2)
% for k = 1:n.projectiles 
%     if impactM(k,1) < input.z_steps(1)
%         %plot(Vmassr(k,1:end-1),Valt(k,1:end-1)/1000, 'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on; 
%         plot(Emassr(k,1:end-1),Ealt(k,1:end-1)/1000,'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on;  
%         %plot(Mmassr(k,1:end-1),Malt(k,1:end-1)/1000,'Color',color_rhop(k,:), 'LineWidth',width_Lp(k)); hold on;
%         ylabel('Altitude (km)'); xlabel('Mass Remaining (%)'); xlim([0 100]); ylim([0 200]);
%         %title("Vel:"+flag.vel(1)+"   Ang:"+flag.ang(1)+"   Mass:"+flag.mass(1)+"   Sf:"+flag.shape(1)+"   C_D:"+flag.drag(1)+"   C_L:"+flag.lift(1)+"   C_h_t:"+flag.htc(1)+"   C_h_a:"+flag.hoa(1)+"   Atmos_m_a_x:"+atmos.max(1));
%         title('Mass Remaining')
%         %legend({'Venus','Earth','Mars'},'Location','northwest'); 
%         grid on; grid minor; 
%     end
% end

%%
% figure(5); boxplot(craterL(:,3)/1000,FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Lunar'); ylim([-10 350]);
% figure(6); boxplot(craterM(:,3)/1000,FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Mars'); ylim([-10 350]);
% figure(7); boxplot(craterE(:,3)/1000,FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Earth'); ylim([-10 350]);
% figure(8); boxplot(craterV(:,3)/1000,FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Venus'); ylim([-10 350]);

% for k = 1:n.projectiles
% figure(9); plot(FLAVRS(k,2)/1000,craterL(k,3)/1000,'.','Color',color_rhop(k,:)); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Lunar'); ylim([-10 350]); hold on;
% figure(10); plot(real(FLAVRS(k,2)/1000),real(craterV(k,3)/1000),'.','Color',color_rhop(k,:)); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Lunar'); ylim([-10 350]); hold on;
% end
% for k = 1:n.projectiles
%     figure(5); semilogx(FLAVRS(k,2)/1000,craterL(k,3)/1000,'*','Color',color_rhop(k,:)); xlabel('Final Crater Diameter (km)'); ylabel('Initial Projectile Diameter  (km)'); grid on; title('Lunar'); ylim([-10 350]); xlim([0.001 100]); hold on;
% end
% for k = 1:n.projectiles
%     figure(6); semilogx(FLAVRS(k,2)/1000,craterM(k,3)/1000,'*','Color',color_rhop(k,:)); xlabel('Final Crater Diameter (km)'); ylabel('Initial Projectile Diameter  (km)'); grid on; title('Mars'); ylim([-10 350]); xlim([0.001 100]); hold on;
% end
% for k = 1:n.projectiles
%     figure(7); semilogx(FLAVRS(k,2)/1000,craterE(k,3)/1000,'*','Color',color_rhop(k,:)); xlabel('Final Crater Diameter (km)'); ylabel('Initial Projectile Diameter  (km)'); grid on; title('Earth'); ylim([-10 350]); xlim([0.001 100]); hold on;
% end

% subplot(2,3,1)
% plot(Vvel(4,:)/1000,Valt(4,:)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2'); hold on; 
% plot(Evel(4,:)/1000,Ealt(4,:)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2'); hold on; 
% plot(Mvel(4,:)/1000,Malt(4,:)/1000,'Color',[0.6350 0.0780 0.1840],'LineWidth',2'); hold on;
% xlim([0 22]); ylim([0 100]); ylabel('Altitude (km)'); xlabel('Velocity (km/s)'); title('Ries Impactor');
% grid on; grid minor;
% legend({'Venus','Earth','Mars'},'Location','northwest'); 

%Subplot of Earth impactors to V/E/M
% subplot(2,3,1)
% plot(Vvel(4,:)/1000,Valt(4,:)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2'); hold on; 
% plot(Evel(4,:)/1000,Ealt(4,:)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2'); hold on; 
% plot(Mvel(4,:)/1000,Malt(4,:)/1000,'Color',[0.6350 0.0780 0.1840],'LineWidth',2'); hold on;
% xlim([0 22]); ylim([0 100]);  xlabel('Velocity (km/s)'); ylabel('Altitude (km)'); title('Ries Impactor');
% legend({'Venus','Earth','Mars'},'Location','northwest'); 
% grid on; grid minor;
% subplot(2,3,2)
% plot(Vvel(3,:)/1000,Valt(3,:)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2'); hold on; 
% plot(Evel(3,:)/1000,Ealt(3,:)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2'); hold on; 
% plot(Mvel(3,:)/1000,Malt(3,:)/1000,'Color',[0.6350 0.0780 0.1840],'LineWidth',2'); hold on;
% xlim([0 22]); ylim([0 100]);  xlabel('Velocity (km/s)'); title('Barringer Impactor');
% grid on; grid minor;
% subplot(2,3,3)
% plot(Vvel(2,:)/1000,Valt(2,:)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2'); hold on; 
% plot(Evel(2,:)/1000,Ealt(2,:)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2'); hold on; 
% plot(Mvel(2,:)/1000,Malt(2,:)/1000,'Color',[0.6350 0.0780 0.1840],'LineWidth',2'); hold on;
% xlim([0 22]); ylim([0 100]); xlabel('Velocity (km/s)'); title('Chelyabinsk Fireball');
% grid on; grid minor;
% subplot(2,3,4)
% plot(Vmassr(4,:),Valt(4,:)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2'); hold on; 
% plot(Emassr(4,:),Ealt(4,:)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2'); hold on; 
% plot(Mmassr(4,:),Malt(4,:)/1000,'Color',[0.6350 0.0780 0.1840],'LineWidth',2'); hold on;
% xlim([0 105]); ylim([0 100]); ylabel('Altitude (km)'); xlabel('Mass Remaining (%)');
% grid on; grid minor;
% subplot(2,3,5)
% plot(Vmassr(3,:),Valt(3,:)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2'); hold on; 
% plot(Emassr(3,:),Ealt(3,:)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2'); hold on; 
% plot(Mmassr(3,:),Malt(3,:)/1000,'Color',[0.6350 0.0780 0.1840],'LineWidth',2'); hold on;
% xlim([0 105]); ylim([0 100]);  xlabel('Mass Remaining (%)');
% grid on; grid minor;
% subplot(2,3,6)
% plot(Vmassr(2,:),Valt(2,:)/1000,'Color',[0.8500 0.3250 0.0980],'LineWidth',2'); hold on; 
% plot(Emassr(2,:),Ealt(2,:)/1000,'Color',[0 0.4470 0.7410],'LineWidth',2'); hold on; 
% plot(Mmassr(2,:),Malt(2,:)/1000,'Color',[0.6350 0.0780 0.1840],'LineWidth',2'); hold on;
% xlim([0 105]); ylim([0 100]); xlabel('Mass Remaining (%)');
% grid on; grid minor; 

% figure(5); boxplot(craterL(:,3),FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Lunar'); ylim([-10 350]);
% figure(6); boxplot(craterM(:,3),FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Mars'); ylim([-10 350]);
% figure(7); boxplot(craterE(:,3),FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Earth'); ylim([-10 350]);
% figure(8); boxplot(craterV(:,3),FLAVRS(:,2)/1000,'Notch','on'); xlabel('Initial Projectile Diameter (km)'); ylabel('Final Crater Diameter (km)'); grid on; title('Venus'); ylim([-10 350]);

%%

% subplot(2,2,1)
% plot(Vmassr.',Valt.'/1000,'Color',[0.8500 0.3250 0.0980]); hold on; 
% plot(Emassr.',Ealt.'/1000,'Color',[0 0.4470 0.7410]); hold on; 
% plot(Mmassr.',Malt.'/1000,'Color',[0.6350 0.0780 0.1840]); hold on;
% xlim([0 105]); ylim([0 100]); ylabel('Altitude (km)'); xlabel('Mass Remaining (%)'); 
% grid on; grid minor;
% 
% subplot(2,2,2)
% plot(Vvel.'/1000,Valt.'/1000,'Color',[0.8500 0.3250 0.0980]); hold on; 
% plot(Evel.'/1000,Ealt.'/1000,'Color',[0 0.4470 0.7410]); hold on; 
% plot(Mvel.'/1000,Malt.'/1000,'Color',[0.6350 0.0780 0.1840]); hold on;
% xlim([0 25]); ylim([0 100]); ylabel('Altitude (km)'); xlabel('Velocity (km/s)');
% grid on; grid minor;
% 
% subplot(2,2,3)
% plot(Vang.',Valt.'/1000,'Color',[0.8500 0.3250 0.0980]); hold on; 
% plot(Eang.',Ealt.'/1000,'Color',[0 0.4470 0.7410]); hold on; 
% plot(Mang.',Malt.'/1000,'Color',[0.6350 0.0780 0.1840]); hold on;
% ylim([0 100]); ylabel('Altitude (km)'); xlabel('Angle of Trajectory (deg.)'); 
% grid on; grid minor;