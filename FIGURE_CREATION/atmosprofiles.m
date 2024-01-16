%% ATMOSPROFILES
 
% Developed by: Marissa Dudek (1)
% (1) University of North Carolina at Chapel Hill, Geological Sciences
 
% Clean workspace
clc; format long; 
set(0,'DefaultFigureWindowStyle','docked')
counter = 0;
  
% Import MATALB workspace
load('Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.mat'); % Terrestrial analog popultion (est.2024)
 
% ----------------------------------------------------------------------------

% Establish zeros matricies
proj_alt = (1:1:200000).';
proj_ztemp = zeros(200000,1);
proj_zrho = zeros(200000,1);
proj_zpress = zeros(200000,1);
flag.zrho(1) = 2;

% Call target variables
% ... g         gravity                 m/s2
% ... H         scale height            m
% ... R         radius                  m
% ... Dsc       s-c crater              m
% ... rho_t     surf. density           kg/m3
% ... diam_lim  ab/frag limit           m    
% ... zrho      atmos.density @ surf    kg/m3
% ... zpress    atmos.pressure @ surf   Pa
% ... Venus
    g.Venus = 8.87;
    H.Venus = 15900;
    R.Venus = 6052000;
    Dsc.Venus = 40000; 
    zrho.Venus = 67;
    zpress.Venus = 93.22*(10^5);
% ... Earth 
    g.Earth = 9.81;
    H.Earth = 8500;
    R.Earth = 6378000;
    Dsc.Earth = 3200;
    zrho.Earth = 10.13;
    zpress.Earth = 10.13*(10^4);
% ... Mars
    g.Mars = 3.79; 
    H.Mars = 11100;
    R.Mars = 1737500;
    Dsc.Mars = 9000;
    zrho.Mars = 0.02;
    zpress.Mars = 636;
 
% ----------------------------------------------------------------------------
 
% Venus
for j = 1:200000
% ---------------------------
% Blumenthal, Kay, Palen, Smith (2012). Understanding Our Universe. New York: W.W. Norton & Company. p. 167. ISBN 9780393912104.
% ... temperature (C)
if proj_alt(j) >= 150000
    proj_ztemp(j) = 127;
elseif proj_alt(j) >= 110000
    proj_ztemp(j) = ((proj_alt(j)-150000)/((150000-110000)/((127)-(-112))))+(127);    
elseif proj_alt(j) >= 100000
    proj_ztemp(j) = -112;
elseif proj_alt(j) >= 90000
    proj_ztemp(j) = ((proj_alt(j)-100000)/((100000-90000)/((-112)-(-104))))+(-112);
elseif proj_alt(j) >= 80000
    proj_ztemp(j) = ((proj_alt(j)-90000)/((90000-80000)/((-104)-(-76))))+(-104);
elseif proj_alt(j) >= 70000
    proj_ztemp(j) = ((proj_alt(j)-80000)/((80000-70000)/((-76)-(-43))))+(-76);
elseif proj_alt(j) >= 60000
    proj_ztemp(j) = ((proj_alt(j)-70000)/((70000-60000)/((-43)-(-10))))+(-43);
elseif proj_alt(j) >= 55000
    proj_ztemp(j) = ((proj_alt(j)-60000)/((60000-55000)/((-10)-(27))))+(-10);        
elseif proj_alt(j) >= 50000
    proj_ztemp(j) = ((proj_alt(j)-55000)/((55000-50000)/((27)-(75))))+(27);
elseif proj_alt(j) >= 45000
    proj_ztemp(j) = ((proj_alt(j)-50000)/((50000-45000)/((75)-(110))))+(75);  
elseif proj_alt(j) >= 40000
    proj_ztemp(j) = ((proj_alt(j)-45000)/((45000-40000)/((110)-(143))))+(110);
elseif proj_alt(j) >= 30000
    proj_ztemp(j) = ((proj_alt(j)-40000)/((40000-30000)/((143)-(222))))+(143);
elseif proj_alt(j) >= 20000
    proj_ztemp(j) = ((proj_alt(j)-30000)/((30000-20000)/((222)-(306))))+(222);
elseif proj_alt(j) >= 10000
    proj_ztemp(j) = ((proj_alt(j)-20000)/((20000-10000)/((306)-(385))))+(306);
elseif proj_alt(j) >= 0
    proj_ztemp(j) = ((proj_alt(j)-10000)/((10000-0)/(385-(462))))+(385);
end
% ---------------------------
% ... convert temperature from celcius to kelvin (K)
proj_ztemp(j) = (proj_ztemp(j) + 273.15);
% ---------------------------
% ... density (kg/m3)
if flag.zrho(1) == 1
    proj_zrho(j) = 67 * (exp(-proj_alt(j)/H.target));
else
    % Herrick and Phillips (pg. 260-261)
    if proj_alt(j) >= 70000
        proj_zrho(j) = 0.0789 * (exp(-23.3*((proj_alt(j)-70000)/100000)));
        %proj_zrho(j) = 67 * (exp(-proj_alt(j)/H.target));
    elseif proj_alt(j) < 70000
        proj_zrho(j) =  67 * ((1-(proj_alt(j)/94800))^5);    
    end
end
% ---------------------------
% ... pressure (Pa)
proj_zpress(j) = (93.22*(10^5)) * (exp((-g.Venus*0.04334*proj_alt(j))/(8.314*737))); 
% ---------------------------
end
Vztemp = proj_ztemp;
Vzrho = proj_zrho;
Vzpress = proj_zpress;
 
% Earth
for j = 1:200000
% ---------------------------
% ... temperature (C)
if proj_alt(j) >= 110000
    proj_ztemp(j) = 200;
elseif proj_alt(j) >= 90500
    proj_ztemp(j) = ((proj_alt(j)-110000)/((110000-90500)/(200-(-83))))+(200);
elseif proj_alt(j) >= 80000
    proj_ztemp(j) = -83;
elseif proj_alt(j) >= 60500
    proj_ztemp(j) = ((proj_alt(j)-80000)/((80000-60500)/((-83)-(-17))))+(-83);
elseif proj_alt(j) >= 52500
    proj_ztemp(j) = ((proj_alt(j)-60500)/((60500-52500)/((-17)-(-3))))+(-17);
elseif proj_alt(j) >= 48500
    proj_ztemp(j) = -3;
elseif proj_alt(j) >= 32000
    proj_ztemp(j) = ((proj_alt(j)-48500)/((48500-32000)/((-3)-(-45))))+(-3);
elseif proj_alt(j) >= 20000
    proj_ztemp(j) = ((proj_alt(j)-32000)/((32000-20000)/((-45)-(-57))))+(-45);
elseif proj_alt(j) >= 11000
    proj_ztemp(j) = -57;
elseif proj_alt(j) >= 0
    proj_ztemp(j) = ((proj_alt(j)-11000)/((11000-0)/(14+(-57))))+(-57);
end
% ---------------------------
% ... convert temperature from celcius to kelvin
proj_ztemp(j) = (proj_ztemp(j) + 273.15);
% ---------------------------
% ... density (kg/m3)
proj_zrho(j) = (1.20) * (exp(-proj_alt(j)/H.Earth)); 
% ---------------------------
% ... pressure (Pa)
proj_zpress(j) = (10.13*(10^4)) * (exp((-g.Earth*0.02896*proj_alt(j))/(8.314*288)));
% --------------------------- 
end
Eztemp = proj_ztemp;
Ezrho = proj_zrho;
Ezpress = proj_zpress;
 
% Mars
for j = 1:200000
% ---------------------------
% ... temperature (K)
if proj_alt(j) >= 150000
    proj_ztemp(j) = 300;           
elseif proj_alt(j) >= 120000
    proj_ztemp(j) = ((proj_alt(j)-150000)/((150000-120000)/(300-(175))))+(300);         
elseif proj_alt(j) >= 100000
    proj_ztemp(j) = ((proj_alt(j)-120000)/((120000-100000)/(175-(125))))+(175); 
elseif proj_alt(j) >= 70000
    proj_ztemp(j) = ((proj_alt(j)-90000)/((90000-70000)/(125-(125))))+(125);   
elseif proj_alt(j) >= 40000
    proj_ztemp(j) = ((proj_alt(j)-70000)/((70000-40000)/(125-(150))))+(125);   
elseif proj_alt(j) >= 30000
    proj_ztemp(j) = ((proj_alt(j)-40000)/((40000-30000)/(150-(170))))+(150);
elseif proj_alt(j) >= 20000
    proj_ztemp(j) = ((proj_alt(j)-30000)/((30000-20000)/(170-(190))))+(170);
elseif proj_alt(j) >= 0
    proj_ztemp(j) = ((proj_alt(j)-20000)/((20000-0)/(190-(209))))+(190);
end
% ---------------------------
% ... density (kg/m3)
proj_zrho(j) = (0.02) * (exp(-proj_alt(j)/H.Mars)); 
% ---------------------------
% ... pressure (Pa)
proj_zpress(j) = (636) * (exp((-g.Mars*0.04334*proj_alt(j))/(8.314*210)));
% --------------------------- 
end
Mztemp = proj_ztemp;
Mzrho = proj_zrho;
Mzpress = proj_zpress;

% Venus
    msx = Vztemp; msy = proj_alt/1000;
    lx = (1000/1e7)*(Vzpress); ly = proj_alt/1000;
    evx = (1000/100)*Vzrho; evy = proj_alt/1000;
    figure(1)
    a = axes('units','normalized','position',[.1 .35 .7 .6],'xtick',0:100:1000);
    plot(msx,msy,'color','k','LineWidth',2); hold on; xlim([0 1000]); 
    plot(lx,ly,'color','#0072b2','LineWidth',2); hold on
    plot(evx,evy,'color','#d55e00','LineWidth',2);
    title('Venus'); ylabel('Altitude [km]'); 
    yline(50,'--'); yline(75,'--'); yline(110,'--');
    text(610,25,'Troposphere'); text(500,62.5,'Lower Stratosphere'); text(500,92.5,'Upper Statosphere'); text(500,145,'Thermosphere');
    xlabel(a,'Temperature [K]'); hold on;
    b = axes('units','normalized','position',[.1 .25 .7 0.000001],'xlim',[0 1e7],'color','#0072b2','xtick',0:1000000:1e7,'XColor','#0072b2');
    xlabel(b,'Pressure [Pa]'); b.FontSize = 10;
    c = axes('units','normalized','position',[.1 .15 .7 0.000001],'xlim',[0 100],'color','#d55e00','xtick',0:10:100,'XColor','#d55e00');
    xlabel(c,'Density [kg/m^3]'); c.FontSize = 10;
    
% Earth
    msx = Eztemp; msy = proj_alt/1000;
    lx = (1000/5e5)*(Ezpress); ly = proj_alt/1000;
    evx = (1000/10)*Ezrho; evy = proj_alt/1000;
    figure(2)
    a = axes('units','normalized','position',[.1 .35 .7 .6],'xtick',0:50:1000);
    plot(msx,msy,'color','k','LineWidth',2); hold on; xlim([0 1000]);
    plot(lx,ly,'color','#0072b2','LineWidth',2); hold on
    plot(evx,evy,'color','#d55e00','LineWidth',2); 
    title('Earth'); ylabel('Altitude [km]'); 
    yline(12,'--'); yline(50,'--'); yline(85,'--');
    text(610,6,'Troposphere'); text(550,30,'Stratosphere'); text(550,67.5,'Mesosphere'); text(550,145,'Thermosphere');
    xlabel(a,'Temperature [K]'); hold on;
    b = axes('units','normalized','position',[.1 .25 .7 0.000001],'xlim',[0 5e5],'color','#0072b2','xtick',0:50000:5e5,'XColor','#0072b2');
    xlabel(b,'Pressure [Pa]'); b.FontSize = 10;
    c = axes('units','normalized','position',[.1 .15 .7 0.000001],'xlim',[0 10],'color','#d55e00','xtick',0:1:10,'XColor','#d55e00');
    xlabel(c,'Density [kg/m^3]'); c.FontSize = 10;
    
% Mars
    msx = Mztemp; msy = proj_alt/1000;
    lx = (1000/10000)*(Mzpress); ly = proj_alt/1000;
    evx = (1000/1)*Mzrho; evy = proj_alt/1000;
    figure(3)
    a = axes('units','normalized','position',[.1 .35 .7 .6],'xtick',0:50:1000);
    plot(msx,msy,'color','k','LineWidth',2); hold on; xlim([0 1000]);
    plot(lx,ly,'color','#0072b2','LineWidth',2); hold on
    plot(evx,evy,'color','#d55e00','LineWidth',2);
    title('Mars'); ylabel('Altitude [km]'); 
    yline(40,'--'); yline(100,'--');
    text(500,25,'Troposphere'); text(500,70,'Stratosphere'); text(500,145,'Thermosphere');
    xlabel(a,'Temperature [K]'); hold on;
    b = axes('units','normalized','position',[.1 .25 .7 0.000001],'xlim',[0 10000],'color','#0072b2','xtick',0:2000:10000,'XColor','#0072b2');
    xlabel(b,'Pressure [Pa]'); b.FontSize = 10;
    c = axes('units','normalized','position',[.1 .15 .7 0.000001],'xlim',[0 1],'color','#d55e00','xtick',0:0.1:1,'XColor','#d55e00');
    xlabel(c,'Density [kg/m^3]'); c.FontSize = 10;
 