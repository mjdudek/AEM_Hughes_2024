%% Atmospheric Effects Model (AEM)

% Developed by: Marissa Dudek (1), R Shane McGary (2)
% (1) University of North Carolina at Chapel Hill, Geological Sciences
% (2) James Madison University, Geology and Environmental Sciences

% Contact information
% ... Marissa Dudek, mjdudek@email.unc.edu

% Last udpate
% ... January 8, 2024

% ---------------------------------------------------------

% Clean workspace
clear all; clc; format long; 
counter = 0;

% Citations following this convention
% ... Author(s) (year) Brief title
% ... cite(length(cite),:) = {'Author(s)', 'Year', 'Brief title'}
cite = {'Author(s)', 'Year', 'Brief title'};

% Required MATLAB App Library
% ... Deep Learning Toolbox (MathWorks): all combinations
% ... Curve Fitting Toolbox (MathWorks): curve fitting
% ... Violin Plot (Antoine Legouhy): al_goodplot - boxblot & violin plot 
set(0,'DefaultFigureWindowStyle','docked')

%% 1A. User input variables
% In this section, the user will input the location for the source of the projectile population[s] within the asteroid belt [outer, average, inner],
% the target planetary body within the inner solar system [Mercury, Venus, Earth, Mars], the number of timesteps as the projectile travels from the
% asteroid belt to the target body, the number of timesteps as the projectile travels from the top of the target's armosphere - if applicable - to the
% target surface, as well as the timespan over which impact cratering is occuring [Gyr]. The user will also select any combination of effects to occur 
% during atmospheric travel, including ablation [change in mass, velocity, angle] and and break-up [change in the number of fragments] of the
% projectile. 

% ---------------------------------------------------------

% Import FLAVRS matrix (projectile population)
% ... load in with FLAVRS format (n x 11) as .mat 
% ... column = characterstic, row = projectile
%FLAVRS = load('FLAVRS_TerrestrialAnalogPopulation.mat'); % Terrestrial analog popultion (est.2024)
FLAVRS = load('FLAVRS_TheoreticalPopulation.mat'); % Theoretical Pop (est.2024)

% Read FLAVRS strucutre
FLAVRS = FLAVRS.FLAVRS;

% Projectile population size
% ... dim       dimensions of FLAVRS matrix 
% ... n.proj    number of projectile total
dim = size(FLAVRS);
n.projectiles = dim(1,1);   

% Define projectile densities 
% ... ch    chondritic
% ... ir    iron
% ... si    stony iron
% ... st    stony
rhop.ch(1) = 2700; %3110, 2700
rhop.ir(1) = 7800; %7600, 7800
rhop.si(1) = 4500; %4396, 4500
rhop.st(1) = 3300; %3320, 3300

% Input timesteps
% ... step_mod  number of time steps each second is broken into (e.g. 10=tenth, 100=hundreth, 1000=milisecond)
% ... z_steps   number of time steps from the target's upper atmosphere to the surface (e.g. 60=1m, 600=10m, 900=15m)
% ... z_scale   when multiplied by the target's atmospheric scale height (H), is the maximum atmospheric altitude
input.step_mod(1) = 10;
input.z_steps(1) = 70*input.step_mod(1); 
input.z_scale(1) = 1;

% Choose any or no atmospheric effects 
% ... ['0' off \ '1' on]
% ... ATMOS     atmospheric effects 
% ... DVDT      projectile velocity
% ... DADT      projectile angle of trajectory
% ... DMDT      effect of ablation on the projectile
% ... BREAK     effect of breakup on the projectile 
flag.ATMOS(1) = 1;
flag.DVDT(1) = 1;
flag.DADT(1) = 1;
flag.DMDT(1) = 1;
flag.BREAK(1) = 1;

% Flag options
% ... vel       change in velocity equation     [1958, 1971, 1980, 1995, 2003, 2004, 2005]
% ... shape     projectile shape factor         [1.21, 1.919] 
% ... drag      drag coefficent                 [0.5, 1, 2, 1971, 2016]
% ... ang       change in angle equations       [1971, 1980, 1993]
% ... lift      lift coeffecient                [0, 0.001, 0.01, 0.1]
% ... mass      change in mass equation         [1958, 1971, 1980, 1981]
% ... htc       heat transfer coefficney        [0.02]
% ... hao       heat of ablation                [1.89e6, 5e6, 8.01e6]
flag.vel(1) = 1971; 
flag.shape(1) = 1.21; 
flag.drag(1) = 0.47; 
flag.ang(1) = 1971; 
flag.lift(1) = 0.01;
flag.mass(1) = 1971;
flag.htc(1) = 0.02;
flag.hoa(1) = 5e6; 
flag.zrho(1) = 1; 

%% 1B. Planetary target variables
% In this section the target characteristics are defined
% ---------------------------------------------------------

% Planetary target loop
for m = 1:4
    
% Call target variables
% ... g         gravity                 m/s2
% ... H         scale height            m
% ... R         radius                  m
% ... Dsc       s-c crater              m
% ... rho_t     surf. density           kg/m3
% ... diam_lim  ab/frag limit           m    
% ... zrho      atmos.density @ surf    kg/m3
% ... zpress    atmos.pressure @ surf   Pa
if m == 1 % Venus
    g.target = 8.87;
    H.target = 15900;
    R.target = 6052000;
    Dsc.target = 40000; 
    rho_t = 2830; 
    diam_lim = 10000;
    zrho.target = 67;
    zpress.target = 93.22*(10^5);
elseif m == 2 % Earth 
    g.target = 9.81;
    H.target = 8500;
    R.target = 6378000;
    Dsc.target = 3200;
    rho_t = 2830;
    diam_lim = 1000;
    zrho.target = 10.13;
    zpress.target = 10.13*(10^4);
elseif m == 3 % Mars
    g.target = 3.79; 
    H.target = 11100;
    R.target = 1737500;
    Dsc.target = 9000;
    rho_t = 2830;
    diam_lim = 100;
    zrho.target = 0.02;
    zpress.target = 636;
elseif m == 4 % Moon 
    g.target = 1.62; 
    R.target = 1737500;
    Dsc.target = 15000; % m, Pike (1980)
    rho_t = 2830; 
end

%% 4.1 Initial conditions

for i = 1

% Define global values and matrices to be populated, in alphabetical order
% ... atmos.max m       maximum altitude
% ... alt       m       altitude (z)
% ... ang       deg.    angle of trajectory relative to horizontal
% ... cfa       m2      circualr frontal area
% ... dist      m       distance (x)
% ... drag      -       drag coefficient
% ... frag      #       number of projectile fragments
% ... htc       -       heat tranfer coeffient
% ... hoa       J/kg    heat of ablation 
% ... lift      -       lift coefficient
% ... Ld        m       diameter of projectile
% ... Lz        m       diameter of frontal area
% ... mass      kg      projectile mass
% ... massr     -       ratio of projectile mass to initial mass
% ... panfact   -       pancake factor
% ... Ps        Pa      stagnation pressure
% ... rad       m       radius of projectile fragments
% ... reyn      ?       reynold's number
% ... time      s       travel time
% ... vel       m/s     velocity
% ... velb      m/s     velocity at breakup
% ... velt      m/s     terminal velocity
% ... xarea     m2      effective cross sectional area
% ... Yi        Pa      yield strength
% ... z_break   m       altitude at breakup    
% ... z_burst   m       altitude of airburst
% ... zrho      kg/m3   atmospheric density
% ... ztemp     K       atmospheric temperature
% ... zpress    Pa      atmospheric pressure

% Set maximum altitude 
%atmos.max = (input.z_scale*H.target*25); 
atmos.max = 200000; 

% Set global values and matrices to be populated, in alphabetical order
% ... variables changing with altitude
proj_alt = NaN(n.projectiles, input.z_steps);
proj_ang = zeros(n.projectiles, input.z_steps);
proj_dist = zeros(n.projectiles, input.z_steps);
proj_drag = zeros(n.projectiles, input.z_steps);
proj_frag = zeros(n.projectiles, input.z_steps);
proj_htc = zeros(n.projectiles, input.z_steps);
proj_hoa = zeros(n.projectiles, input.z_steps);
proj_lift = zeros(n.projectiles, input.z_steps);
proj_mass = NaN(n.projectiles, input.z_steps);
proj_massr = zeros(n.projectiles, input.z_steps);
proj_cfa = zeros(n.projectiles, input.z_steps);
proj_diam = zeros(n.projectiles, input.z_steps);
proj_reyn = zeros(n.projectiles, input.z_steps);
proj_shape = zeros(n.projectiles, input.z_steps);
proj_time = zeros(n.projectiles, input.z_steps);
proj_vel = zeros(n.projectiles, input.z_steps);
proj_velcr = zeros(n.projectiles, input.z_steps);
proj_velt = zeros(n.projectiles, input.z_steps);
proj_xarea = zeros(n.projectiles, input.z_steps);
proj_zrho = zeros(n.projectiles, input.z_steps);
proj_ztemp = zeros(n.projectiles, input.z_steps);
proj_zpress = zeros(n.projectiles, input.z_steps);

% Set breakup variables
Ps = zeros(n.projectiles, input.z_steps);
Yi = zeros(n.projectiles, input.z_steps);
If = zeros(n.projectiles, input.z_steps);
Ld = zeros(n.projectiles, input.z_steps);
proj_velb = zeros(n.projectiles, input.z_steps);
z_break = zeros(n.projectiles,input.z_steps);
Lz = zeros(n.projectiles, input.z_steps);
panfact = ones(n.projectiles, input.z_steps);
z_burst = NaN(n.projectiles, input.z_steps);

% Calcualte initial conditions
for k = 1:n.projectiles % outer loop
    % populated in order of appearance 
    % ---------------------------------------------------------
    % atmosphere
    proj_alt(k,1) = atmos.max; 
    % temperature at altitude (K)
    % pressure at altitude (Pa)
    % density at altitude (kg/m3)
    % ... see next section for citations and equation validation
    if m == 1 % Venus
        proj_ztemp(k,1) = 108 + 273.15; 
        proj_zpress(k,1) = zpress.target * (exp((-g.target*0.04334*proj_alt(k,1))/(8.314*737)));  
        proj_zrho(k,1) = zrho.target * (exp(-proj_alt(k,1)/H.target));
    elseif m == 2 % Earth
        proj_ztemp(k,1) = 200 + 273.15; 
        proj_zpress(k,1) = zpress.target * (exp((-g.target*0.02896*proj_alt(k,1))/(8.314*288)));
        proj_zrho(k,1) = zrho.target  * (exp(-proj_alt(k,1)/H.target));
    elseif m == 3 % Mars
        proj_ztemp(k,1) = 300; 
        proj_zpress(k,1) = zpress.target * (exp((-g.target*0.04334*proj_alt(k,1))/(8.314*210)));
        proj_zrho(k,1) = zrho.target * (exp(-proj_alt(k,1)/H.target));
    end
    % ---------------------------------------------------------
    % projectile characteristics
    % ... angle of trajectroy (deg. from horizontal)
    proj_ang(k,1) = FLAVRS(k,11);  
    % ... distance (x)
    proj_dist(k,1) = 0;   
    % ... number of initial fragments (1)
    proj_frag(k,1) = 1;    
    % ... mass (kg)
    proj_mass(k,1) = FLAVRS(k,7);
    % ... mass ratio (%)
    proj_massr(k,1) = (proj_mass(k,1)/proj_mass(k,1))*100;
    % ---------------------------------------------------------
    % projectile area
    proj_xarea(k,1) = flag.shape(1)*((proj_mass(k,1)/FLAVRS(k,5))^(2/3));  
    % ---------------------------------------------------------
    % projectile diameter
    proj_diam(k,1) = (nthroot((FLAVRS(k,7)/((4/3)*(pi)*FLAVRS(k,5))),3))*2;
    % ---------------------------------------------------------
    % circular frontal area (*corrected from radius 220425)
    if k == 1; cite(length(cite)+1,:) = {'Baldwin and Shaeffer', '1971', 'Crater Field Formation (Appedix.2)'}; end    
    proj_cfa(k,1) = (pi^(-1/2))*(flag.shape(1)^(1/2))*(FLAVRS(k,5)^(-1/3))*(proj_mass(k,1)^(1/3))*(proj_frag(k,1)^(-1/3));
    % ---------------------------------------------------------
    % projectile velocity
    proj_vel(k,1) = FLAVRS(k,4);
    % ---------------------------------------------------------
    % Renold's number
    proj_reyn(k,1) = ( (proj_zrho(k,1)*proj_vel(k,1)*proj_diam(k,1)) / (2.791e-7*(proj_ztemp(k,1)^(0.7355))) );
    % ---------------------------------------------------------
    % drag
    if flag.drag(1) == 1971
        
         % ... calculated using Baldwin and Shaeffer (1971)
        if k == 1; cite(length(cite)+1,:) = {'Baldwin and Shaeffer', '1971', 'Ablation and Breakup of Meteoroids (Appedix.1)'}; end
        proj_drag(k,1) = 1 + exp(-1.837E6*(proj_zrho(k,1)/zrho.target(1))*proj_cfa(k,1));   
    
        elseif flag.drag(1) == 2011
        % ... calculated using Polezhaev and Chirov (2011)
        if k == 1; cite(length(cite)+1,:) = {'Polezhaev and Chirov', '2011', 'Thermopedia, Drag Coefficient [https://www.thermopedia.com/content/707/]'}; end
        if proj_reyn(k,1) <= 0.02
            proj_drag(k,1) = 24/proj_reyn(k,1);
            
        else
            proj_drag(k,1) = ((21.12/proj_reyn(k,1))+(6.3/sqrt(proj_reyn(k,1)))+0.25);
        end
        
    elseif flag.drag(1) == 2016
        % ... calculated using Morrison
        if k == 1; cite(length(cite)+1,:) = {'Morrison', '2013', 'Data Correlation for Drag Coefficient for Sphere '}; end
        if proj_reyn(k,1) <= 2
            proj_drag(k,1) = (24/proj_reyn(k,1));
        else
            proj_drag(k,1) = (24/proj_reyn(k,j) + ...
                        ( (2.6*(proj_reyn(k,1)/5)) / (1+(proj_reyn(k,1)/5)^1.52) ) + ...
                        ( (0.411*(proj_reyn(k,1)/2.63e5)^7.94) / (1+(proj_reyn(k,1)/2.63e5)^8.00) ) + ...
                        ( (0.25*(proj_reyn(k,1)/10e6)) / (1+(proj_reyn(k,1)/10e6)) ));
        end
        
    else
        proj_drag(k,1) = flag.drag(1);
        
    end  
    % ---------------------------------------------------------
    % lift   
    proj_lift(k,1) = flag.lift(1); 
    % ---------------------------------------------------------
    % heat transfer coefficent    
    proj_htc(k,1) = flag.htc(1);   
    % ---------------------------------------------------------
    % heat of ablation   
    proj_hoa(k,1) = flag.hoa(1);
    % ---------------------------------------------------------
    % terminal velocity
    proj_velt(k,1) = (sqrt((2*proj_mass(k,1)*g.target) / (proj_zrho(k,1)*proj_xarea(k,1)*proj_drag(k,1)) ));
    % ---------------------------------------------------------
    % critical velocity
    if k == 1; cite(length(cite)+1,:) = {'Passey and Melosh', '1980', 'Crater Field Formation (Eq.8)'}; end
    proj_velcr(k,1) = 3000;
    % ---------------------------------------------------------
    % time
    proj_time(k,1) = 0;
    % ---------------------------------------------------------
    % shape factor
    proj_shape(k,1) = flag.shape(1);
    % ---------------------------------------------------------
    % pancake factor
    panfact(k,1) = 1;
    % ---------------------------------------------------------
    
end % k = 1:n.projectiles % outer loop

end % Initial conditoins

%% 4.2 Atmospheric effects

% Calculate atmospheric effects
% =========================================================
% Start looping

% Outer loop to run through every projectile (row)
count_out = 0;
for k = 1:n.projectiles
   
% Inner loop to run through every time step (column)
count_in = 0;
for j = 2:input.z_steps 
    
if proj_alt(k,j-1) <= 0
    proj_mass(k,j) = proj_mass(k,j-1);
    proj_massr(k,j) = proj_massr(k,j-1);
elseif isnan(proj_alt(k,j-1))
    proj_alt(k,j-1) = 0;
    proj_alt(k,j) = 0;
end % if proj_vel(k,j-1) <= 0

% if proj_velt(k,j-1) <= proj_vel(k,j-1)
%     proj_vel(k,j) = g.target;
% end % if proj_vel(k,j-1) <= 0
    
% =========================================================
% Change in altitude

for i = 1

    if k == 1; cite(length(cite)+1,:) = {'Passey and Melosh', '1980', 'Crater Field Formation (Eq.4)'}; end       
    % ... previous altitude minus (the velocity at the previous step divided by the time step mod) multiplied by the sin(previous theta)
    %proj_alt(k,j) = proj_alt(k,j-1) - (( (proj_vel(k,j-1)) * (sind(proj_ang(k,j-1))) )/input.step_mod(1)); 
    proj_alt(k,j) = proj_alt(k,j-1) - ( (proj_vel(k,j-1) * (sind(proj_ang(k,j-1))) )/input.step_mod(1)); 
    
    % ... if the altitude of the step is 0, the altitude remains zero and does not go negative 
    if proj_alt(k,j) < 0 
        proj_alt(k,j) = 0;
    end
    
end

% ---------------------------------------------------------
% Change in distance

for i = 1 
    
    if k == 1; cite(length(cite)+1,:) = {'Passey and Melosh', '1980', 'Crater Field Formation (Eq.9)'}; end
    % ... previous distance plus (the previous velocity multiplied by sin(previous theta)) divided by ( 1 plus (height / radius of the target))
    proj_dist(k,j) = proj_dist(k,j-1) + (( ((proj_vel(k,j-1))*cosd(proj_ang(k,j-1))) / (1+(proj_alt(k,j)/R.target)) )/input.step_mod(1));
    
    % ... if the altitude of the step is 0, the distance traveled is unchanged
    if proj_alt(k,j-1) <=0 || proj_alt(k,j)<=0
        proj_dist(k,j) = proj_dist(k,j-1);
    end    

end

% ---------------------------------------------------------
% Change in atmospheric temperature, pressure, and density

for i = 1
    
% Venus atmospheric temperature profiles
if m == 1

    % ---------------------------
    % Blumenthal, Kay, Palen, Smith (2012). Understanding Our Universe. New York: W.W. Norton & Company. p. 167. ISBN 9780393912104.
    % ... temperature (C)
    if proj_alt(k,j) >= 150000
        proj_ztemp(k,j) = 127;
    elseif proj_alt(k,j) >= 110000
        proj_ztemp(k,j) = ((proj_alt(k,j)-150000)/((150000-110000)/((127)-(-112))))+(127);    
    elseif proj_alt(k,j) >= 100000
        proj_ztemp(k,j) = -112;
    elseif proj_alt(k,j) >= 90000
        proj_ztemp(k,j) = ((proj_alt(k,j)-100000)/((100000-90000)/((-112)-(-104))))+(-112);
    elseif proj_alt(k,j) >= 80000
        proj_ztemp(k,j) = ((proj_alt(k,j)-90000)/((90000-80000)/((-104)-(-76))))+(-104);
    elseif proj_alt(k,j) >= 70000
        proj_ztemp(k,j) = ((proj_alt(k,j)-80000)/((80000-70000)/((-76)-(-43))))+(-76);
    elseif proj_alt(k,j) >= 60000
        proj_ztemp(k,j) = ((proj_alt(k,j)-70000)/((70000-60000)/((-43)-(-10))))+(-43);
    elseif proj_alt(k,j) >= 55000
        proj_ztemp(k,j) = ((proj_alt(k,j)-60000)/((60000-55000)/((-10)-(27))))+(-10);        
    elseif proj_alt(k,j) >= 50000
        proj_ztemp(k,j) = ((proj_alt(k,j)-55000)/((55000-50000)/((27)-(75))))+(27);
    elseif proj_alt(k,j) >= 45000
        proj_ztemp(k,j) = ((proj_alt(k,j)-50000)/((50000-45000)/((75)-(110))))+(75);  
    elseif proj_alt(k,j) >= 40000
        proj_ztemp(k,j) = ((proj_alt(k,j)-45000)/((45000-40000)/((110)-(143))))+(110);
    elseif proj_alt(k,j) >= 30000
        proj_ztemp(k,j) = ((proj_alt(k,j)-40000)/((40000-30000)/((143)-(222))))+(143);
    elseif proj_alt(k,j) >= 20000
        proj_ztemp(k,j) = ((proj_alt(k,j)-30000)/((30000-20000)/((222)-(306))))+(222);
    elseif proj_alt(k,j) >= 10000
        proj_ztemp(k,j) = ((proj_alt(k,j)-20000)/((20000-10000)/((306)-(385))))+(306);
    elseif proj_alt(k,j) >= 0
        proj_ztemp(k,j) = ((proj_alt(k,j)-10000)/((10000-0)/(385-(462))))+(385);
    end
    % ---------------------------
    % ... convert temperature from celcius to kelvin (K)
    proj_ztemp(k,j) = (proj_ztemp(k,j) + 273.15);
    % ---------------------------
    % ... density (kg/m3)
    if flag.zrho(1) == 1
        proj_zrho(k,j) = zrho.target * (exp(-proj_alt(k,j)/H.target));
    else
        % Herrick and Phillips (pg. 260-261)
        if proj_alt(k,j) >= 70000
            proj_zrho(k,j) = 0.0789 * (exp(-23.3*((proj_alt(k,j)-70000)/100000)));
        else 
            proj_zrho(k,j) = zrho.target * ((1-(proj_alt(k,j)/94800))^5);    
        end
    end
    % ---------------------------
    % ... pressure (Pa)
    proj_zpress(k,j) = zpress.target * (exp((-g.target*0.04334*proj_alt(k,j))/(8.314*737))); 
    % ---------------------------
end

% Earth atmospheric profiles
if m == 2
    % ---------------------------
    % ... temperature (C)
    if proj_alt(k,j) >= 110000
        proj_ztemp(k,j) = 200;
    elseif proj_alt(k,j) >= 90500
        proj_ztemp(k,j) = ((proj_alt(k,j)-110000)/((110000-90500)/(200-(-83))))+(200);
    elseif proj_alt(k,j) >= 80000
        proj_ztemp(k,j) = -83;
    elseif proj_alt(k,j) >= 60500
        proj_ztemp(k,j) = ((proj_alt(k,j)-80000)/((80000-60500)/((-83)-(-17))))+(-83);
    elseif proj_alt(k,j) >= 52500
        proj_ztemp(k,j) = ((proj_alt(k,j)-60500)/((60500-52500)/((-17)-(-3))))+(-17);
    elseif proj_alt(k,j) >= 48500
        proj_ztemp(k,j) = -3;
    elseif proj_alt(k,j) >= 32000
        proj_ztemp(k,j) = ((proj_alt(k,j)-48500)/((48500-32000)/((-3)-(-45))))+(-3);
    elseif proj_alt(k,j) >= 20000
        proj_ztemp(k,j) = ((proj_alt(k,j)-32000)/((32000-20000)/((-45)-(-57))))+(-45);
    elseif proj_alt(k,j) >= 11000
        proj_ztemp(k,j) = -57;
    elseif proj_alt(k,j) >= 0
        proj_ztemp(k,j) = ((proj_alt(k,j)-11000)/((11000-0)/(14+(-57))))+(-57);
    end
    % ---------------------------
    % ... convert temperature from celcius to kelvin
    proj_ztemp(k,j) = (proj_ztemp(k,j) + 273.15);
    % ---------------------------
    % ... density (kg/m3)
    proj_zrho(k,j) = zrho.target * (exp(-proj_alt(k,j)/H.target)); 
    % ---------------------------
    % ... pressure (Pa)
    proj_zpress(k,j) = zpress.target * (exp((-g.target*0.02896*proj_alt(k,j))/(8.314*288)));
    % --------------------------- 
end

% Mars atmospheric profiles
if m == 3
    % ---------------------------
    % ... temperature (K)
    if proj_alt(k,j) >= 150000
        proj_ztemp(k,j) = 300;           
    elseif proj_alt(k,j) >= 120000
        proj_ztemp(k,j) = ((proj_alt(k,j)-150000)/((150000-120000)/(300-(175))))+(300);         
    elseif proj_alt(k,j) >= 100000
        proj_ztemp(k,j) = ((proj_alt(k,j)-120000)/((120000-100000)/(175-(125))))+(175); 
    elseif proj_alt(k,j) >= 70000
        proj_ztemp(k,j) = ((proj_alt(k,j)-90000)/((90000-70000)/(125-(125))))+(125);   
    elseif proj_alt(k,j) >= 40000
        proj_ztemp(k,j) = ((proj_alt(k,j)-70000)/((70000-40000)/(125-(150))))+(125);   
    elseif proj_alt(k,j) >= 30000
        proj_ztemp(k,j) = ((proj_alt(k,j)-40000)/((40000-30000)/(150-(170))))+(150);
    elseif proj_alt(k,j) >= 20000
        proj_ztemp(k,j) = ((proj_alt(k,j)-30000)/((30000-20000)/(170-(190))))+(170);
    elseif proj_alt(k,j) >= 0
        proj_ztemp(k,j) = ((proj_alt(k,j)-20000)/((20000-0)/(190-(209))))+(190);
    end
    % ---------------------------
    % ... density (kg/m3)
    proj_zrho(k,j) = zrho.target * (exp(-proj_alt(k,j)/H.target)); 
    % ---------------------------
    % ... pressure (Pa)
    proj_zpress(k,j) = zpress.target * (exp((-g.target*0.04334*proj_alt(k,j))/(8.314*210)));
    % --------------------------- 
end  

end

% ---------------------------------------------------------
% Reynold's number

for i = 1
    
    % Calculate using standard Reynolds number calculation w/ dynamic viscosity of sphere through air 
    % ... proj_reyn(k,j) = [(density of fluid)*(flow speed)*(characteristic linear dimension)] / (dynamic viscosity of fluid)
    % ... dynamic viscosity of air, [10e-6 - 10e?5] Lower/Upper range of gaseous viscosity
    %proj_reyn(k,j) = ( (proj_zrho(k,j)*(proj_vel(k,j-1)/input.step_mod(1))*proj_rad(k,j-1)) / (2.791e-7*(proj_ztemp(k,j)^(0.7355))) );
    proj_reyn(k,j) = ( (proj_zrho(k,j)*proj_vel(k,j-1)*proj_diam(k,j-1)) / (2.791e-7*(proj_ztemp(k,j)^(0.7355))) );
    
end
% ---------------------------------------------------------
% Drag coefficient

for i = 1
    
    if flag.drag(1) == 1971
         % ... calculated using Baldwin and Shaeffer (1971)
        if k == 1; cite(length(cite)+1,:) = {'Baldwin and Shaeffer', '1971', 'Ablation and Breakup of Meteoroids (Appedix.1)'}; end
        proj_drag(k,j) = 1 + exp(-1.837E6*(proj_zrho(k,j)/zrho.target(1))*proj_cfa(k,j-1));   
    
    elseif flag.drag(1) == 2011
        % ... calculated using Polezhaev and Chirov (2011)
        if k == 1; cite(length(cite)+1,:) = {'Polezhaev and Chirov', '2011', 'Thermopedia, Drag Coefficient [https://www.thermopedia.com/content/707/]'}; end
        if proj_reyn(k,j) <= 0.02
            proj_drag(k,j) = 24/proj_reyn(k,j);
        else
            proj_drag(k,j) = ((21.12/proj_reyn(k,j))+(6.3/sqrt(proj_reyn(k,j)))+0.25);
        end
        
    elseif flag.drag(1) == 2016
        % ... calculated using Morrison
        if k == 1; cite(length(cite)+1,:) = {'Morrison', '2013', 'Data Correlation for Drag Coefficient for Sphere '}; end
        if proj_reyn(k,j) <= 2
            proj_drag(k,j) = (24/proj_reyn(k,j));
        else
            proj_drag(k,j) = (24/proj_reyn(k,j) + ...
                        ( (2.6*(proj_reyn(k,j)/5)) / (1+(proj_reyn(k,j)/5)^1.52) ) + ...
                        ( (0.411*(proj_reyn(k,j)/2.63e5)^7.94) / (1+(proj_reyn(k,j)/2.63e5)^8.00) ) + ...
                        ( (0.25*(proj_reyn(k,j)/10e6)) / (1+(proj_reyn(k,j)/10e6)) ));
        end
        
    else
        proj_drag(k,j) = flag.drag(1);
    end 
    
end

% ---------------------------------------------------------
% Shape factor

for i = 1
    
if flag.shape(1) == 1995
    if k == 1; cite(length(cite)+1,:) = {'Svetsov, Nemtchinov, and Teterev', '1995', 'Disintegration of Large Meteoroids in Earth Atmos'}; end
    proj_shape(k,j) = (pi*((proj_diam(k,j-1)/2)^2)*(proj_mass(k,j-1)^(-2/3))*(FLAVRS(k,5)^(2/3)));
else
    proj_shape(k,j) = flag.shape(1);
end

end
    
% ---------------------------------------------------------
% Change in velocity

if flag.DVDT == 1
    
% Change in velocity

    if flag.vel(1) == 1958
        % ... calculate according to Opik (1958)
        if k == 1; cite(length(cite)+1,:) = {'Opik', '1958', 'Meteor Physics'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - ((proj_drag(k,j)*proj_zrho(k,j)*flag.shape(1)*((proj_vel(k,j-1))^2) ... 
                                            / ( ((proj_mass(k,j-1))^(1/3))*(FLAVRS(k,5)^(2/3)) ) )/input.step_mod(1));
                                  
    elseif flag.vel(1) == 1971
        % ... calculate according to Baldwin and Shaeffer (1971)
        if k == 1; cite(length(cite)+1,:) = {'Badlwin and Shaeffer', '1971', 'Ablation and Breakup of Large Meteoroids'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - (((1/2)*proj_drag(k,j)*proj_zrho(k,j)*((proj_vel(k,j-1))^2)* ...
                                          (proj_xarea(k,j-1)/proj_mass(k,j-1))+(g.target*cosd(90-proj_ang(k,j-1))))/input.step_mod(1));
                                      
    elseif flag.vel(1) == 1980
        % ... calculate according to Passey and Melosh (1980)
        if k == 1; cite(length(cite)+1,:) = {'Passey and Melosh', '1980', 'Crater Field Formation (Eq.1)'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - (((( (proj_drag(k,j)*proj_zrho(k,j)*proj_xarea(k,j-1)*((proj_vel(k,j-1))^2))) / proj_mass(k,j-1)) ...
                                          + (g.target*sind(proj_ang(k,j-1))))/input.step_mod(1));
                                      
    elseif flag.vel(1) == 1993
        % ... calculate according to Chyba, Thomas, and Zahnle (1993)
        if k == 1; cite(length(cite)+1,:) = {'Chyba, Thomas, and Zahnle', '1993', '[TBD]'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - (( (((1/2)*proj_drag(k,j)*proj_zrho(k,j)*proj_xarea(k,j-1)*(proj_vel(k,j-1)^2)) + ...
                                           ((g.target(1)/proj_mass(k,j-1))*(sind(proj_ang(k,j-1)))))) / input.step_mod(1)) / proj_mass(k,j-1);
                                      
    elseif flag.vel(1) == 1995
        % ... calculate according to Svetsov, Nemtchinov, and Teterev (1995)
        if k == 1; cite(length(cite)+1,:) = {'Svetsov, Nemtchinov, and Teterev', '1995', 'Disintegration of Large Meteoroids in Earth Atmos'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - ((((1/2)*proj_drag(k,j)*proj_zrho(k,j)*(((proj_vel(k,j-1))^2)*(proj_xarea(k,j-1))) ...
                                          /(proj_mass(k,j-1))))/input.step_mod(1));
                                      
     elseif flag.vel(1) == 2003
        % ... calculate according to Korcansky and Zahnle (2003)
        if k == 1; cite(length(cite)+1,:) = {'Korcansky & Zahnle', '2003', '[TBD]'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - ((proj_drag(k,j)*proj_zrho(k,j)*(((proj_vel(k,j-1))^2)*(proj_xarea(k,j-1))) ...
                                          /(proj_mass(k,j-1)))/input.step_mod(1));
                                    
     elseif flag.vel(1) == 2004
        % ... calculate according to Campbell-Brown & Koschny (2004)
        if k == 1; cite(length(cite)+1,:) = {'Campbell-Brown & Koschny', '2004', '[TBD]'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - ((((proj_drag(k,j)*proj_zrho(k,j)*((proj_vel(k,j-1))^2))/(proj_mass(k,j-1))) ...
                                            * (proj_xarea(k,j-1))*( ((proj_mass(k,j-1))/(FLAVRS(k,5)))^(2/3) )) / input.step_mod(1));
                                                                               
    elseif flag.vel(1) == 2005
        % ... calculate according to Collins, Melosh, and Marcus (2005), Equ.6
        if k == 1; cite(length(cite)+1,:) = {'Collins, Melosh, and Marcus', '2005', 'Earth Impact Effects Program'}; end
        proj_vel(k,j) = proj_vel(k,j-1) - (( ((3*proj_zrho(k,j)*proj_drag(k,j)) ...
                                            / (4*FLAVRS(k,5)*(proj_diam(k,j-1)))) ...
                                            * ((proj_vel(k,j-1))^2) )/input.step_mod(1)) ;

    end
    
    % Correct velocity equation to not be negative
    if proj_vel(k,j) <= 0
        proj_vel(k,j) = 0;
    elseif proj_vel(k,j) > proj_vel(k,j-1)
        proj_vel(k,j) = proj_vel(k,j-1);
    end
    
else

% if flag.DVDT == 0
    proj_vel(k,j) = proj_vel(k,j-1) + 0;
    
end

% ---------------------------------------------------------
% Terminal velocity

for i = 1
    
    % ... calculate terminal velocity
    if k == 1; cite(length(cite)+1,:) = {'NASA Glenn Research Center', 'Retreived 2021', 'Termlinal Velocity (gravity and drag)'}; end
    proj_velt(k,j) = (sqrt((2*proj_mass(k,j-1)*g.target)/(proj_zrho(k,j)*proj_xarea(k,j-1)*proj_drag(k,j))));
    
end

% ---------------------------------------------------------
% Lift coefficient

for i = 1

    proj_lift(k,j) = flag.lift(1);
    
end

% ---------------------------------------------------------
% Change in angle

if flag.DADT == 1
    
    if flag.ang(1) == 1971
        % ... calcualte according to Badlwin & Shaeffer, equation converted from relative to vertical to horiztonal
        if k == 1; cite(length(cite)+1,:) = {'Badlwin and Shaeffer', '1971', 'Ablation and Breakup of Large Meteoroids'}; end
        proj_ang(k,j) = proj_ang(k,j-1) + abs(-sind(( (g.target/(proj_vel(k,j))) * (1-(cosd(90-proj_ang(k,j-1)))^2))) / input.step_mod(1));

    elseif flag.ang(1) == 1980
    if k == 1; cite(length(cite)+1,:) = {'Passey and Melosh', '1980', 'Crater Field Formation (Eq.8)'}; end
    a = (( (( (proj_mass(k,j-1)*g.target*cosd(proj_ang(k,j-1))) ...
                -((1/2)*(proj_lift(k,j)*proj_zrho(k,j)*proj_xarea(k,j-1)*((proj_vel(k,j))^2))) ) ...
                / (proj_mass(k,j-1)*(proj_vel(k,j)))) ...
                - ( ((proj_vel(k,j))*cosd(proj_ang(k,j-1))) / (R.target + proj_alt(k,j))))...
                / input.step_mod(1));
    proj_ang(k,j) = proj_ang(k,j-1) + abs(a);  
    
    elseif flag.ang(1) == 1993
    % ... calcualte according to Chyba, Thomas, and Zahnle
    if k == 1; cite(length(cite)+1,:) = {'Chyba, Thomas, and Zahnle', '1993', '[TBD]'}; end
    a = ( ((g.target*cosd(proj_ang(k,j-1)))/proj_vel(k,j)) ...
                - ( (proj_lift(k,j)*proj_zrho(k,j)*proj_xarea(k,j-1)*(proj_vel(k,j)) ) / (2*proj_mass(k,j-1))) ...
                - ( ((proj_vel(k,j))*cosd(proj_ang(k,j-1))) / (R.target + proj_alt(k,j)))) ...
                / input.step_mod(1);
    proj_ang(k,j) = proj_ang(k,j-1) + abs(a);  
    
    end
    
    % Correct for overestimating angle and decreasing angle
    if proj_ang(k,j) >= 90
        proj_ang(k,j) = 90;
    elseif proj_ang(k,j) < proj_ang(k,j-1)
        proj_ang(k,j) = proj_ang(k,j-1);
    end
    
else

% if flag.DADT == 0
proj_ang(k,j) = proj_ang(k,j-1) + 0;
    
end

% ---------------------------------------------------------
% Change in mass

for i = 1

if flag.DMDT == 1
    
% For Collins, Melosh, and Marcus
if proj_diam(k,j) <= diam_lim(1)

    % Heat transfer coefficent
        proj_htc(k,j) = flag.htc(1);

    % Heat of ablation 
        proj_hoa(k,j) = flag.hoa(1);

    % Critical velocity
        proj_velcr(k,j) = proj_velcr(k,j-1);

    % Change in mass
    if flag.mass(1) == 1958
        % Calculate according to Opik
        proj_mass(k,j) = proj_mass(k,j-1) - (((proj_htc(k,j)*proj_zrho(k,j)*flag.shape(1)*((proj_vel(k,j))^3)*(proj_mass(k,j-1)^(2/3)))) ...
                                            / (2*proj_hoa(k,j)*(FLAVRS(k,5)^(2/3))) )/input.step_mod(1);  
                                        
    elseif flag.mass(1) == 1971
        % Calculate acoording to Baldwin & Shaeffer
        proj_mass(k,j) = proj_mass(k,j-1) - ((1/2)*((proj_htc(k,j)*proj_zrho(k,j) ... 
                                            *((proj_vel(k,j))^3))) ...
                                            *(proj_xarea(k,j-1)/proj_hoa(k,j))) ...
                                            /input.step_mod(1); 

    elseif flag.mass(1) == 1980
        if k == 1; cite(length(cite)+1,:) = {'Passey and Melosh', '1980', 'Crater Field Formation (Eq.8)'}; end
            proj_mass(k,j) = proj_mass(k,j-1) - (((1/2)*((proj_htc(k,j)*proj_zrho(k,j)*proj_xarea(k,j-1)*((proj_vel(k,j-1))^3)) ...
                                            / (proj_hoa(k,j)) ))/input.step_mod(1));
%         proj_mass(k,j) = proj_mass(k,j-1) - (((1/2)*((proj_htc(k,j)*proj_zrho(k,j)*proj_xarea(k,j-1)*((proj_vel(k,j-1))^2)) ...
%                                             / (proj_hoa(k,j)) ) * ( ((proj_vel(k,j)^2) - (proj_velcr(k,j)^2)) ...
%                                             / (proj_vel(k,j)^2) ) )/input.step_mod(1));

    elseif flag.mass(1) == 1981
        if k == 1; cite(length(cite)+1,:) = {'Bronshten et al.', '1981', '[TBD]'}; end
            proj_mass(k,j) = proj_mass(k,j-1) - (((proj_htc(k,j))*((proj_xarea(k,j-1)*(proj_vel(k,j-1)^3))/(2*proj_hoa(k,j)))) ...
                                                /input.step_mod(1));
                                            
    elseif flag.mass(1) == 1993
        if k == 1; cite(length(cite)+1,:) = {'Chyba, Thomas & Zahnle', '1993', '[TBD]'}; end
            proj_mass(k,j) = proj_mass(k,j-1) - (((1/2)*(proj_htc(k,j)*proj_xarea(k,j-1)*(proj_vel(k,j-1)^3)))/(proj_hoa(k,j))) ...
                                                /input.step_mod(1);

    elseif flag.mass(1) == 2004
        if k == 1; cite(length(cite)+1,:) = {'Campbell-Brown & Koschny ', '2004', '[TBD]'}; end
            proj_mass(k,j) = proj_mass(k,j-1) - ((((proj_htc(k,j)*proj_xarea(k,j-1))/(2*proj_diam(k,j))) ...
                                                *((proj_mass(k,j-1)/FLAVRS(k,5))^(2/3)) ...
                                                *(proj_zrho(k,j)*(proj_vel(k,j)^3)))/input.step_mod(1));

    else 
        proj_mass(k,j) = proj_mass(k,j-1); 
    end
    % Correct for imaginary numbers
    proj_mass(k,j) = real(proj_mass(k,j)); 
    
end % diam_lin
    
else

% if flag.DMDT == 0
    proj_mass(k,j) = proj_mass(k,j-1); 
    
end 

% Ratio of mass to initial mass
proj_massr(k,j) = (proj_mass(k,j)/proj_mass(k,1))*100;

end

% ---------------------------------------------------------  
% Change in diameter & circular frontal area

for i = 1
    
% Projectile fragments CFA
    if k == 1; cite(length(cite)+1,:) = {'Baldwin and Shaeffer', '1971', 'Ablation and Breakup (Appedix.1)'}; end
    proj_cfa(k,j) = (pi^(-1/2))*(flag.shape(1)^(1/2))*(FLAVRS(k,5)^(-1/3))*(proj_mass(k,j)^(1/3))*(proj_frag(k,j-1)^(-1/3));
    
% Projectile fragment diameter    
    proj_diam(k,j) = (nthroot((proj_mass(k,j)/((4/3)*(pi)*FLAVRS(k,5))),3))*2;
    if proj_diam(k,j) <= 0
        proj_diam(k,j) = 0;
    end

end

% ---------------------------------------------------------
% Break up

if flag.BREAK == 1
    
% For Collins, Melosh, and Marcus
if proj_diam(k,j) <= diam_lim(1)
    
    % Set initial frags
    proj_frag(k,j) = proj_frag(k,j-1);
    
    % Yeild strength & fragmentation
    % ... stagnation pressure
    % ... yeild strength of impactor (equatino works for 1000-8000 kg/m3)
    Ps(k,j) = (proj_zrho(k,j)*(proj_vel(k,j)^2)); 
    Yi(k,j) = 10^(2.107 + ( 0.0624*(sqrt(FLAVRS(k,5))) ) );

    % Breakup if stagnation pressure exceeds yield strength
    if Ps(k,j) > Yi(k,j)

        % Altitude of breakup
        z_break(k,j) = proj_alt(k,j-1);
        proj_frag(k,j) = (proj_frag(k,j)+1);
        Ld(k,j) = (proj_diam(k,j)*(sind(proj_ang(k,j)))*sqrt((FLAVRS(k,5))/(proj_drag(k,j)*(proj_zrho(k,j)) )) );
        proj_diam(k,j) = (proj_diam(k,j)*sqrt(1 + (((2*H.target)/Ld(k,j))^2) * (( (exp((z_break(k)-proj_alt(k,j))/(2*H.target(1)))) -1) ^2)));
        
        % Airburst
        panfact(k,j) = proj_diam(k,j)/FLAVRS(k,2);
        if panfact(k,j) >= 5
            if panfact(k,j) > panfact(k,j-1)
                z_burst(k,j) = proj_alt(k,j);
            end
        end
        
    end % Ps(k,j) > Yi(k,j)
    
else % proj_diam <= 1km
    
% No fragmentation 
panfact(k,j) = panfact(k,j-1); 
proj_frag(k,j) = proj_frag(k,j-1);    

end % proj_diam <= 1km

else

% if flag.BREAK == 0
    proj_frag(k,j) = proj_frag(k,j-1);
    
end % flag.BREAK

% ---------------------------------------------------------
% Effective area

   proj_xarea(k,j) = (flag.shape(1)*(proj_frag(k,j)^(1/3)))*((proj_mass(k,j)/FLAVRS(k,5))^(2/3));  
        
% ---------------------------------------------------------
% Time

    proj_time(k,j) = (j-1)*input.step_mod(1);
    
% ---------------------------------------------------------

count_in = count_in + 1;
end % j = 1:input.z_steps % inner loop

% % Correct for plotting purposes to show impact to surface when mass>0
% if proj_alt(k,input.z_steps) > 0 && proj_massr(k,input.z_steps) > 0
%     proj_alt(k,j) = 0;
% end

count_out = count_out + 1;
end % k = 1:n.projectiles % outer loop

% =========================================================
% Crater formation

% Set global values and matrices to be populated, in alphabetical order
% ... (1) column    #    	time step of altitude = 0
% ... (2) ang       deg.    angle of trajectory relative to horizontal
% ... (3) dist      m       distance (x)
% ... (4) frag      #       number of projectile fragments
% ... (5) mass      kg      projectile mass
% ... (6) massr     %       projectile mass ratio to initial
% ... (7) vel       m/s     velocity
% ... (8) rad       m       projectile fragment radius
% ... (9) xarea     m2      effective cross sectional area
% ... (10) time     s       time to reach surface
impact = zeros(n.projectiles, 11);

% Set variables to be populated in 'crater' matrix
% ... (1) Dt        m       transient crater diameter
% ... (2) SCid      ID      simple or complex crater ID (1/simple, 2/complex)
% ... (3) Df        m       final crater diameter
% ... (4) df        m       final crater depth
% ... Calcualte using variables from 'impact' matrix
% ...... (2) ang       deg.    angle of trajectory relative to horizontal
% ...... (4) frag      #       number of projectile fragments
% ...... (7) rad       m       projectile fragment radius
% ...... (8) vel       m/s     velocity
crater = zeros(n.projectiles, 6);

for k = 1:n.projectiles
    
% Pull final projectile characteristics
% --------------------------------------------
if m ~= 4 
    % If m = (1/Venus, 2/Earth, 3/Mars)  
    % Determine time step of impact
    % ... if the projectile alt has stagnated
    if (proj_alt(k,input.z_steps(1)-1)) > 0
        impact(k,1) = input.z_steps(1);
    else
        impact(k,1) = find(proj_alt(k,:) == 0,1);
        % Pull impact characteristics
        impact(k,2) = proj_ang(k,impact(k,1));
        impact(k,3) = proj_dist(k,impact(k,1));
        impact(k,4) = proj_frag(k,impact(k,1));
        impact(k,5) = proj_mass(k,impact(k,1));
        impact(k,6) = proj_massr(k,impact(k,1));
        impact(k,7) = proj_cfa(k,impact(k,1));
        impact(k,8) = proj_vel(k,impact(k,1));
        impact(k,9) = proj_xarea(k,impact(k,1));
        impact(k,10) = proj_time(k,impact(k,1));
        impact(k,11) = proj_diam(k,impact(k,1));
    end

% --------------------------------------------
else   
    % If m == 4 (Moon)  
    % Pull impact characteristics
    impact(k,1) = 1;
    impact(k,2) = FLAVRS(k,11);
    impact(k,3) = 0;
    impact(k,4) = 1;
    impact(k,5) = FLAVRS(k,7);
    impact(k,6) = 100;
    impact(k,7) = FLAVRS(k,2)/2;
    impact(k,8) = FLAVRS(k,4);
    impact(k,9) = FLAVRS(k,1);
    impact(k,10) = 1;
    impact(k,11) = FLAVRS(k,2);
end
% --------------------------------------------

% Pull real numbers only
impact = real(impact);
impact(isnan(impact)) = 0;

% Caclulate transient crater diameter
% ... (1) Dt        m       transient crater diameter
% ... (1) Dt        km      transient crater diameter (converted for depth calcs)
% ....... Calculation walkthrough
% ....... Dt = 1.161 * [ (rho_p/rho_t)^(1/3) * L^0.78 * V^0.44 * g^-0.22 * sin(theta)^(1/3) ] 
% ....... crater(1) = 1.161 * ( (FLAVRS(k,5)/rho_t)^(1/3) * ((impact(7)*2)^0.78) * (impact(8)^0.44) * (g.target^-0.22) * sin(impact(2))^(1/3) );

    % calcualte transient diameter
    crater(k,1) = 1.161 * ( ((FLAVRS(k,5)/rho_t)^(1/3)) * (impact(k,11)^0.78) * (impact(k,8)^0.44) * (g.target^-0.22) * sind(impact(k,2))^(1/3) ) ;
    crater(k,1) = crater(k,1)/1000;
    
% Determine simple or complex crater
% ... SCid = simple or complex crater identification
% ... (2) SCid      ID      simple or complex crater ID (1/simple, 2/complex)

    % Is it simple or complex?
    if crater(k,1) < (Dsc.target(1)/1000)
        crater(k,2) = 1;
    else
        crater(k,2) = 2;
    end % crater(k,1) < Dsc.target(1)

% Caclulate final crater diameter
% ... Df = Calculate transient crater diameterviolinplot
% ... (3) Df        km       final crater diameter
% ... (4) df        km       final crater depth
% ... (5) panfact   -        maximum pancake factor
% ... (6) zburst    km       altitude of airburst
% ... D = 1.17 * ( Dt^1.13 / Dsc.target(1)^0.13 )
   
    if crater(k,2) == 1
        
        % Simple crater depth
        crater(k,3) = (1.25 * crater(k,1));
        a = ( crater(k,1) / (2*sqrt(2)) );
        b = ( 0.07 * ( (crater(k,1)^4) / (1.25*(crater(k,1)))^3) );
        c = ( 2.8 * (0.032 * (1.25*(crater(k,1)))^3) );
        d = ( ( a + ( 0.04 * ( ((crater(k,1)^4)) / ((1.25*(crater(k,1)))^3) ) ) ) / ( a * (1.25*(crater(k,1)))^2 ) );
        crater(k,4) = a+b-(c*d);
        
    else % crater(k,2) == 2
        
        % Complex crater depth
        crater(k,3) = (1.17 * ( crater(k,1)^1.13 / (Dsc.target(1)/1000)^0.13 ));
        crater(k,4) = 1.04*((crater(k,3))^0.301);
        
    end % crater(k,2) == 1 
    
    crater(k,5) = max(panfact(k,:));
    crater(k,6) = min(z_burst(k,:));

% end % impact(k,1) ~= input.z_steps(1)
    
end % n_projectiles

% =========================================================
% Pull statistics

% Pull initial and final characteristics
% ... (1) rho_p         kg/m3   projectile density
% ... (2) init_diam     m       initial diameter
% ... (3) init_mass     kg      initial mass
% ... (4) init_vel      km/s    initial velocity
% ... (5) init_angle    deg     initial degree
% ... (6) final_diam    m       final diameter
% ... (7) final_mass    kg      final mass 
% ... (8) final_massr   %       final mass remain
% ... (9) final_vel     km/s    final velocity
% ... (10) final_angle  deg     final degree
% ... (11) trans_diam   m       transient crater diameter
% ... (12) crater_diam  m       final crater diameter
% ... (13) crater_depth m       final crater depth

stats = horzcat(FLAVRS(:,5),FLAVRS(:,2),FLAVRS(:,7),FLAVRS(:,4),FLAVRS(:,3),...
                impact(:,11), impact(:,5), impact(:,6), impact(:,8)/1000, impact(:,2),...
                (crater(:,1)*1000),(crater(:,3)*1000),(crater(:,4)*1000),...
                crater(:,5),(crater(:,6)));
        
% =========================================================
% Save all target variables

if m == 1 % Venus
    Valt = proj_alt;      Vang = proj_ang;          Vdist = proj_dist;      Vreyn = proj_reyn;
    Vdrag = proj_drag;    Vfrag = proj_frag;        Vlift = proj_lift;      Vmassr = proj_massr;
    Vmass = proj_mass;    Vrad = proj_cfa;          Vvel = proj_vel;        Vtime = proj_time;
    Vxarea = proj_xarea;  Vzrho = proj_zrho;        VFLAVRS2 = FLAVRS;      Vvelt = proj_velt; 
    Vztemp = proj_ztemp;  Vzpress = proj_zpress;    Vdiam = proj_diam;      Vcfa = proj_cfa;
    VPs = Ps;             VYi = Yi;                 VLd = Ld;               VLz = Lz; 
    Vzbreak = z_break;    Vzburst = z_burst;        Vshape = proj_shape;    Vpanfact = panfact; 
    craterV = crater;     impactV = impact;         statsV = stats;
elseif m == 2 % Earth
    Ealt = proj_alt;      Eang = proj_ang;          Edist = proj_dist;      Ereyn = proj_reyn;
    Edrag = proj_drag;    Efrag = proj_frag;        Elift = proj_lift;      Emassr = proj_massr;
    Emass = proj_mass;    Erad = proj_cfa;          Evel = proj_vel;        Etime = proj_time; 
    Exarea = proj_xarea;  Ezrho = proj_zrho;        EFLAVRS2 = FLAVRS;      Evelt = proj_velt;   
    Eztemp = proj_ztemp;  Ezpress = proj_zpress;    Ediam = proj_diam;      Ecfa = proj_cfa;
    EPs = Ps;             EYi = Yi;                 ELd = Ld;               ELz = Lz; 
    Ezbreak = z_break;    Ezburst = z_burst;        Eshape = proj_shape;    Epanfact = panfact;
    craterE = crater;     impactE = impact;         statsE = stats;
elseif m == 3 % Mars
    Malt = proj_alt;      Mang = proj_ang;        Mdist = proj_dist;        Mreyn = proj_reyn;
    Mdrag = proj_drag;    Mfrag = proj_frag;      Mlift = proj_lift;        Mmassr = proj_massr;
    Mmass = proj_mass;    Mrad = proj_cfa;        Mvel = proj_vel;          Mtime = proj_time;
    Mxarea = proj_xarea;  Mzrho = proj_zrho;      MFLAVRS2 = FLAVRS;        Mvelt = proj_velt;
    Mztemp = proj_ztemp;  Mzpress = proj_zpress;  Mdiam = proj_diam;        Mcfa = proj_cfa;  
    MPs = Ps;             MYi = Yi;               MLd = Ld;                 MLz = Lz; 
    Mzbreak = z_break;    Mzburst = z_burst;      Mshape = proj_shape;      Mpanfact = panfact;
    craterM = crater;     impactM = impact;       statsM = stats;
elseif m == 4 % Moon
    Lalt = proj_alt;      Lang = proj_ang;        Ldist = proj_dist;        Lreyn = proj_reyn;
    Ldrag = proj_drag;    Lfrag = proj_frag;      Llift = proj_lift;        Lmassr = proj_massr;
    Lmass = proj_mass;    Lrad = proj_cfa;        Lvel = proj_vel;          Ltime = proj_time;
    Lxarea = proj_xarea;  Lzrho = proj_zrho;      LFLAVRS2 = FLAVRS;        Lvelt = proj_velt;
    Lztemp = proj_ztemp;  Lzpress = proj_zpress;  Ldiam = proj_diam;        Lcfa = proj_cfa;      
    craterL = crater;     impactL = impact;       statsL = stats;
end

% =========================================================
% End looping

end % planetary target loop
