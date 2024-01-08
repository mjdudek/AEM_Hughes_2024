%% FLAVRS Builder (for AEM)
% Last update: 2024-01-05

% Developed by: Marissa Dudek (1)
% (1) University of North Carolina at Chapel Hill, Geological Sciences

% Input projectile parameters
% ... Diameter (m)
% ... Density (kg/m3) [2700,3300,4500,7800] [ch,st,si,ir]
% ... Mass (kg) 
% ... Velocity (km/sec)
% ... Angle (deg. from horizontal)
% ... Combo? (1/on, 0/off)
Diameter = [14000,50,50,20]; 
Density = [2700,7800,2700,3300]; 
Mass = [0,0,0,0];
Velocity = [17,17,17,17];
Angle = [45,45,45,20]; 
Combo = [0];

% Testing ------------------------------

% % Theortical popultion ---------- 
% Diameter = [1,10,15,20,25,30,35,40,45,50,60,75,100,150,200,250,350,500,750,1000,1500,2500,3500,5000,7500,10000]; 
% Density = [2700,3300,4500,7800]; 
% Mass = 0;
% Velocity = [1,5,10,15,20,25,30,40];
% Angle = [15,30,45,60,75,90]; 
% Combo = [1];

% % Terrestrial Analogs (x4) ---------- 
% Chicxulub, Meteor (Canyon Diablo), Tunguska, Chelyabinsk
% Diameter = [14000,50,50,20]; 
% Density = [2700,7800,2700,3300]; 
% Mass = [0,0,0,0];
% Velocity = [17,17,17,17];
% Angle = [45,45,45,20]; 
% Combo = [0];

% % Terrestrial Analogs (x6) ---------- 
% Tunguska, Chelyabinsk, Barringer, Ries Crater, Chicxulub, South Pole-Aitken Basin
% Diameter = [60,22,40,1500,17500,315000]; 
% Density = [2700,3300,7800,2700,2700,2700]; 
% Mass = [0,0,0,0,0,0];
% Velocity = [20,19,20,20,20,20];
% Angle = [45,20,45,30,45,45]; 
% Combo = [0];

% % Large & Weak, slow/avg/fast ---------- 
% % Output = [78605171.3093685,10000,0,10000,2700,223562.579416875,1.41371669411541e+15,0,0,0,45;78605171.3093685,10000,0,20000,2700,223562.579416875,1.41371669411541e+15,0,0,0,45;78605171.3093685,10000,0,30000,2700,223562.579416875,1.41371669411541e+15,0,0,0,45]
% Diameter = 10000; 
% Density = 2700; 
% Mass = 0;
% Velocity = [10, 20, 30];
% Angle = 45; 
% Combo = [1];

% % Small & Strong, slow/avg/fast -------
% % Output = [786051.713093685,1000,0,10000,8000,48778188.9687869,4188790204786.39,0,0,0,45;786051.713093685,1000,0,20000,8000,48778188.9687869,4188790204786.39,0,0,0,45;786051.713093685,1000,0,30000,8000,48778188.9687869,4188790204786.39,0,0,0,45]
% Diameter = 1000; 
% Density = 8000; 
% Mass = 0;
% Velocity = [10, 20, 30];
% Angle = 45; 
% Combo = [1];

% % Average & Average, slow/avg/fast ---
% % Output = [1965.12928273421,50,0,10000,3300,491594.586164026,215984494.934298,0,0,0,45;1965.12928273421,50,0,20000,3300,491594.586164026,215984494.934298,0,0,0,45;1965.12928273421,50,0,30000,3300,491594.586164026,215984494.934298,0,0,0,45]
% Diameter = 50; 
% Density = 3300; 
% Mass = 0;
% Velocity = [10, 20, 30];
% Angle = 45; 
% Combo = [1];

%% Populate projectile population

% Build matrix of projectiles
if Combo == 1
    % Option 1: build projectiles w/ combination of all characteristics
    projectile = combvec(Diameter, Density, Mass, Velocity, Angle).';
    n.projectiles = size(projectile,1); 
else % If Combo == 0
    % Option 2: build projectiles w/ only user inputs
    projectile = [Diameter; Density; Mass; Velocity; Angle].';
    n.projectiles = size(projectile,1); 
end

% Build FLAVRS
FLAVRS = zeros(n.projectiles, 11);
for k = 1:n.projectiles  
    % ---------------------------------------------------------
    % (Mass and Density) or (Diameter and Density)?
    if projectile(k,1) == 0
        projectile(k,1) = nthroot((projectile(k,3)/((4/3)*(pi)*projectile(k,2))),3);
    end
    if projectile(k,3) == 0
        projectile(k,3) = ((4/3)*(pi)*((projectile(k,1)/2)^3)*projectile(k,2)); 
    end
    % ---------------------------------------------------------
    % (:,2) L = projectile diameter (m)
    FLAVRS(k,2) = projectile(k,1); 
    % ---------------------------------------------------------
    % (:,3) A = planetary position
    FLAVRS(k,3) = 0;
    % ---------------------------------------------------------
    % (:,4) V = velocity at top of atmosphere (m/s)
    FLAVRS(k,4) = projectile(k,4)*1000;   
    % ---------------------------------------------------------
    % (:,5) R = projectile density (rho)  
    FLAVRS(k,5) = projectile(k,2);
    % ---------------------------------------------------------
    % (:,6) S = projectile strength   
    FLAVRS(k,6) = 10^(2.107+(0.0624*sqrt(FLAVRS(k,5)))); 
    % ---------------------------------------------------------
    % (:,7) m = projectile mass
    FLAVRS(k,7) = ((4/3)*(pi)*((FLAVRS(k,2)/2)^3)*FLAVRS(k,5)); 
    % ---------------------------------------------------------
    % (:,1) F = effective cross-sectional area
    % Note: Make sure shape.f matches the shape.f in AEM
    shape.f = 1.21; 
    FLAVRS(k,1) = shape.f*((FLAVRS(k,7)/FLAVRS(k,5))^(2/3));  
    % ---------------------------------------------------------
    % (:,8) # = tracking number from MDLPVCA
    FLAVRS(k,8) = 0; 
    % ---------------------------------------------------------
    % (:,9) p = projectile probability of occurance
    FLAVRS(k,9) = 0; 
    % ---------------------------------------------------------
    % (:,10) c = LPF crater bin
    FLAVRS(k,10) = 0;
    % ---------------------------------------------------------
    % (:,11) a = angle of trajectory
    FLAVRS(k,11) = projectile(k,5);
    % ---------------------------------------------------------
end

FLAVRS

% Display properties
%sprintf("Diameter = " +FLAVRS(2)/1000+ " km \nDensity = " +FLAVRS(5)+ " kg/m3 \nMass = " +FLAVRS(7)+ " kg \nVelocity = " +FLAVRS(4)/1000+ " km/s \nAngle = " +FLAVRS(11)+ " deg")

