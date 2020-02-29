%%
%filename: decarb_Monte_Carlo_SNB.m

%author: Evan J. Ramos
%date:   27 Aug 2019

%Description: Much like decarb_Monte_Carlo_script.m, except using
%             information specific to the Sierra Nevada Batholith over a
%             narrower 30 Myr time period.


clc

set(0, 'DefaultAxesFontSize', 15,...
       'DefaultLineLinewidth', 2,...
       'DefaultAxesXColor','k',...
       'DefaultAxesYColor','k')

%% Parameters for the computation

%dimensions for sedimentary outcrops for N. America
%Total arc
arc_length   = [2000 3500]; %[km] North America
arc_width    = [300 400];   %[km] arc widths
arc_depth    = 7;           %[km] total stratigraphic thickness
rho_sed      = 2700;        %[Tg/km^3] sedimentary rock density
phi          = 0.01;        %[1] sedimentary rock porosity

%Paleozoic Tinemahah Reservoir section
Pz.depth     = 3.5;         %[km] Paleozoic stratigraphic thickness
Pz.carb      = 0.23;        %[1] proportion of carbonate
Pz.silc      = 0.77;        %[1] proportion of siliciclastic
Pz.mixd      = 0;           %[1] proportion of mixed

%Triassic-Jurassic sections
Mz.depth     = 3.5;
Mz.carb      = [.14 .34];
Mz.mixd      = [.15 .43];

%surface area addition rates for N America
arc_flux   = [1000 3000]; %[km^2/Myr] Cretaceous addition rate
dt         = 40;          %[Myr] time period of flux

%% Decarbonation parameters

%each value is merely an assumption, but could be better justified if a
%reaction stoichiometry is given
f_carb = [0 .9];  %percentage range devolatilized for carbonates (max is the wollastonite forming reaction)
f_silc = [0 .2]; %percentage range devolatilized for siliciclastics (max is arbitrarily small)
f_mixd = [0 .7]; %percentage range devolatilized for mixed carbonate-siliciclastics
f_magm = [.08 .5]; %percentage range devolatilized after assimilation (estimates from Carter and Dasgupta, 2016)

f_stoic(4).min = f_magm(1);
f_stoic(4).max = f_magm(2);
f_stoic(3).min = f_mixd(1);
f_stoic(3).max = f_mixd(2);
f_stoic(2).min = f_silc(1);
f_stoic(2).max = f_silc(2);
f_stoic(1).min = f_carb(1);
f_stoic(1).max = f_carb(2);

depth  = 7; %[km] depth range over which decarbonation occurs

%% Load empirical intrusion volume-aureole volume estimator

load intrusion_aureole_volume

%% Monte Carlo computation

n_iterations = 50000; %number of iterations for Monte Carlo simulation
sed_vols = zeros(n_iterations,1);
mag_vols = zeros(n_iterations,1);
f_intr = zeros(n_iterations,1);
f_aur = zeros(n_iterations,1);
F_CO2_assim = zeros(n_iterations,1);
F_CO2_meta = zeros(n_iterations,1);

Mz_fcarb = zeros(n_iterations,1);
Pz_fcarb = zeros(n_iterations,1);

for i = 1:n_iterations
    
    %step 1: generate range of rock volumes, within their min and max
    sed_vols(i) = ((arc_length(2)-arc_length(1))*rand + arc_length(1))*...
                  ((arc_width(2)-arc_width(1))*rand + arc_width(1))*...
                  arc_depth;
    
    %step 2: generate range of magma volumes, within their max and min
    mag_vols(i) = (arc_flux(2)-arc_flux(1)*rand + arc_flux(1))*arc_depth*...
                  dt;
              
    %step 3: determine volume fractions
    f_intr(i) = mag_vols(i)/sed_vols(i);
    f_aur(i)  = interp1(f_intrusion,f_aureole,f_intr(i));
    
    %step 3: generate random stoichiometric decarbonation fractions
    f_carb_i = (f_stoic(1).max - f_stoic(1).min)*rand + f_stoic(1).min;
    f_silc_i = (f_stoic(2).max - f_stoic(2).min)*rand + f_stoic(2).min;
    f_mixd_i = (f_stoic(3).max - f_stoic(3).min)*rand + f_stoic(3).min;
    f_asim_i = (f_stoic(4).max - f_stoic(4).min)*rand + f_stoic(4).min;
    
    %step 4: generate average proportions for the sedimentary rock volume
    Mz_carb_i = (Mz.carb(2)-Mz.carb(1))*rand + Mz.carb(1);
    Mz_mixd_i = (Mz.mixd(2)-Mz.mixd(1))*rand + Mz.mixd(1);
    Mz_silc_i = 1 - (Mz_carb_i + Mz_mixd_i);
    
    carb_prop_i = (Pz.depth*Pz.carb + Mz.depth*Mz_carb_i)/arc_depth;
    silc_prop_i = (Pz.depth*Pz.silc + Mz.depth*Mz_silc_i)/arc_depth;
    mixd_prop_i = (Pz.depth*Pz.mixd + Mz.depth*Mz_mixd_i)/arc_depth;
    
    %Compute max CO2 volume and associated fluxes
    V_CO2_i = (1-phi)*sed_vols(i)*(f_carb_i*carb_prop_i + ...
                                   f_silc_i*silc_prop_i + ...
                                   f_mixd_i*mixd_prop_i); %[km^3] max CO2 volume
    Pz_fcarb(i) = f_carb_i*.5*Pz.carb + f_mixd_i*.5*Pz.mixd + f_silc_i*.5*Pz.silc;
    Mz_fcarb(i) = f_carb_i*.5*Mz_carb_i + f_mixd_i*.5*Mz_mixd_i + f_silc_i*.5*Mz_silc_i;
    tot = Pz_fcarb(i) + Mz_fcarb(i);
    Pz_fcarb(i) = Pz_fcarb(i)/tot;
    Mz_fcarb(i) = Mz_fcarb(i)/tot;
    
    F_CO2_assim(i) = f_asim_i*rho_sed*f_intr(i)*V_CO2_i/dt/1e6;
    F_CO2_meta(i) = rho_sed*f_aur(i)*V_CO2_i/dt/1e6;
                    %[Mt/yr] decarbonation flux
    
    fprintf('\n\nIteration %d of %d complete.\n\n',i,n_iterations)
end
                                                       
%% Print results:

clc
fprintf('Results: \n\n')
fprintf('Intrusion volume fraction: %.2f +- %.2f \n',...
         mean(f_intr),std(f_intr))
fprintf('Aureole volume fraction: %.2f +- %.2f \n',...
         mean(f_aur),std(f_aur))
fprintf('Assimilation flux: %.1f +- %.1f [Tg/yr]\n',...
         mean(F_CO2_assim),std(F_CO2_assim))
fprintf('Aureole decarbonation flux: %.1f +- %.1f [Tg/yr]\n',...
         mean(F_CO2_meta),std(F_CO2_meta))
fprintf('Net continental arc metamorphic flux: %.1f +- %.1f [Tg/yr]\n\n',...
         mean(F_CO2_meta) + mean(F_CO2_assim),...
         sqrt(mean(F_CO2_assim)^2 + mean(F_CO2_meta)^2))
fprintf('Percentage of CO2 from Mesozoic units: %.1f +- %.1f\n\n',...
         mean(Mz_fcarb)*100,std(Mz_fcarb)*100)
fprintf('Percentage of CO2 from Paleozoic units: %.1f +- %.1f\n\n',...
         mean(Pz_fcarb)*100,std(Pz_fcarb)*100)