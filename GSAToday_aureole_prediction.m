%file: GSAToday_aureole_prediction.m

%author:  Evan J. Ramos
%date:    20 Aug 2019

%Description: This script using a 2D reactive transport model to predict
%the volume of a metamorphic aureole as a function of host rock
%permeability and intrusion volumes.

close all
clear
clc

%% Simulation parameters

filename                   = 'mineralisotopeconstant.txt';

%% Dimensional constants

%Space
%%%%
Depth    = 5000;     %[m]
Length   = 5000;    %[m] 8000
Width    = 1;        %[m] z-direction width

%Time
%%%%
yrs2s    = 3600*24*365; %conversion factor to seconds from years
dt       = 1e11; %1e11 %1e10      %[s] time step
times    = 0:dt:100*dt; %40%320 %[s] vector of total steps


%Darcy
%%%%
Mu       = 8.9e-4;      %[Pa*s] Dynamic viscosity of water at 20 deg C %%%%%%%%%%
rho_f    = 1000;        %[kg/m^3] Density of water
rho_s    = 2800;        %[kg/m^3] Density of solid
rho_m    = 1500;        %[kg/m^3] Density of magmatic fluid
grav     = 9.81;        %[m/s^2] Acceleration due to gravity on Earth
tensile  = 20e6;        %[Pa] average tensile strength of limestone (Mark Zoback)
phi_min  = 1e-4;        %[m^3/m^3] minimum imposed porosity
phi_0    = 5e-3;        %[m^3/m^3] reference porosity
k_0      = 1e-17;       %[m^2] reference permeability
n        = 3;           %integer exponent for power law calc
kvkh     = 1;


%Temperature
%%%%
T_s      = 293;         %[K]
qt_s     = .06;         %[W/m^2] surface heat flux 
A        = 2e-6;        %[W/m^3] radiogenic heat production: average
kf       = 0.65;        %[W/m/K] thermal conductivity of water
ksol     = 2.08;        %[W/m/K] thermal conductivity of rock
cp_f     = 1046;        %[J/(kg K)] specific heat of water
cp_s     = 908;         %[J/(kg K)] specific heat of solid

%pluton
phi_p    = .0025;             %[m^3/m^3] porosity of pluton
rho_p    = 2800;              %[kg/m^3] pluton density
cp_p     = 1450;              %[J/kg/K] specific heat 
Tp       = 1355;%1173;              %[K] pluton temperature
Tliq     = 1355;%1173;              %[K] liquidus temperature
Tsol     = 983;               %[K] solidus temperature
t_inj    = dt;                %[s] timescale over which pluton injects its fluids
wt_perc  = .0089;             %[1] weight percent volatiles in magma
Q_gamma  = wt_perc*rho_m/(rho_f*t_inj); %[kg/m^3/s] fluid source injection rate 
Q_latent = 190e3;             %[J/kg] total latent heat of crystallization
dT       = Tliq - Tsol;       %[K] temperature range of crystallization (Hanson, 1995)

dHdT_l   = Q_latent/dT;       %[J/kg/K] latent heat of crystallization pluton
melt     = @(T) (T - Tsol)/dT;%[1] function that computes melt fraction
T_crit   = 773;               %[K] temperature at which pluton becomes permeable
ho_f     = 117390;            %[J/kg] specific enthalpy of water at STP

%Stable isotopes
%%%%
tau      = sqrt(2);        %[m/m] tortuosity
vphi_m   = 1;              %[m^3/m^3] volume fraction of mineral in solid part
Do       = 1e-9;           %[m^2/s] solute diffusivity in water (Marc 2016)
delta_s  = -10;            %[per mil] surface delta 18O value
delta_p  = 6.4;              %[per mil] plutonic delta 18O value (D'Errico 2010)
R_std    = .0020052;       %[1] 18O/16O of VSMOW
R_gamma  = delta_p/1000*R_std + R_std;
xw       = .89;            %[kg/kg] mass ratio of oxygen to total molar mass

%convert deltas to ratios
Rwi      = (delta_s/1000)*R_std + R_std;

%% Build Grid and operators

Nx = 200; Ny = 50;
% Defining the computational grid and conductivities
Grid.xmin = 0; Grid.xmax = Length;  Grid.Nx = Nx;
Grid.ymin = 0; Grid.ymax = Depth;   Grid.Ny = Ny; 
Grid.Nz = 1;
Grid.dz = Width;
Grid.psi_dir = 'xy';
Grid      = build_grid(Grid);
[Xc,Yc]   = meshgrid(Grid.xc,Grid.yc);

%Operators
[D,G,I]           = build_ops(Grid);

%% Characteristic scaling terms

%Topography driven flow
H      = Depth;         %[m]
dh     = .1*H;          %[m]

%% Boundary conditions

% Parameters: pressure
Param.p.dof_dir   = Grid.dof_ymax;
Param.p.dof_f_dir = Grid.dof_f_ymax;
slope             = linspace(0,dh,Grid.Nx/2)';
Param.p.g         = rho_f*grav*(0*Depth + [slope; flipud(slope)]);
Param.dof_out     = Param.p.dof_dir;
Param.p.dof_neu   = [];
Param.p.dof_f_neu = [];
Param.p.qb        = [];

%% Build boundaries

%Pressure
[Bp,Np,fn_p] = build_bnd(Param.p,Grid);

%% Create Medium
Mineral                    = build_mineral_new(filename,Grid);     

%area, permeability, and porosity ranges
areas                      = linspace(100,5000,20).^2;
ks                         = logspace(-15,-17,3);
phis                       = poro_from_perm(k_0,ks,phi_min,phi_0,n);
FS_PS                      = zeros(length(times),length(areas) + length(ks));

%% Loop through permeabilities
%decarbonation amounts and rates for constant pluton size

Fracs.frac = cell(length(areas),length(ks));
Fracs.plut = cell(length(areas),length(ks));

for j = 1:length(ks)
    for kk = 1:length(areas)
    %% Create the rest of the medium
    %define pluton location
    pluton_len   = sqrt(areas(kk));
    pluton_wid   = pluton_len;
    Pluton.dof   = Grid.dof( Xc > (Grid.xmax/2 - pluton_wid/2) & ...
                             Xc < (Grid.xmax/2 + pluton_wid/2) & ...
                             Yc > Grid.ymin                    & ...
                             Yc < pluton_len);
    Da           = D(Pluton.dof,:);
    Pluton.dof_f = Grid.dof_f(abs(sum(Da,1)) > eps);
    
    %define mineral attributes
    Mineral.cal.phi  = .75*ones(Grid.N,1); Mineral.cal.phi(Pluton.dof)  = 0;
    Mineral.dol.phi  = .25*ones(Grid.N,1); Mineral.dol.phi(Pluton.dof)  = 0;
    Mineral.qtz.phi  = 0*ones(Grid.N,1);   Mineral.qtz.phi(Pluton.dof)  = .20;
    Mineral.An40.phi = 0*ones(Grid.N,1);   Mineral.An40.phi(Pluton.dof) = .58;
    Mineral.mus.phi  = 0*ones(Grid.N,1);   Mineral.mus.phi(Pluton.dof)  = 0;
    Mineral.bte.phi  = 0*ones(Grid.N,1);   Mineral.bte.phi(Pluton.dof)  = .15;
    Mineral.hbd.phi  = 0*ones(Grid.N,1);   Mineral.hbd.phi(Pluton.dof)  = .07;
    
    %define iteration porosity and permeability
    phi              = phis(j)*ones(Grid.N,1);
    phi(Pluton.dof)  = phi_p;
    k                = ks(j)*ones(Grid.N,1);
    k(Pluton.dof)    = perm_from_poro(k_0,phi_p,phi_min,phi_0,n);
    k                = reshape(k,Grid.Ny,Grid.Nx);
    Kd               = comp_mean(k,-1,kvkh,Grid);
    phi_m            = vphi_m*(1-phi);
    phi_lim_d        = phi(Grid.dof_xmin);
    
    %Thermal properties
    %Heat transport
    kbar   = phi/tau*kf + (1-phi)*ksol;
    cp_s_i = cp_s + dHdT_l;
    rhobar = phi*rho_f*cp_f + (1-phi)*rho_s*cp_s;

    rhobar(Pluton.dof) = phi(Pluton.dof)*rho_f*cp_f + ...
                         (1-phi(Pluton.dof))*rho_p*cp_s_i;

    qt_b   = qt_s - A*H;                             %[W/m^2] basal heat flux
    Teq    = @(z) (A/kbar(Grid.dof_ymin(1)))*(H*z - .5*(H^2 + z.^2)) + ...
                  qt_s/kbar(Grid.dof_ymin(1))*(H-z) + T_s; %[K]
    
    %% Initialize pressure, temperature, and fluid stable isotope values
    %pressure
    p_litho          = reshape(rho_s*grav*(H-Yc),Grid.N,1);
    p                = p_litho/rho_s*rho_f;
    fs_ps            = zeros(length(times),1);

    %Heat transport operators
    theta_a          = 1;                          %advection solved explicitly
    theta_d          = 0;                          %diffusion solved implicitly
    Ts               = zeros(Grid.N,length(times));
    Ts(:,1)          = reshape(repmat(Teq(Grid.yc),1,Grid.Nx),Grid.N,1);
    Ts(Pluton.dof,1) = Tp;
    fs_T             = A*ones(Grid.N,1);
    fs_T(Pluton.dof) = 0;


    %Stable isotope transport operators, initialized
    theta              = 0;
    [~,~,~,del_w,del_m]= make_avg_rxn_terms(Grid,Mineral,Ts(:,1));
    Rs                 = zeros(2*Grid.N,length(times) + 1);
    Rs(1:Grid.N,1)     = del_w/1000*R_std + R_std;       %Rw
    Rs(Grid.N+1:end,1) = del_m/1000*R_std + R_std;       %Rm
    delta_ms           = zeros(Grid.N,length(times));
    delta_ms(:,1)      = del_m;
    fs_O               = zeros(2*Grid.N,1);

    %diagonalized matrices
    Phi                = comp_mean(reshape(phi,Grid.Ny,Grid.Nx),-1,kvkh,Grid);
    Phi_diag           = spdiags(phi,0,Grid.N,Grid.N);
    Phi_m              = spdiags(phi_m,0,Grid.N,Grid.N);
    Rhobar             = spdiags(rhobar,0,Grid.N,Grid.N);%%%%%%%%%%
    Kbar               = comp_mean(reshape(kbar,Grid.Ny,Grid.Nx),1,kvkh,Grid); %%%%%%%%%%% 10 Aug 2017: arithmetic mean instead of harmonic mean
    Xw                 = spdiags(xw*ones(Grid.N,1),0,Grid.N,Grid.N);
    Rho_f              = rho_f;
    KdMu               = Kd/Mu;

        for i = 2:length(times)
            %% Solve steady state pressure --> Darcy Flux
            %pressure
            Lp             = -D*KdMu*G;
            fs_grav        = D*(KdMu*Rho_f*grav)*G*Yc(:);
            fs_p           = zeros(Grid.N,1);
            %Gamma
            if i > 2
                %tracking how much crystallization has occurred between the
                %subsequent time steps so that the appropriate scale to the source
                %term can be applied to the magmatic fluid input.
                diff_melt_i        = melt(Ts(Pluton.dof,i-1));
                diff_melt_ii       = melt(Ts(Pluton.dof,i-2));
                diff_melt_i(diff_melt_i < 0) = 0;
                diff_melt_ii(diff_melt_ii < 0) = 0;
                diff_melt = diff_melt_ii - diff_melt_i;
                diff_melt(diff_melt<0) = 0;
                fs_p(Pluton.dof) = fs_p(Pluton.dof) + Q_gamma*diff_melt;
                fs_T(Pluton.dof) = fs_p(Pluton.dof).*new_rho(Pluton.dof)*ho_f;
                fs_O(Pluton.dof) = fs_p(Pluton.dof)*R_gamma*xw;
                fs_ps(i-1) = sum(fs_p);
            end

            fp             = fs_p + fn_p + fs_grav;
            p             = solve_lbvp(Lp,fp,Bp,Param.p.g,Np); 

            q              = comp_flux_p(D,KdMu,Rho_f*grav*G*Yc(:),G,p,fs_p,Grid,Param.p);

            qtop           = q(Grid.dof_f_ymax);
            qdown_dof      = Grid.dof_ymax(qtop < 0);
            qdown_dof_f    = Grid.dof_f_ymax(qtop < 0);

            %1st order upwinding of fluxes
            Aq = flux_upwind(q,Grid);

            %% Build parameters

            %Molten regions identified
            molten_dof        = find(Ts(:,i-1) > (Tp - dT));
            molten_dof        = intersect(molten_dof,Pluton.dof);
            Da                = D(molten_dof,:);
            molten_dof_f      = (Grid.dof_f(abs(sum(Da,1)) > eps))';
            % Parameters: temperature
            Param.T.dof_dir   = qdown_dof;
            Param.T.dof_f_dir = qdown_dof_f;
            Param.T.g         = T_s*ones(length(qdown_dof),1);
            Param.T.dof_neu   = Grid.dof_ymin(~ismember(Grid.dof_ymin,Pluton.dof));
            Param.T.dof_f_neu = Grid.dof_f_ymin(~ismember(Grid.dof_f_ymin,Pluton.dof_f));
            Param.T.qb        = qt_b;
            %Parameters: stable isotope
            Param.O_w.dof_dir   = [qdown_dof; molten_dof];
            Param.O_w.dof_f_dir = [qdown_dof_f; molten_dof_f];
            Param.O_w.g         = [Rwi*ones(length(qdown_dof),1); R_gamma*ones(length(molten_dof),1)];
            Param.O_w.dof_neu   = [];
            Param.O_w.dof_f_neu = [];
            Param.O_w.qb        = [];

            %Stable isotope
            [BO,NO,fn_O] = build_bnd(Param.O_w,Grid,'coupled');
            %Temperature
            [BT,NT,fn_T] = build_bnd(Param.T,Grid);

            %% Solve for temperature
            % Build heat operators
            Im_T             = @(theta_a,theta_d) ...
                             Rhobar + dt*(1-theta_a)*D*(rho_f*cp_f*Aq) ...
                                    - dt*(1-theta_d)*D*(Kbar*G); %%%%%%
            Ex_T             = @(theta_a,theta_d) ...
                             Rhobar - dt*theta_a*D*(rho_f*cp_f*Aq)...
                                    + dt*theta_d*D*(Kbar*G);     %%%%%%

            Ts(:,i)    = solve_lbvp(Im_T(theta_a,theta_d),...
                                    Ex_T(theta_a,theta_d)*Ts(:,i-1) + ...
                                    dt*(fs_T + fn_T),BT,Param.T.g,NT);

            %% Solve for stable isotopes
            %update diagonalized matrices
            [Alpha,Rk,Xm]  = make_avg_rxn_terms(Grid,Mineral,Ts(:,i));

            %create stable isotope operators
            %dRw/dt PDE
            Im_ww          = @(theta) Phi_diag*Xw...
                             + (1-theta)*dt*...
                             (D*(Aq*xw - Phi/tau*xw*Do*G) + ...
                             (Phi_m*Xm*Rk*Alpha));
            Im_wm          = @(theta) -dt*(1-theta)*(Phi_m*Xm*Rk);
            Ex_ww          = @(theta) Phi_diag.*Xw...
                             - theta*dt*...
                             (D*(Aq*xw - Phi/tau*xw*Do*G) - ...
                             (Phi_m*Xm*Rk*Alpha));
            Ex_wm          = @(theta) theta*dt*(Phi_m*Xm*Rk);

            %dRm/dt ODE
            Im_mw          = @(theta) -dt*(1-theta)*Rk*Alpha;
            Im_mm          = @(theta) I + dt*(1-theta)*Rk;
            Ex_mw          = @(theta) dt*theta*Rk*Alpha;
            Ex_mm          = @(theta) I - dt*theta*Rk;

            %coupled
            Im             = [Im_ww(theta) Im_wm(theta); Im_mw(theta) Im_mm(theta)];
            Ex             = [Ex_ww(theta) Ex_wm(theta); Ex_mw(theta) Ex_mm(theta)];

            %solve stable isotope values
            Rs(:,i)        = solve_lbvp(Im,Ex*Rs(:,i-1) + dt*(fs_O),...
                                        BO,Param.O_w.g,NO);
            delta_ms(:,i)  = 1000*(Rs(Grid.N+1:end,i)-R_std)/R_std;
            
            %% Update T-dependent variables
            %update fluid density
            [new_rho,~,~] = steam_table(Ts(:,i),p,Grid);
            Rho_f         = comp_mean(reshape(new_rho,Grid.Ny,Grid.Nx),1,kvkh,Grid);

            %update source terms, diagonalized matrices
            solid            = find(Ts(:,i) < T_crit);
            solid            = intersect(solid,Pluton.dof);
            fs_frac          = pluton_source(Ts(Pluton.dof,i));
            solidus          = find(Ts(:,i) < (750 + 273));
            rhobar(solidus)  = phi(solidus)*rho_f*cp_f + (1-phi(solidus))*rho_s*cp_s;
            Rhobar           = spdiags(rhobar,0,Grid.N,Grid.N);
            %change permeability of cooling pluton
            [~,solid_loc]    = ismember(Yc(solid),Grid.yc);
            phi(solid)       = phi_lim_d(solid_loc); %porosity of limestone
            k(solid)         = perm_from_poro(1e-17,phi(solid),1e-4,5e-3,3);
            Kd               = comp_mean(k,-1,kvkh,Grid);
            KdMu             = Kd/Mu;
            Phi              = comp_mean(reshape(phi,Grid.Ny,Grid.Nx),-1,1,Grid);
            Phi_diag         = spdiags(phi,0,Grid.N,Grid.N);

            %STORE ASSIMILATION FLUX
            fprintf('\n%d of %d iterations\n',i-1,length(times)-1)
        end
    FS_PS(:,j+length(areas)) = fs_ps;

    %% Percent skarn decarbonated

    aureole = Grid.dof(~ismember(Grid.dof,Pluton.dof));

    %finding d18O_m values below the maximum value for the decarbonated
    %region (think d13C vs. d18O plot from Valley 1986)
    max_d18O_m     = 8; %[permil VSMOW]
    M_CaCO3        = .1; %[kg] molar mass of CaCO3
    M_C            = .012; %[kg] molar mass of C
    
    d_bounds.skarn = [max_d18O_m max(delta_ms(:))];
    d_bounds.therm = [500 700] + 273;
    d_styles.skarn = 'kinetic';
    d_styles.therm = 'linear';
    [skarn_frac, ~, aur_decarb] = decarb_perc_aureole(delta_ms,...
                                                      Ts,...
                                                      d_styles,...
                                                      d_bounds);
    Fracs.frac{kk,j} = aur_decarb;
    Fracs.plut{kk,j} = Pluton;
    
    fprintf('\n\nLoop %d of %d complete. \n\n',...
           (j-1)*length(areas) + kk,...
           length(ks)*length(areas))
    end
end