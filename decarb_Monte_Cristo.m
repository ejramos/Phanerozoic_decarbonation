function [sed_vols, f_decarb, f_aureole, F_CO2_assim, F_CO2_meta] = ...
                                       decarb_Monte_Cristo(areal_addition,...
                                                           ages,...
                                                           unit_vols,...
                                                           depth,...
                                                           time_binwidth,...
                                                           props,...
                                                           f_stoic,...
                                                           burial_rte,...
                                                           n_iterations,...
                                                           t_range,...
                                                           f_int,...
                                                           f_aur)
%%

%author: Evan J. Ramos
%date:   22 Aug 2019

%description: Using rock information gathered from Macrostrat and areal
%addition rates of granitoids from Cao et al. (2017) EPSL, this function
%performs a Monte Carlo analysis to compute a range of sediment volumes and
%concomitant volumes rock reacted and CO2 fluxes.

%Inputs:
%       areal_addition: Nbins x 1 vector of areal addition rates
%                       [1000 km^2/Myr]

%       ages: N x 1 vector of rock unit ages [Myr]

%       unit vols: N x 2 vector of erosion-corrected rock unit volumes 
%                  col 1: minimum vol [km^3], col: 2 maximum vol [km^3]

%       depth: scalar regarding depth over which magma is emplaced [km]

%       time_binwidth: scalar for time binwidth related to areal_addition
%                      [Myr]

%       props: structures containing proportion of rock types for each rock 
%              unit
%              Fields: carb (carbonate), silc (siliciclastic)

%       f_stoic: 3 x 1 vector of structures containing the range of
%                stoichiometric fractions that each rock type emits
%                element 1: carbonate; element 2: siliciclastic;
%                element 3: mixed;
%                Fields: min, max

%       burial_rte: 2 x 1 vector containing burial rates
%                   over which the magma intrudes [Myr]
%                   element 1: min rate [km/Myr]; element 2: max lag [km/Myr]

%       n_iterations: scalar for number of Monte Carlo iterations

%       t_range: 2 x 1 vector containing range of times over which data 
%                applies
%                element 1: t_bot [Myr]; element 2: t_top [Myr]

%       f_int/f_aur: vectors containing empirical estimates that relate
%                    intrusion volume to aureole volume

%Ouputs:
%       sed_vols: Structures containing Nbins x n_iterations matrix of sed 
%                 volumes for each rock type (each field) [km^3]
%                 Fields: carb (carbonate), silc (siliciclastic), mixd
%                         (mixed), total

%       f_decarb: Nbins x n_iterations matrix of volume/domain ratios

%       f_aureole: Nbins x n_iterations matrix of aureole/domain ratios

%       F_CO2_assim: Nbins x n_iterations matrix of assimilation CO2 fluxes 
%                    [Tg/yr]

%       F_CO2_meta: Nbins x n_uterations matrix of contact metamorphic CO2
%                   fluxes [Tg/yr]

%% Create bin width time intervals and preallocate ouputs

%time range for histogram estimates
range = (t_range(2)+time_binwidth/2):time_binwidth:(t_range(1)-time_binwidth/2);
Nbins = length(range);

%preallocated output matrices
sed_vols.carb = zeros(Nbins,n_iterations);
sed_vols.silc = zeros(Nbins,n_iterations);
sed_vols.mixd = zeros(Nbins,n_iterations);
sed_vols.plut = zeros(Nbins,n_iterations);
sed_vols.totl = zeros(Nbins,n_iterations);
f_decarb      = zeros(Nbins,n_iterations);
f_aureole     = zeros(Nbins,n_iterations); 
F_CO2_assim   = zeros(Nbins,n_iterations);
F_CO2_meta    = zeros(Nbins,n_iterations);

%number of individual rock units
Nunits = length(ages);

%% Monte Carlo approximation performed at each time bin

for i = 1:n_iterations
    %step 0: generate a random burial rate w/in range of rates + compute
    %        time lags, given depth of intrusion and burial rate
    rte_burial = (burial_rte(2)-burial_rte(1))*rand + burial_rte(1);
    min_depth  = 3; %[km] minimum depth of skarn formation
    max_depth  = min_depth + depth;
    
    time_lags  = [min_depth max_depth]/rte_burial;
    
    %step 1: generate range rock volumes, within their min and max
    rock_vol_i = (unit_vols(:,2)-unit_vols(:,1)).*rand(Nunits,1) + ...
                  unit_vols(:,1);
    carb_vol_i = props.carb.*rock_vol_i;
    silc_vol_i = props.silc.*rock_vol_i;
    volc_vol_i = props.volc.*rock_vol_i;
    plut_vol_i = props.plut.*rock_vol_i;
    meta_vol_i = props.plut.*rock_vol_i;
    
    %step 2: generate cumulative volumes for each time bin
    carb_vol_cum_i = zeros(Nbins,1);
    silc_vol_cum_i = zeros(Nbins,1);
    mixd_vol_cum_i = zeros(Nbins,1);
    volc_vol_cum_i = zeros(Nbins,1);
    plut_vol_cum_i = zeros(Nbins,1);
    meta_vol_cum_i = zeros(Nbins,1);
    totl_vol_cum_i = zeros(Nbins,1);
    
    for j = 1:Nbins
        select = ages >= range(j)-time_binwidth/2 & ages < range(j)+time_binwidth/2;
        carb_sel = props.carb' > 0 & props.silc' == 0;
        silc_sel = props.silc' > 0 & props.carb' == 0;
        mixd_sel = props.carb' > 0 & props.silc' > 0;
        volc_sel = props.volc' > 0;
        plut_sel = props.plut' > 0;
        meta_sel = props.meta' > 0;
        
        carb_vol_cum_i(j) = sum(carb_vol_i(select & carb_sel));
        silc_vol_cum_i(j) = sum(silc_vol_i(select & silc_sel));
        mixd_vol_cum_i(j) = sum(carb_vol_i(select & mixd_sel)) + ...
                            sum(silc_vol_i(select & mixd_sel));
        volc_vol_cum_i(j) = sum(volc_vol_i(select & volc_sel));
        plut_vol_cum_i(j) = sum(plut_vol_i(select & plut_sel));
        meta_vol_cum_i(j) = sum(meta_vol_i(select & meta_sel));
        totl_vol_cum_i(j) = sum(rock_vol_i(select)) - ...
                            (volc_vol_cum_i(j) + plut_vol_cum_i(j) + meta_vol_cum_i(j));
    end
    
    %step 3: generate random stoichiometric decarbonation fractions
    f_carb_i = (f_stoic(1).max - f_stoic(1).min)*rand + f_stoic(1).min;
    f_silc_i = (f_stoic(2).max - f_stoic(2).min)*rand + f_stoic(2).min;
    f_mixd_i = (f_stoic(3).max - f_stoic(3).min)*rand + f_stoic(3).min;
    f_asim_i = (f_stoic(4).max - f_stoic(4).min)*rand + f_stoic(4).min;
    
    %step 4: compare rock volumes for intruded magma vs. sediment + time
    %        lag, then compute original CO2 flux estimate
    f_decarb_i    = zeros(Nbins,1);
    f_aureole_i   = zeros(Nbins,1);
    F_CO2_assim_i = zeros(Nbins,1);
    F_CO2_meta_i  = zeros(Nbins,1);
    
    rho_sed = 2700; %[Tg/km^3] average density of sedimentary rock
    for j = 1:Nbins
        
        %update select to look for rocks with ages within time lag 
        select = ages >= (range(j) + time_lags(1)) & ...
                 ages <= (range(j) + time_lags(2));
             
        %note: carb, silc, and mixd sel remain the same
        carb_vol_cum_j = sum(carb_vol_i(select & carb_sel));
        silc_vol_cum_j = sum(silc_vol_i(select & silc_sel));
        mixd_vol_cum_j = sum(carb_vol_i(select & mixd_sel)) + ...
                            sum(silc_vol_i(select & mixd_sel));
        totl_vol_cum_j = sum(rock_vol_i(select));
        
        %compute ratio of intrusion vol vs. sed vol
        f_decarb_i(j)  = areal_addition(j)*1000*depth*time_binwidth/...
                        totl_vol_cum_j;
        %correct if f_decarb > 1 (impossible)
        if f_decarb_i(j) > 1
            f_decarb_i(j) = 1;
        end
        f_aureole_i(j) = interp1(f_int,f_aur,f_decarb_i(j));
        %correct if f_decarb + f_aureole > 1 (impossible)
        if f_decarb_i(j) + f_aureole_i(j) > 1
            f_aureole_i(j) = 1 - f_decarb_i(j);
        end
                   
        %compute decarbonation flux
        phi     = 0.01;                    %[km^3/km^3] porosity
        V_CO2_j = (1-phi)*(f_carb_i*carb_vol_cum_j + ...
                           f_silc_i*silc_vol_cum_j + ...
                           f_mixd_i*mixd_vol_cum_j); %[km^3] max CO2 volume
        
        F_CO2_assim_i(j) = f_asim_i*rho_sed*f_decarb_i(j)*V_CO2_j/time_binwidth/1e6;
        F_CO2_meta_i(j) = rho_sed*f_aureole_i(j)*V_CO2_j/time_binwidth/1e6;
                    %[Tg/yr] decarbonation flux
    end

    %step 5: store results from monte carlo iteration into matrix
    sed_vols.carb(:,i) = carb_vol_cum_i;
    sed_vols.silc(:,i) = silc_vol_cum_i;
    sed_vols.mixd(:,i) = mixd_vol_cum_i;
    sed_vols.plut(:,i) = plut_vol_cum_i;
    sed_vols.totl(:,i) = totl_vol_cum_i;
    f_decarb(:,i)      = f_decarb_i;
    f_aureole(:,i)     = f_aureole_i;
    F_CO2_assim(:,i)   = F_CO2_assim_i;
    F_CO2_meta(:,i)    = F_CO2_meta_i;
    
    fprintf('\n\nIteration %d of %d complete.\n\n',i,n_iterations)
end
