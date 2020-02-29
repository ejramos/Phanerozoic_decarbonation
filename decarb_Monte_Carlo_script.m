%%
%filename: decarb_Monte_Carlo_script.m

%author: Evan J. Ramos
%date:   22 Aug 2019

%Description: This script amasses sedimentary rock information from
%             Macrostrat, magma flux information from Cao et al. 2017 EPSL,
%             and monte carlo analyses to assess metamorphic decarbonation
%             fluxes.


close all
clear
clc

set(0, 'DefaultAxesFontSize', 15,...
       'DefaultLineLinewidth', 2,...
       'DefaultAxesXColor','k',...
       'DefaultAxesYColor','k')

%% Upload data from Macrostrat

options = weboptions('Timeout', 30); %it takes a long time to query these
                                     %data online, so we establish a longer
                                     %timeout time (30 seconds)
url = 'https://macrostrat.org/api/v2/units?age_top=0&age_bottom=600&response=long';

Phan = webread(url,options);
Phan = Phan.success.data;

%% Toggle computation of N. America or Global flux

opt = 1;
if opt == 1 %N. America
    filename = 'Cao_et_al_2017_surf_area_NAmerica_uncorrected.txt';
    filestring = 'Granitoid: N. America';
elseif opt == 2 %Global
    filename = 'Cao_et_al_2017_surf_area_Global_uncorrected.txt';
    filestring = 'Granitoid: Global';
end

%% Compute proportion of carbonate vs. siliclastic

sed_type = {Phan.lith};
carbs = 'carbonate';
silcs = 'siliciclastic';
pluts = 'plutonic';
volcs = 'volcanic';
metas = 'meta';

carb_prop = zeros(length(sed_type),1);
silc_prop = zeros(length(sed_type),1);
volc_prop = zeros(length(sed_type),1);
plut_prop = zeros(length(sed_type),1);
meta_prop = zeros(length(sed_type),1);

classes = cell(69130,1); %arbitrarily long cell area to store rock type
count = 1;

for i = 1:length(sed_type)
    type = {sed_type{i}.type};
    prop = [sed_type{i}.prop];
    carb_prop(i) = sum(prop.*strcmpi(type,carbs));
    silc_prop(i) = sum(prop.*strcmpi(type,silcs));
    volc_prop(i) = sum(prop.*strcmpi(type,volcs));
    plut_prop(i) = sum(prop.*strcmpi(type,pluts));
    
    for j = 1:length(type)
        meta_prop(i) = meta_prop(i) + prop(j)*~isempty(strfind(type{j},metas));
        classes{count} = type{j};
        count = count + 1;
    end
end
classes = unique(classes);

%remove excess proportion (1 is max)
carb_prop(carb_prop > 1) = 1;
silc_prop(silc_prop > 1) = 1;
volc_prop(volc_prop > 1) = 1;
plut_prop(plut_prop > 1) = 1;
meta_prop(meta_prop > 1) = 1;

%store to structure
props.carb = carb_prop;
props.silc = silc_prop;
props.volc = volc_prop;
props.plut = plut_prop;
props.meta = meta_prop;

%% Isolate time of interest + making discrete time intervals

t_top = 0;
t_bot = 600;
ages = ([Phan.b_age] + [Phan.t_age])/2;

%Binning cumulative volumes for stacked bar chart
time_binwidth = 10; %[Myr]
range = (t_top+time_binwidth/2):time_binwidth:(t_bot-time_binwidth/2);

%Times for decarbonation calculation
t_range(2) = t_top;
t_range(1) = t_bot;

%% Compute volumes of carbonate and siliciclastic

Col_areas = [Phan.col_area]'; %[km^2]
Min_thick = [Phan.min_thick]'/1000; %[km]
Max_thick = [Phan.max_thick]'/1000; %[km]

%All in units are [km^3]
%col 1: minimum volume
%col 2: maximum volume
unit_vols = [Min_thick.*Col_areas Max_thick.*Col_areas];
%volumes corrected for erosion, following compensational model
tau = 400;%Inf; %[Myr] half-life for erosion
erosion_correction = [2.^(ages'/tau) 2.^(ages'/tau)];

unit_vols = unit_vols .* erosion_correction;

%% Load areal addition rates + interpolate

%Loading  surface area addition rates as reported in Cao et al., 2017 EPSL
%opts: Global uncorrected for erosion, N. America uncorrected for erosion
surf_area_file = filename;
surf_area_mat = load(surf_area_file);
surf_area_times = surf_area_mat(:,1); %[Mya]
surf_area_areas = surf_area_mat(:,2); %[1000 km^2/Myr]
%correct for erosion in accordance with the Macrostrat compilation
surf_area_areas = surf_area_areas.*(2.^(surf_area_times/tau));
%isolating surface area addition rates at the specified time intervals
surf_area_times = surf_area_times(...
    surf_area_times <= max(range) & surf_area_times >= min(range));
surf_area_areas = surf_area_areas(...
    surf_area_times <= max(range) & surf_area_times >= min(range));

%interpolating the values and replacing NaNs with zeros
areal_addition = interp1(surf_area_times,surf_area_areas,range)'; 
                 %[10^3 km^2/Myr] surface area addition rate of magma
areal_addition(isnan(areal_addition)) = 0;

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

depth  = 10; %[km] depth range over which decarbonation occurs

%Imposed burial rates
burial_rte = [.01 1];%[.033 .33]; %[km/Myr] range of burial rates for sediments
%min rate: akin to a passive margin
%max rate: akin to a foreland basin

%% Load empirical intrusion volume-aureole volume estimator

load intrusion_aureole_volume

%% Monte Carlo computation

n_iterations = 10000; %number of iterations for Monte Carlo simulation

[sed_vols, f_intr, f_aur, F_CO2_assim, F_CO2_meta] = ...
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
                                                           f_intrusion,...
                                                           f_aureole);
                                                       
%% Compute means + std

%volumes
mean_carb_vol = mean(sed_vols.carb,2); sd_carb_vol = std(sed_vols.carb,0,2);
mean_silc_vol = mean(sed_vols.silc,2); sd_silc_vol = std(sed_vols.silc,0,2);
mean_mixd_vol = mean(sed_vols.mixd,2); sd_mixd_vol = std(sed_vols.mixd,0,2);
mean_totl_vol = mean(sed_vols.totl,2); sd_totl_vol = std(sed_vols.totl,0,2);

mean_rest_vol = mean_totl_vol - ...
                (mean_carb_vol + mean_silc_vol + mean_mixd_vol);

%volume fractions
mean_f_intr = mean(f_intr,2); sd_f_intr = std(f_intr,0,2);
mean_f_aur = mean(f_aur,2); sd_f_aur = std(f_aur,0,2);

mean_f_intr(areal_addition == 0) = NaN; sd_f_intr(areal_addition == 0) = NaN;
mean_f_aur(areal_addition == 0) = NaN; sd_f_aur(areal_addition == 0) = NaN;

%decarbonation fluxes
mean_F_CO2_assim = mean(F_CO2_assim,2); sd_F_CO2_assim = std(F_CO2_assim,0,2);
mean_F_CO2_meta = mean(F_CO2_meta,2); sd_F_CO2_meta = std(F_CO2_meta,0,2);

mean_F_CO2_assim(areal_addition == 0) = NaN; sd_F_CO2_assim(areal_addition == 0) = NaN;
mean_F_CO2_meta(areal_addition == 0) = NaN; sd_F_CO2_meta(areal_addition == 0) = NaN;

%% Visualization: decarb rate

%% left side of figure: volumes
figure

%colormaps
colors = colormap(pink(7));

subplot(4,2,1)
b = bar([range' range' range' range'],...
        [mean_carb_vol mean_silc_vol mean_mixd_vol mean_rest_vol],...
        'Stacked');
for i = 1:length(b)
    b(i).FaceColor = colors(i,:);
end

grid on
legend('Carbonate','Siliciclastic','Mixed','Other','Interpreter','Latex')
ylabel('Volume [km$^3$]','Interpreter','Latex')
title('Rock and granitoid volumes','Interpreter','Latex','FontSize',20)
set(gca,'TickLabelInterpreter','Latex','XTickLabel',[])

subplot(4,2,3)
bar(range,areal_addition,'FaceColor',[1 1 1],...
    'FaceAlpha',1,'EdgeColor','k')
grid on
ylabel('F$_{\mathrm{granitoid}}$ [1000 $\times$ km$^2$/Myr]','Interpreter','Latex')
ylim([0 15])
set(gca,'TickLabelInterpreter','Latex','XTickLabel',[])

subplot(4,2,[5 7])
plot(range,cumsum(mean_carb_vol,'reverse'),'-','Color',colors(1,:)); hold on;
plot(range,cumsum(mean_silc_vol,'reverse'),'-','Color',colors(2,:))
plot(range,cumsum(mean_mixd_vol,'reverse'),'-','Color',colors(3,:))
plot(range,cumsum(mean_rest_vol,'reverse'),'-','Color',colors(4,:))
plot(range,cumsum(1000*time_binwidth*areal_addition*depth,'reverse'),'k-.')
grid on
xlabel('Time [Mya]','Interpreter','Latex')
ylabel('Cumulative volume [km$^3$]','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex')
legend('Carbonate','Siliciclastic','Mixed','Other',filestring,'Interpreter','Latex')

%% right side of figure: fluxes
subplot(4,2,2)
errorbar(range,mean_f_intr,sd_f_intr,'ko-','MarkerSize',8,'MarkerFaceColor',colors(2,:))
hold on
errorbar(range,mean_f_aur,sd_f_aur,'kd-','MarkerSize',8,'MarkerFaceColor',colors(4,:))
grid on
ylim([0 1])
ylabel('V$_{i}$/V$_{\mathrm{sediment}}$ [km$^3$/km$^3$]',...
       'Interpreter','Latex')
title('Decarbonation predictions','Interpreter','Latex','FontSize',20)
legend('Granitoid: Global','Aureole: Global','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex','XTickLabel',[])

subplot(4,2,[4 6])
patch_fix = ~isnan(mean_F_CO2_meta);
%assimilation flux
%errorbar(range,mean_F_CO2_assim,sd_F_CO2_assim,'ko-','MarkerSize',8,'MarkerFaceColor','k')
patch([range(patch_fix)'; fliplr(range(patch_fix))'],...
      [mean_F_CO2_assim(patch_fix) + sd_F_CO2_assim(patch_fix); ...
      flipud(mean_F_CO2_assim(patch_fix) - sd_F_CO2_assim(patch_fix))],...
      colors(2,:),'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
hold on
plot(range,mean_F_CO2_assim,'ko-','MarkerSize',8,'MarkerFaceColor',colors(2,:))

%aureole flux
patch([range(patch_fix)'; fliplr(range(patch_fix))'],...
      [mean_F_CO2_meta(patch_fix) + sd_F_CO2_meta(patch_fix); ...
      flipud(mean_F_CO2_meta(patch_fix) - sd_F_CO2_meta(patch_fix))],...
      colors(4,:),'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
hold on
plot(range,mean_F_CO2_meta,'kd-','MarkerSize',8,'MarkerFaceColor',colors(4,:))

plot(range,mean_F_CO2_assim + mean_F_CO2_meta,'k--','LineWidth',3)
grid on
ylabel('F$_{\mathrm{CO_2}}$ [Mt/yr]','Interpreter','Latex')
legend('F$_{\mathrm{assimilation}}$','F$_{\mathrm{aureole}}$','Mean F$_{\mathrm{CO_2}}$: Global','Location','Northwest','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex','XTickLabel',[])
ylim([0 200])
xlim([0 600])

subplot(4,2,8)
load('Royer2004pCO2.txt')
load('Foster2017pCO2.txt');
smooth_Royer_time = linspace(min(-Royer2004pCO2(:,1)),...
                              max(-Royer2004pCO2(:,1)),...
                              200);
smooth_Royer_pCO2 = interp1(-Royer2004pCO2(:,1),300*Royer2004pCO2(:,2),...
                             smooth_Royer_time);
h = plot(smooth_Royer_time(smooth_Royer_time>420),...
     smooth_Royer_pCO2(smooth_Royer_time>420),'k--','LineWidth',2);
hold on
patch([Foster2017pCO2(:,1); flipud(Foster2017pCO2(:,1))],...
      [Foster2017pCO2(:,3); flipud(Foster2017pCO2(:,4))],colors(3,:),...
      'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
plot(Foster2017pCO2(:,1),Foster2017pCO2(:,2),'k-','LineWidth',2)

ylabel('Atm. pCO$_2$ [ppm]','Interpreter','Latex')
xlabel('Time [Ma]','Interpreter','Latex')
legend('Royer et al., 2004','Foster et al., 2017',...
       'Interpreter','Latex','Location','Northwest')
set(gca,'TickLabelInterpreter','Latex')
grid on
ylim([0 6000])
xlim([0 600])