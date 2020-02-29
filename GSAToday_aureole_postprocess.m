%file: GSAToday_aureole_postprocess.m

%author:  Evan J. Ramos
%date:    20 Aug 2019

%Description: This file post processes results from 
%GSAToday_aureole_prediction.m so that a relationship between intrusion 
%size and aureole size (over a range of crustal permeabilities) can be 
%determined.

close all
clear
clc

set(0, 'DefaultAxesFontSize', 25,...
       'DefaultLineLinewidth', 2,...
       'DefaultAxesXColor','k',...
       'DefaultAxesYColor','k')

%% Select file

opt = 1;
if opt == 1
    %Each file contains these variables: Grid, Fracs, areas, ks
    filename = 'GSA_Today_2019_FINAL';
    load(filename)
elseif opt == 2
    %Run the script that performs the numerical solutions
    GSAToday_aureole_prediction
end

%% Declaring parameters

%Domain parameters
Area   = Grid.Lx*Grid.Ly; %[m^2] domain area

%Aureole parameters
frac_min = .2; %[m] aureole size limited at 20% decarbonated (arbitrary)

%% Computing aureole areas

%Constant area
aureole_areas = zeros(size(Fracs.frac));
for i = 1:numel(Fracs.frac)
    
    %isolate cells around pluton
    aureole = Grid.dof(~ismember(Grid.dof,Fracs.plut{i}.dof));
    frac = Fracs.frac{i}(aureole,end);
    
    %compute area by multiplying the area of an individual grid cell by the
    %number of cells that meet the aureole criterion
    aureole_areas(i) = Grid.dx*Grid.dy*sum(frac);
end
aur_avg = mean(aureole_areas,2);
aur_std = std(aureole_areas,0,2);

%% Plotting

figure
colors = colormap(pink(5));

p = plot(areas/Area,aureole_areas/Area,'-');
for i = 1:size(aureole_areas,2)
    p(i).Color = colors(i,:);
end
hold on
errorbar(areas/Area,aur_avg/Area,aur_std/Area,'kd-.','MarkerFaceColor','k',...
    'MarkerSize',9,'LineWidth',2)

grid on
xlabel('V$_{\mathrm{intrusion}}$/(V$_{\mathrm{host\:rock}}$ + V$_{\mathrm{intrusion}}$)',...
       'interpreter','latex')
ylabel('V$_{\mathrm{aureole}}$/(V$_{\mathrm{host\:rock}}$ + V$_{\mathrm{intrusion}}$)',...
       'interpreter','latex')
set(gca,'TickLabelInterpreter','Latex')

l = legend('-15','-16','-17','Mean','interpreter','latex','location','northeast');
title(l,'log$_{10}$ permeability [m$^2$]','interpreter','latex')
title('Predicted aureole volume from intrusion volume','interpreter','latex')