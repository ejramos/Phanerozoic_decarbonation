function [skarn_frac, therm_frac, net_decarb] = decarb_perc_aureole(delta_ms,...
                                                                    Ts,...
                                                                    d_styles,...
                                                                    d_bounds)

%author: Evan J. Ramos
%date:   14 Apr 2018

%Description: This function computes the percentage of decarbonation that
%the contact aureole has undergone. The user specifies the temperature
%range over which thermal decarbonation flux operates and the d18O range
%over which the skarn decarbonation flux operates. The user then specifies
%the style of decarbonation for  each (linear, kinetic, or instantaneous).
%The computation for the percentage decarbonated is based off an analytical
%solution for the steady state advection-diffusion equation.

%Inputs:
% delta_ms: array of d18O mineral values [permil VSMOW]
%       Ts: array of temperatures [K]
% d_styles: structure where each field contains strings that specify the 
%           style of decarbonation. Fields: therm, skarn.
% d_bounds: structure where each field contains the minimum and maximum
%           values over which a decarbonation rate is computed. Fields:
%           therm, skarn.

%Outputs:
% skarn_frac: array of fraction decarbonated via skarn formation [0,1]
% therm_frac: array of fraction decarbonated via thermal decomp [0,1]
% net_decarb: the sum of skarn_frac and therm_frac

%Decarbonation equation:
% fraction = (exp(Pe*x) - 1)/(exp(Pe) - 1), if Pe ~= 0
% fraction = x, if Pe = 0
% Pe: Peclet number (dimensionless number which will scale differently for
%     decarbonation style)

%% Decarbonation equations

%skarn formation
dO       = d_bounds.skarn(2) - d_bounds.skarn(1);
O_nondim = @(O) 1 - (O - d_bounds.skarn(1))/dO;
if strcmp(d_styles.skarn,'linear')
    skarn_decarb = @(O) O_nondim(O);
elseif strcmp(d_styles.skarn,'kinetic')
    PeO = 7;
    skarn_decarb = @(O) (exp(PeO*O_nondim(O)) - 1)/(exp(PeO) - 1);
elseif strcmp(d_styles.skarn,'instantaneous')
    PeO = -10;
    skarn_decarb = @(O) (exp(PeO*O_nondim(O)) - 1)/(exp(PeO) - 1);
end

%%%%% NOTE: thermal decomposition can be included as a decarbonation flux,
%%%%% but we are excluding this from our model calculations.
%thermal decomposition
dT       = d_bounds.therm(2) - d_bounds.therm(1);
T_nondim = @(T) (T - d_bounds.therm(1))/dT;
if strcmp(d_styles.therm,'linear')
    therm_decarb = @(T) T_nondim(T);
elseif strcmp(d_styles.therm,'kinetic')
    PeT = 7;
    therm_decarb = @(T) (exp(PeT*T_nondim(T)) - 1)/(exp(PeT) - 1);
elseif strcmp(d_styles.therm,'instantaneous')
    PeT = -10;
    therm_decarb = @(T) (exp(PeT*T_nondim(T)) - 1)/(exp(PeT) - 1);
end

%% Compute decarbonation fractions

skarn_frac = skarn_decarb(delta_ms);
therm_frac = therm_decarb(Ts);

%% Correct decarbonation fractions

%correct for temperatures below the minimum decarbonation temperature
therm_frac(Ts < d_bounds.therm(1)) = 0; %not decarbonated whatsoever
%take cumulative sum through each thermal decomp fraction
therm_frac = cumsum(therm_frac,2);
therm_frac(therm_frac > 1) = 1;

%correct for d18O mineral values below the minimum decarbonation d18O value
skarn_frac(delta_ms < d_bounds.skarn(1)) = 1; %totally decarbonated
%take cumulative sum through each skarn decarb fraction
skarn_frac = cumsum(skarn_frac,2);
skarn_frac(skarn_frac > 1) = 1;

%% Net fraction
net_decarb = therm_frac + skarn_frac;
net_decarb(net_decarb > 1) = 1;

end