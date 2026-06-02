%close all;
clc;
hold on;

%% Received Power
P_rx_dBm = (-35: 1: 0);

%% Parameters: Changed Frequently
ER_dB = 8;
TIA_Noise_SD = 2e-6; % Amps
Rl = 50;             % Ohms
dark_current = 3e-9; % Amps
RIN = -150;          % dB/Hz

% option = 0 => i-decision/i-threshold is average power  
% option = 1 => i-decision/i-threshold is chosen as shown in Agrawal
option = 1;

% graphtype = 0 => log scale 
% graphtype = 1 => loglog scale
graphtype = 0;

%% ISI/MPN Parameters
LengthOfFiber = 1; % meters
BW_m = 4700; % modal BW for fiber MHz*km
lamb_c = 850; % laser center wavelength nm
lamb_0 = 1300; % zero dispersion wavelngth nm
sig_lambda = 0.3; % RMS spectral width of the laser spectrum in nm
laser_rise_time = 10; % Ts in ps
BW_receiver = 30e3; % in MHz -3dB BW of optical receiver electrical
BitRate = 30e3; % in MHz
S_0 = 0.11; % dispersion slope parameter ps/(nm^2 * km)
k_mpn = 0.1; % ratio
DCD = 0; %ps Duty-cycle-Distortion

%% Parameters: Changed Less Frequently
Resp = 0.3; % responsivity in A/W
RINNOISE_BW = 15.67e9; % noise of RIN in Hz
delta_f_noise = 19.9588e9; % effective noise BW in HZ
T_kelvin = 293.15; % temperature in K

%% SystemVue Values
P_rx_dBm_SysVue = [
-30
-27
-24
-21
-18
-15
-12
-9
-6
-3
-0
];

BER_Test = [
0.48125
0.4638
0.4289
0.35965
0.23967
0.09157
0.01357
1.05E-03
9.63E-05
2.75E-05
1.38E-05
];

%% Constants
kB = 1.3806488e-23;   % Boltzmann constant in J/K
q = 1.602176565e-19;  % Elementary charge in coulombs

%% Convert to Linear Terms
P_rx = 10.^(P_rx_dBm / 10)/1000;
ER_rat = 10^(ER_dB/10); 

%% Calulate Power of binary 1 and binary 0
P1 = 2 * P_rx .* (ER_rat / (ER_rat + 1));    
P0 = 2 * P_rx .* (1/(ER_rat+1));

%% Calculate Current Received
I1 = P1 .* Resp;
I0 = P0 .* Resp;

%% Dispersion Related Penalties

if LengthOfFiber > 0
    LengthOfFiber_km = LengthOfFiber/1000;
    D1 = 0.25 * S_0 * (lamb_c - lamb_0^4 / lamb_c^3); % ps/(nm-km) 
    D2 = 0.7 * S_0 * sig_lambda; % ps/(nm-km)
    Dispersion = sqrt(D1^2 + D2^2); % ps/(nm-km)
    BW_CD = 0.187*1e6 / (LengthOfFiber_km * sig_lambda * Dispersion); % MHz
    BW_m_eff = BW_m / LengthOfFiber_km; % MHz

    T_fiber_exit = sqrt(laser_rise_time^2 ...
        + 1e6*((480/BW_m_eff)^2 + (480/BW_CD)^2)); %ps

    T_receiver = 329e3/BW_receiver; %ps
    T_composite = sqrt(T_fiber_exit^2 + T_receiver^2); %ps
    T_eff = (1/(BitRate*10^6) - DCD*1e-12)*10^12; %pseconds
    BW_eff = 1/(T_eff); %THz = 1/ps
    
    %ISI power penalty
    T_b = 1/(BitRate)*1e6; % ps
    
    %%% Worst case eye closure used in Chalmers "4-PAM for High-Speed"
    E_m = 1.425 * exp(-1.28*(T_eff/T_composite)^2); 
    
    %%% Worst case eye closure used in IEEE
%     ht0 = 0.5*(erf(2.563 * T_eff/(2*sqrt(2)*T_composite))...
%         - erf(-2.563*T_eff/(2*sqrt(2)*T_composite))); 
%     P_ISI_IEEE_dB = 10*log10(1/(2*ht0 - 1));
    %%%
    
    if E_m > 1
        error('Eye Closed Due to ISI');
    end

    P_ISI = (1 - E_m); % worst case eye penalty Chalmers
    P_ISI_dB = 10*log10(P_ISI);
    ISI_oma = (I1-I0) * (1 - P_ISI);
    I1 = I1 - ISI_oma/2;
    I0 = I0 + ISI_oma/2;
    sigma_mpn = k_mpn/sqrt(2) * ...
        (1 - exp(-pi*BW_eff*Dispersion*LengthOfFiber_km*sig_lambda));
else 
    sigma_mpn = 0;
end

%% Noise Sources
thermalNoise_variance = 4*kB*T_kelvin*delta_f_noise/Rl;              
shot_noise1_variance = 2.*q.*(I1 + dark_current)*delta_f_noise;    
shot_noise0_variance = 2.*q.*(I0 + dark_current)*delta_f_noise;    
RIN_noise1_SD = sqrt(RINNOISE_BW * 10^(RIN/10))*I1;
RIN_noise0_SD = sqrt(RINNOISE_BW * 10^(RIN/10))*I0;
MPN_noise1_SD = sigma_mpn * I1;
MPN_noise0_SD = sigma_mpn * I0;

%% Add Up All Noise 
sd1_noise_tot = sqrt(TIA_Noise_SD.^2 + shot_noise1_variance + ...
    thermalNoise_variance + RIN_noise1_SD.^2 + MPN_noise1_SD.^2);
sd0_noise_tot = sqrt(TIA_Noise_SD.^2 + shot_noise0_variance + ...
    thermalNoise_variance + RIN_noise0_SD.^2 + MPN_noise1_SD.^2);

%% Calculate Threshold Current
if option == 0
    % Use average power
    I_thres = P_rx*Resp;
else 
    % Use Agrawal threshold
    I_thres = (I1.*sd0_noise_tot+I0.*sd1_noise_tot)./ ...
        (sd1_noise_tot+sd0_noise_tot);
end

%% Calculate BER
BER = (1/4) .* (erfc((I1-I_thres)./(sqrt(2).*sd1_noise_tot)) + ...
    erfc((I_thres-I0)./(sqrt(2).*sd0_noise_tot)));
BER = BER';

%% Plot Graphs
%hold on
if (graphtype == 0)
    % log plot
    h = plot(P_rx_dBm, -log10(BER), '--+b', ...
        P_rx_dBm_SysVue, -log10(BER_Test),...
        '-xr','LineWidth', 2.5, 'MarkerSize', 13);

    ylim([0.2, 9]);
    xlim([-35,0]);
    %set(gca,'YTick',[0,1,2,3,4,5,6,7],'FontSize', 2,'fontweight','b');

else
    % loglog plot
    h = semilogy(P_rx_dBm, -log10(BER), '--+b', ...
        P_rx_dBm_SysVue, -log10(BER_Test),...
        '-xr','LineWidth', 2.5, 'MarkerSize', 13);
    ylim([0.2, 11]);
    xlim([-35,0]);
    set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
end

grid on
title('BER Vs Rx Power','FontSize', 24,'fontweight','b');
set(gca,'Ydir','reverse')
set(h,'linewidth', 3);
set(gca,'FontSize', 21,'fontweight','b')

xlabel('P_r_x(dBm)','fontweight','b')
ylabel('-log(BER)','fontweight','b')
legend(...
    [' Analytical model:  ER=' num2str(ER_dB) 'dB, RIN=' num2str(RIN) 'dB/Hz, TIA Noise=' num2str(TIA_Noise_SD*1e6) 'uA' ],...
    ['Test Data: ER=' num2str(ER_dB)  'dB, RIN=' num2str(RIN) 'dB/Hz, TIA Noise=' num2str(TIA_Noise_SD*1e6) 'uA'],...
    'Location', 'SouthEast');