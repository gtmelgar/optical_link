rng('shuffle') 
OpAttn_dBm = -8:1:-2;
numOfBits = 1e5; 

plotData = true;
numData = length(OpAttn_dBm);
%% Received Power
avePowerTX = 1-3;      % Watts, this doesn't matter as long as you're looking for RX power
q = 1.6021766208e-19;
kB = 1.38064852e-23;
%% Parameters: 
ER_dB = 4;            % dB
noiseTIA = 3.0e-6;        % Amps
Rl = 50;                % Ohms
darkCurrent = 3e-9;     % Amps
RIN = -145;             % dB/Hz
bitRate = 25e9;         % in Hz
laserBW = 18e9;         % in Hz
rinBW = 18e9;           % in Hz
resp = 0.35;             % Responsivity in A/W
pinBW = 30e9;           % BW of photoDetector in Hz
tiaBW = 30e9;           % BW of TIA in Hz
distance = 100;         % fiber distance in meters
       
samplesPerBit = 16;
tempCel = 20;
dcBlockBW = 40e3;
% Ratio of power in each mode must add to 1 or less, must be a 1xn matrix
modePowerDist = [0.1 0.1 0.1 0.1];
modePowerDist = modePowerDist/sum(modePowerDist);
% Wavelengths in the signal, must be a 1xn matrix
modeLambda = [855.9 855.1 854.8 854.22 854.15 853.9];
% Used for prbs pattern generator
% PRBS7:    g = [7 6]    init = [1 1 0 0 0 0 1]
% PRBS15:   g = [15 14]  init = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 1]

init = [1 1 0 0 0 0 1];
g = [7 1 0];

avePowerdBm = 10*log10(1000*avePowerTX); % convert watts to dBm
tempKelvin = tempCel + 273.15;           % covenrt degC to kelvin
sampleRate = bitRate*samplesPerBit;      % total BW of system
ER_LIN = 10^(ER_dB/10);                  % convert extinction ratio dB to linear

% calculate power of bit 1 and bit 0 given extinction ratio and average power
P0 = 2 * avePowerTX*(1/(ER_LIN + 1));
P1 = 2 * avePowerTX*(ER_LIN/(ER_LIN + 1));

%% Create PRBS signal of logical ones and zeros
% create the PRBS Pattern using initial condition and polynomial
% init must have a vector equal in size to the polynomial degree "g"
% PRBS7:    g = [7 6]    init = [1 1 0 0 0 0 1]
% PRBS15:   g = [15 14]  init = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 1]

prbsPattern = rand(2^11,1); % use rand for now since comm toolbox is needed

%% Upsampling and Pulse Shaping
% create pulse

gaussfilter = gaussdesign(laserBW/bitRate,length(prbsPattern),samplesPerBit); % change arg1
h = fir1(100,(laserBW)/(25e9*16/2));
[signalTX, prbsSignalLogical,~,testStream] = scaleRepeatPAM4(prbsPattern,samplesPerBit,...
    numOfBits,1,P0,P1);

[b,a] = cheby1(5,.1,(laserBW)/(bitRate*16/2));
signalTX=filter(b,a,signalTX);

numSamples = length(signalTX);
% timeVector = (1:numOfBits*samplesPerBit)'*1/sampleRate; % If wanting to plot vs time(s)
%% Add RIN to Tx Pulse (rememeber to shuffle randn)
% RIN filtered by Bessel function 
noiseRIN = sqrt(sampleRate/2 * 10^(RIN/10)).*randn(numSamples,1).*signalTX;

% Create bessel filter, to be later implemented with RC filter
% must be multiplied by 1.514 for -3dB BW to match 
[num,den] = besself(4,2*pi*rinBW*1.514);
[numz, denz] = impinvar(num,den,sampleRate);

%[numz, denz] = butter(2,rinBW/(sampleRate/2));
noiseRIN = filter(numz,denz,noiseRIN);

%% Transmitted Signal
signalTX = signalTX + noiseRIN;

numErrors = zeros(numData,1);
ratioSER = zeros(numData,1);
threshold = zeros(numData,3);
meanRxPower = zeros(numData,1);
for n = 1:numData
    
    %% Attenuation
    signalTemp = signalTX .* 10^((-avePowerdBm+OpAttn_dBm(n))/10);
    
    %% Fiber Model
    % This is the fiber model comment it out if you don't want back to back
%     signalTemp = fiberDelay(signalTemp,modePowerDist,modeLambda,distance,sampleRate);
    
    %% Received Signal
    % convert power signal into electrical signal
    meanRxPower(n) = mean(signalTemp);
    signalRX = signalTemp * resp;

    %% Add PIN noise and Filter signal
    % Shot Noise and Thermal Noise
    shotNoise = sqrt(2.*q.*(signalRX + darkCurrent)*sampleRate/2).*randn(numSamples,1);    
    thermalNoise = sqrt(4*kB*tempKelvin/Rl*sampleRate/2).*randn(numSamples,1);
    % Add noise 
    signalRX = signalRX + shotNoise + thermalNoise;

    % create a bessel filter as a place holder - will prob stay the same
    [num,den] = besself(4,2*pi*pinBW*1.514);
    [numz,denz] = impinvar(num,den,sampleRate);
    signalRX = filter(numz,denz,signalRX);
    
     %% Add TIA Noise
    tiaNoise = noiseTIA .* randn(numSamples,1);
    signalRX = signalRX + tiaNoise;
    % create a bessel filter as a place holder - will prob stay the same
    [num,den] = besself(4,2*pi*tiaBW*1.514);
    [numz,denz] = impinvar(num,den,sampleRate);
    signalRX = filter(numz,denz,signalRX);
    
    %% Scaling
    signalRX = signalRX-mean(signalRX);
    signalRX = sqrt(5)*signalRX/std(signalRX);
    testStream = testStream-mean(testStream);
    testStream = sqrt(5)*testStream/std(testStream);
    
    %% Cross Correlation in order to find delay
    [acor,lag] = xcorr(testStream,signalRX(1:50*samplesPerBit));
    [~,delay] = max(abs(acor));
    signalRX = signalRX(-lag(delay)+1:end);
    
    %% Add DC Block
    p = dcblock(dcBlockBW/(sampleRate/2));
    signalRX = filter([1,-1],[1 -p],signalRX);
    %if (lag(delay) < (1000*16) || lag(delay) > 0)
    %         lag(delay) = -floor(samplesPerBit/2);
    %end
    
    delayOfLag = lag(delay) %#ok<NOPTS>
    
    %% DownSample Signal and Test BER
    downsampleSignal = downsample(signalRX,samplesPerBit,8);
    % threshold(n) = mean(signalRX(62500*16:end));
    %62500*16

    downSampleLength = length(downsampleSignal);
    prbsSignalLength = length(prbsSignalLogical);

    if downSampleLength < prbsSignalLength
        prbsSignalLogical = prbsSignalLogical(1:downSampleLength);
    elseif prbsSignalLength < downSampleLength
        downsampleSignal = downsampleSignal(1:prbsSignalLength);
    end
    
    %% Compute Optimal Thershold
    threshold(n,:) = optThersholdPAM4(downsampleSignal(100000:end), prbsSignalLogical(100000:end));
    % threshold(n) = mean(downsampleSignal(200:end));

    ReceivedBits = zeros(downSampleLength,1);
    
    ReceivedBits(downsampleSignal <= -2.33) = -3;
    ReceivedBits(downsampleSignal > -2.33 & downsampleSignal <= -0.0064) = -1;
    ReceivedBits(downsampleSignal > -0.0064 & downsampleSignal <= 2.24) = 1;
    ReceivedBits(downsampleSignal > 2.24) = 3;

    % Get BER
    [numErrors(n),ratioSER(n)] = symerr(prbsSignalLogical(1:end),ReceivedBits(1:end))
    
end
meanPowerdBm = 10*log10(1000*meanRxPower);
semilogy(meanPowerdBm,ratioSER,'-*');
legend('simulation', 'btb', 'left', 'right');
ax = gca;
set(ax,'YMinorTick','on');
grid on
