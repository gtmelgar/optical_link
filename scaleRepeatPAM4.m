function [prbsSignalFiltered,prbsSignalLogical,prbsSignal,testStream] = scaleRepeatPAM4(prbsPattern,samplesPerBit,...
    numOfBits,gaussfilter,P1,P4)

%{
    input:
        prbsPattern: pattern that will be repeated and scaled and filtered
        samplesPerBit: samples per bit
        numOfBits: number of bits that you want to create, also used to truncated signal
        gaussfilter: the shape of the filter to convolve with
        bit0, bit1: used for scaling
    
    output:
        prbsSignalFiltered: upsampled, repeated, scaled and filtered part of the signal
        prbsSignalLogical: reference signal used to test for bit errors
        prbsSignal: same as prbsSignalFiltered but not filtered
        testStream: used for cross correlation at receiver, find delay  
 %}

% repeat PRBS for n amount of time
repeatPRBS = ceil((numOfBits+1000)/length(prbsPattern));

%% Original signal used for comparison and BER measurements
prbsSignalLogical = prbsPattern(:,ones(repeatPRBS,1));
prbsSignalLogical = prbsSignalLogical(:);
%prbsSignalLogical = prbsSignalLogical(1:numOfBits);

%% Scaled Signal Used for computation
% upsample and hold pattern by the number of samples
prbsPattern = prbsPattern';
prbsPattern = prbsPattern(ones(samplesPerBit,1),:);
prbsPattern = prbsPattern(:);

OMA = P4 - P1;
P2 = P1 + OMA/3;
P3 = P4 - OMA/3;

prbsPattern(prbsPattern == -3) = P1;
prbsPattern(prbsPattern == -1) = P2;
prbsPattern(prbsPattern == 1) = P3;
prbsPattern(prbsPattern == 3) = P4;

%% Repeat upsampled PRBS pattern and Filter signal
prbsSignal = prbsPattern(:,ones(repeatPRBS,1));
% keep signal 
testStream = prbsSignal(1:samplesPerBit*127);
prbsSignal = prbsSignal(:);
prbsSignalFiltered = conv(prbsSignal,gaussfilter,'same'); 
%prbsSignalFiltered = prbsSignalFiltered(1:(numOfBits+1)*samplesPerBit);
% [num,den] = besself(4,2*pi*18e9*1.514);
% [numz, denz] = impinvar(num,den,25e9*16);
% prbsSignalFiltered = filter(numz,denz,prbsSignal);
