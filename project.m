clc
clear all;
close all;


%{
%%Part 1


% common model fiber parameters   
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
Fs = 100e3;   % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 200e-3;  % stimulus duration in seconds
rt = 10e-3;   % rise/fall time in seconds

% PSTH parameters
nrep = 10;                % number of stimulus repetitions


t = 0:1/Fs:T-1/Fs;       % time vector
F0 = 125*2.^[0:1/8:7];   % stimulus frequency vector in Hz
stimdb = -10:10:80;      % stimulus intensity in dB SPL

mxpts = length(t);
irpts = rt*Fs;

%For ANF with BF = 500Hz
BF    = 500;          % BF in Hz;
readings500=zeros(length(stimdb), length(F0));

i=1;
for intensity = stimdb
    j=1;
    for frequency = F0
        %setting up stimulus
        pin = sqrt(2)*20e-6*10^(intensity/20)*sin(2*pi*frequency*t); % unramped stimulus
        pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;    % added onset ramp
        pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts; %added offset ramp

        %running ANF model
        vihc = catmodel_IHC(pin,BF,nrep,1/Fs,T*2,cohc,cihc); 
        [synout,psth] = catmodel_Synapse(vihc,BF,nrep,1/Fs,fiberType,implnt); 
        
        %finding rate from psth
        temp = sum(reshape(psth,length(psth)/2,2))/nrep; %dividing psth into half as reading from
        spike = temp(1);                                 %ANF for 2*T time
       
        
        readings500(i,j)=spike;
        j=j+1;
    end
    i=i+1;
end

figure(1)
title('Rate vs frequency for ANF BF=500Hz, figure 1')
hold on
for i = 1:length(stimdb)
    plot(F0, readings500(i,:), 'DisplayName', num2str(stimdb(i),'%d dB'))
    xscale("log")
    legend
end
hold off


%For ANF with BF = 4000Hz
BF    = 4000;          % BF in Hz;
readings4000=zeros(length(stimdb), length(F0));

i=1;
for intensity = stimdb
    j=1;
    for frequency = F0
        %setting up stimulus
        pin = sqrt(2)*20e-6*10^(intensity/20)*sin(2*pi*frequency*t); % unramped stimulus
        pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;    % added onset ramp
        pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts; %added offset ramp

        %running ANF model
        vihc = catmodel_IHC(pin,BF,nrep,1/Fs,T*2,cohc,cihc); 
        [synout,psth] = catmodel_Synapse(vihc,BF,nrep,1/Fs,fiberType,implnt); 
        
        %finding rate from psth
        temp = sum(reshape(psth,length(psth)/2,2))/nrep; %dividing psth into half as reading from
        spike = temp(1);                                 %ANF for 2*T time
       
        
        readings4000(i,j)=spike;
        j=j+1;
    end
    i=i+1;
end

figure(2)
title('Rate vs frequency for ANF BF=4000Hz, figure 2')
hold on
for i = 1:length(stimdb)
    plot(F0, readings4000(i,:), 'DisplayName', num2str(stimdb(i),'%d dB'))
    xscale("log")
    legend
end
hold off

%For Rate v.s. Intensity
figure(3)

hold on
plot(stimdb, readings500(:,17), 'DisplayName', 'BF=500Hz')
plot(stimdb, readings4000(:,41), 'DisplayName', 'BF=4000Hz')
legend
hold off
title('Rate vs Intensity for ANFs BF=500Hz, 4000Hz, figure 3')

%}


%{
%%Part 2


cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
Fs = 100e3;   % sampling rate in Hz (must be 100, 200 or 500 kHzs
rt = 10e-3;   % rise/fall time in seconds

[y,fs] = audioread('fivewo.wav', [140000,150000]);
y = resample(y, Fs, fs);

y=transpose(y);

dBSPL=20*log10((sqrt(sum(y.^2)/length(y)))/(20*10^(-6)));

y_normal=y/(10^((dBSPL/20)+1)); %steady state stimulus renormalised to -20dBSPL

Y=zeros(length(0:1/4:5),length(y));
for i = 0:1/4:5
Y(4*(i)+1,:)=y_normal*20*10^(i);
end

%For ANF with BF = 600Hz

BF    = 600;          % BF in Hz;
readings600=zeros(1,length(0:1/4:5));

T=length(y)/Fs;
mxpts = length(y);
irpts = rt*Fs;
nrep = 80;


for i = 1:length(Y(:,1))
    
    signal = Y(i,:);
    %setting up stimulus
    signal(1:irpts)=signal(1:irpts).*(0:(irpts-1))/irpts;    % added onset ramp
    signal((mxpts-irpts):mxpts)=signal((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts; %added offset ramp

    %running ANF model
    vihc = catmodel_IHC(signal,BF,nrep,1/Fs,T*2,cohc,cihc);
    [synout,psth] = catmodel_Synapse(vihc,BF,nrep,1/Fs,fiberType,implnt);

    %finding rate from psth
    temp = sum(reshape(psth,length(psth)/2,2))/nrep; %dividing psth into half as reading from
    spike = temp(1);                                 %ANF for 2*T time


    readings600(i)=spike;
end


figure(4)
plot(-20:5:80, readings600)
title('Rate vs intensity for ANF BF=600Hz')



[stimulus,fs] = audioread('fivewo.wav');
stimulus = resample(stimulus, Fs, fs);
stimulus = stimulus(1:end-1);
[sgram,f,t] = spectrogram(stimulus, hanning(2560), 1280);

figure(5)
spectrogram(stimulus, hanning(2560), 1280, 'yaxis')
title('figure 5')

stimulus = transpose(stimulus);

dBSPL = 20*log10((sqrt(sum(stimulus.^2)/length(stimulus)))/(20*10^(-6)));
stimulus_normal = stimulus/(10^((dBSPL/20)+1)); %steady state stimulus renormalised to -20dBSPL
Stimulus = zeros(3,length(stimulus));

Stimulus(1,:) = stimulus_normal*20*10^(4);
Stimulus(2,:) = stimulus_normal*20*10^(8);
Stimulus(3,:) = stimulus_normal*20*10^(16);


BF=125*2.^[0:1/8:6];
T=length(stimulus)/Fs;
mxpts = length(stimulus);
irpts = rt*Fs;
nrep=80;

readings=zeros(3,length(BF), length(stimulus));

for i = 1:3
    for j = 1:length(BF)
        
        signal = Stimulus(i,:);
        % setting up stimulus
        % signal(1:irpts)=signal(1:irpts).*(0:(irpts-1))/irpts;                        % added onset ramp
        % signal((mxpts-irpts):mxpts)=signal((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts; %added offset ramp
    
        %running ANF model
        vihc = catmodel_IHC(signal,BF(j),nrep,1/Fs,T*2,cohc,cihc);
        [synout,psth] = catmodel_Synapse(vihc,BF(j),nrep,1/Fs,fiberType,implnt);
    
        %finding rate from psth
        temp = reshape(psth,2,length(psth)/2)/nrep;                                    %dividing psth into half as reading from
        spike = temp(1,:);                                                             %ANF for 2*T time
        readings(i,j,:)=spike;
    end
end



psthbinwidths = [4e-3, 8e-3, 16e-3, 32e-3, 64e-3, 128e-3];
psthbins = round(psthbinwidths*Fs);   % number of psth bins per psth bin
PsthA = zeros(3, length(BF), 2*length(readings(1,1,:))/psthbins(1));
for i = 1:3
    for j = 1:length(BF)
        %for psthbinwidth = psthbinwidths
            timeout = (1:length(readings(i,j,:)));
            psth=readings(i,j,:);
            psthshift=zeros(1,1,length(psth));
            for k = 1:length(psth)-(psthbins(1)/2)
                psthshift(k) = psth(psthbins(1)/2+k);
            end
            tmp=cat(1, psth(1:length(psthshift)), psthshift);
            tmp=sum(tmp);
            pr = sum(reshape(tmp,psthbins(1)/2,2*length(tmp)/psthbins(1))); % pr of spike in each bin
            PsthA(i,j,:) = pr/psthbinwidths(1); % psth in units of spikes/s
        %end
    end
end

figure(61)
for i = 1:3
    subplot(3,1,i)
    imagesc(squeeze(PsthA(i,:,:)))
    xticks(1:32:2*length(readings(1,1,:))/psthbins(1))
    xticklabels(0:16*psthbinwidths(1):T)
end
title('figure 6A')


PsthB = zeros(3, length(BF), 2*length(readings(1,1,:))/psthbins(2));
for i = 1:3
    for j = 1:length(BF)
        %for psthbinwidth = psthbinwidths
            timeout = (1:length(readings(i,j,:)));
            psth=readings(i,j,:);
            psthshift=zeros(1,1,length(psth));
            for k = 1:length(psth)-(psthbins(2)/2)
                psthshift(k) = psth(psthbins(2)/2+k);
            end
            tmp=cat(1, psth(1:length(psthshift)), psthshift);
            tmp=sum(tmp);
            pr = sum(reshape(tmp,psthbins(2)/2,2*length(tmp)/psthbins(2))); % pr of spike in each bin
            PsthB(i,j,:) = pr/psthbinwidths(2); % psth in units of spikes/s
        %end
    end
end

figure(62)
for i = 1:3
    subplot(3,1,i)
    imagesc(squeeze(PsthB(i,:,:)))
    xticks(1:32:2*length(readings(1,1,:))/psthbins(2))
    xticklabels(0:16*psthbinwidths(2):T)
end
title('figure 6B')


PsthC = zeros(3, length(BF), 2*length(readings(1,1,:))/psthbins(3));
for i = 1:3
    for j = 1:length(BF)
        %for psthbinwidth = psthbinwidths
            timeout = (1:length(readings(i,j,:)));
            psth=readings(i,j,:);
            psthshift=zeros(1,1,length(psth));
            for k = 1:length(psth)-(psthbins(3)/2)
                psthshift(k) = psth(psthbins(3)/2+k);
            end
            tmp=cat(1, psth(1:length(psthshift)), psthshift);
            tmp=sum(tmp);
            pr = sum(reshape(tmp,psthbins(3)/2,2*length(tmp)/psthbins(3))); % pr of spike in each bin
            PsthC(i,j,:) = pr/psthbinwidths(3); % psth in units of spikes/s
        %end
    end
end

figure(63)
for i = 1:3
    subplot(3,1,i)
    imagesc(squeeze(PsthC(i,:,:)))
    xticks(1:32:2*length(readings(1,1,:))/psthbins(3))
    xticklabels(0:16*psthbinwidths(3):T)
end
title('figure 6C')


PsthD = zeros(3, length(BF), 2*length(readings(1,1,:))/psthbins(4));
for i = 1:3
    for j = 1:length(BF)
        %for psthbinwidth = psthbinwidths
            timeout = (1:length(readings(i,j,:)));
            psth=readings(i,j,:);
            psthshift=zeros(1,1,length(psth));
            for k = 1:length(psth)-(psthbins(4)/2)
                psthshift(k) = psth(psthbins(4)/2+k);
            end
            tmp=cat(1, psth(1:length(psthshift)), psthshift);
            tmp=sum(tmp);
            pr = sum(reshape(tmp,psthbins(4)/2,2*length(tmp)/psthbins(4))); % pr of spike in each bin
            PsthD(i,j,:) = pr/psthbinwidths(4); % psth in units of spikes/s
        %end
    end
end

figure(64)
for i = 1:3
    subplot(3,1,i)
    imagesc(squeeze(PsthD(i,:,:)))
    xticks(1:32:2*length(readings(1,1,:))/psthbins(4))
    xticklabels(0:16*psthbinwidths(4):T)
end
title('figure 6D')



PsthE = zeros(3, length(BF), 2*length(readings(1,1,:))/psthbins(5));
for i = 1:3
    for j = 1:length(BF)
        %for psthbinwidth = psthbinwidths
            timeout = (1:length(readings(i,j,:)));
            psth=readings(i,j,:);
            psthshift=zeros(1,1,length(psth));
            for k = 1:length(psth)-(psthbins(5)/2)
                psthshift(k) = psth(psthbins(5)/2+k);
            end
            tmp=cat(1, psth(1:length(psthshift)), psthshift);
            tmp=sum(tmp);
            pr = sum(reshape(tmp,psthbins(5)/2,2*length(tmp)/psthbins(5))); % pr of spike in each bin
            PsthE(i,j,:) = pr/psthbinwidths(5); % psth in units of spikes/s
        %end
    end
end

figure(65)
for i = 1:3
    subplot(3,1,i)
    imagesc(squeeze(PsthE(i,:,:)))
    xticks(1:32:2*length(readings(1,1,:))/psthbins(5))
    xticklabels(0:16*psthbinwidths(5):T)
end
title('figure 6E')

PsthF = zeros(3, length(BF), 2*length(readings(1,1,:))/psthbins(6));
for i = 1:3
    for j = 1:length(BF)
        %for psthbinwidth = psthbinwidths
            timeout = (1:length(readings(i,j,:)));
            psth=readings(i,j,:);
            psthshift=zeros(1,1,length(psth));
            for k = 1:length(psth)-(psthbins(6)/2)
                psthshift(k) = psth(psthbins(6)/2+k);
            end
            tmp=cat(1, psth(1:length(psthshift)), psthshift);
            tmp=sum(tmp);
            pr = sum(reshape(tmp,psthbins(6)/2,2*length(tmp)/psthbins(6))); % pr of spike in each bin
            PsthF(i,j,:) = pr/psthbinwidths(6); % psth in units of spikes/s
        %end
    end
end

figure(66)
for i = 1:3
    subplot(3,1,i)
    imagesc(squeeze(PsthF(i,:,:)))
    xticks(1:32:2*length(readings(1,1,:))/psthbins(6))
    xticklabels(0:16*psthbinwidths(6):T)
end
title('figure 6F')
%}



%%Part 3


cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

Fs = 100e3;   % sampling rate in Hz (must be 100, 200 or 500 kHz)


[stimulus,fs] = audioread('fivewo.wav');
stimulus = resample(stimulus, Fs, fs);
stimulus = stimulus(1:end-1);


signal = transpose(stimulus);

psthbinwidth=0.1e-3;
psthbin=round(psthbinwidth*Fs);


BF=125*2.^[0:1/8:6];
T=length(stimulus)/Fs;
nrep=80;

for j = 1:length(BF)

    %running ANF model
    vihc = catmodel_IHC(signal,BF(j),nrep,1/Fs,T*2,cohc,cihc);
    [synout,psth] = catmodel_Synapse(vihc,BF(j),nrep,1/Fs,fiberType,implnt);

    %finding rate from psth
    temp = reshape(psth,2,length(psth)/2)/nrep;                                  %dividing psth into half as reading from
    spike(j,:) = temp(1,:);                                                           %ANF for 2*T time
end

PSTH = zeros(length(BF), 2*length(spike(1,:))/psthbin);

for j = 1:length(BF)
    %for psthbinwidth = psthbinwidths
    timeout = (1:length(spike(1,:)));
    psth=spike(j,:);
    psthshift=zeros(1,length(psth));
    for k = 1:length(psth)-(psthbin/2)
        psthshift(k) = psth(psthbin/2+k);
    end
    tmp=cat(1, psth(1:length(psthshift)), psthshift);
    tmp=sum(tmp);
    pr = sum(reshape(tmp,psthbin/2,2*length(tmp)/psthbin));                     % pr of spike in each bin
    PSTH(j,:) = pr/psthbinwidth;                                                % psth in units of spikes/s
    %end
end

windowSize = 1280; 
overlap = windowSize / 2; 
nfft = windowSize; 
maxValues = zeros(length(BF),2*length(PSTH)/windowSize-1); 
frequencies = zeros(length(BF),2*length(PSTH)/windowSize-1); 

for i = 1:length(BF)
        maxValue=[];
        frequency=[];
    for start = 1:(windowSize - overlap):(length(PSTH(1,:)) - windowSize + 1)
        window = PSTH(i,start:start + windowSize - 1);
        fftWindow = fft(window, nfft);
        fftMagnitude = abs(fftWindow / windowSize);
        fftMagnitude = fftMagnitude(1:windowSize/2+1);
        fftMagnitude(1) = 0;                                               %removing DC value
        
        [sortedValues, sortedIndices] = sort(fftMagnitude, 'descend');
        
        Value = sortedValues(1); % Maximum FFT value
        Index = sortedIndices(1);
        secondMaxValue = sortedValues(2); % Second maximum FFT value

        threshold = 1.3 * secondMaxValue;

        if Value > threshold
            freq = (Index - 1) * Fs / nfft;
            maxValue = [maxValue, Value];
            frequency = [frequency, freq];
        else
            maxValue = [maxValue, 0];
            frequency = [frequency, 0];
        end
    end
    maxValues(i,:)=maxValue;
    frequencies(i,:)=frequency;
end



[S, F, T, P] = spectrogram(signal, windowSize, overlap, nfft, Fs, 'yaxis');


figure(7);
surf(T,F,10*log10(P),'EdgeColor','none');
axis tight;
view(0,90);
colorbar;
title('Spectrogram with Max Frequency Indicators, figure 7');
xlabel('Time (s)');
ylabel('Frequency (Hz)');


hold on;
for i = 1:8:41
   
    frequency = frequencies(i,:);
    
    timeVector = linspace(min(T), max(T), length(frequency));
    
    nonZeroIndices = find(frequency); 
    plot(timeVector(nonZeroIndices), frequency(nonZeroIndices), '*', 'MarkerSize', 10);
end

hold off;



figure(8);
surf(T,F,10*log10(P),'EdgeColor','none');
axis tight;
view(0,90);
colorbar;
title('Spectrogram with Max Frequency Indicators, figure 8');
xlabel('Time (s)');
ylabel('Frequency (Hz)');


hold on;
for i = 5:8:45
   
    frequency = frequencies(i,:);
    
    timeVector = linspace(min(T), max(T), length(frequency));
    
    nonZeroIndices = find(frequency); 
    plot(timeVector(nonZeroIndices), frequency(nonZeroIndices), '*', 'MarkerSize', 10);
end

hold off;