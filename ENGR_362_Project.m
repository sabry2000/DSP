%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 362 Project
% student name: Ahmed Sabry
% student number: 99025207
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) Import audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputSignal = importdata('ENGR_362_guitar_Fs_is_48000_Hz.txt');
samplingFrequency = 48000;                             % sampling freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) Play sound of raw audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soundsc(inputSignal,samplingFrequency)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) Graph in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplingPeriod = 1/samplingFrequency;           % sampling period       
signalLength = length(inputSignal(:,1));      	% length of signal
time = (0:signalLength-1)*samplingPeriod;    	% time vector
figure,plot(time,inputSignal), axis tight
title('Time Domain Representation of Original Recording')	% labels
xlabel('t (s)')                             % labels
ylabel('y(t)')                              % labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) Graph with DFT/FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = fft(inputSignal);                               % Discrete Fourier transform
F1 = abs(Y/signalLength);                               % frequency
F2 = F1(1:signalLength/2+1);                            % half of frequency
F2(2:end-1) = 2*F2(2:end-1);                        % Discrete Fourier transform
f = samplingFrequency*(0:(signalLength/2))/signalLength;                   % freq vector [Hz]
f_kHz = f/1000;                                     % freq vector [kHz]
figure,plot(f_kHz,F2)                               % plot 
axis([0  max(f_kHz) 0 max(F2)])                     % axis details
title('Frequency Domain Representation of Original Recording') 	% labels
xlabel('F (kHz)')                                   % labels
ylabel('Y(F)')                                      % labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) D major chord frequencies for notes D3, A3, D4, F#4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberOfNotes = 4;
D3 = 146.83;                            % freq of note D3 [Hz]
D3_int = round(D3/max(f)*length(f));     % associated integer to above freq
A3 = 220.00;                            % freq of note A3 [Hz]
A3_int = round(A3/max(f)*length(f));     % associated integer to above freq
D4 = 293.66;                            % freq of note D4 [Hz]
D4_int = round(D4/max(f)*length(f));     % associated integer to above freq
F_sharp_4 = 369.99;                     % freq of note F#4 [Hz]
F_sharp_4_int = ...
    round(F_sharp_4/max(f)*length(f));   % associated integer to above freq

note_freq = [D3 A3 D4 F_sharp_4];       % vector of all note freqs
note_freq_int = ...
  [D3_int A3_int D4_int F_sharp_4_int]; % vector of all int note freqs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (7) Find true maxima and corresponding frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
note_freq_amp = zeros(1,numberOfNotes); % amplitude of each frequency
for i = 1:numberOfNotes                 % loop through all frequencies
    frequencyIndex = note_freq_int(i);
    maxAmp = 0;
    maxAmpIdx = frequencyIndex;
    
    for freqIdx = frequencyIndex-10:frequencyIndex+10 %check 10 frequency points behind and ahead to find true maxima
        amp = F2(freqIdx);
        if amp > maxAmp
            maxAmp = amp;
            maxAmpIdx = freqIdx;
        end
    end
        
    freq = f_kHz(maxAmpIdx);                        %get the frequency, index, and amplitude
    
    note_freq(i) = freq;
    note_freq_int(i) = maxAmpIdx;
    note_freq_amp(i) = maxAmp;
    
    plot(freq,maxAmp,'r*')                          %plot it
end
ylim([0 max(note_freq_amp) * 1.1])
xlim([0 max(note_freq) * 1.1])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8) Plot Frequency Response of each filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterOrders = [3 5 3 4];
passBandRipples = [2 2 5 8];
deltaFs = [2 0.01 1 1];
names = ["D3","A3","D4","F#4"];

for i = 1:numberOfNotes
    freq = note_freq(i) * 1000;
    filterOrder = filterOrders(i);
    passBandRipple = passBandRipples(i);
    deltaF = deltaFs(i);
    name = names(i);

    ttl = sprintf('Magnitude Response (dB) of %s Passband Filter', name);

    maxF = freq + (deltaF * freq / 100);
    minF = freq - (deltaF * freq / 100);

    [b, a] = cheby1(filterOrder,passBandRipple,(maxF)/(samplingFrequency/2),'low');
    [d, c] = cheby1(filterOrder,passBandRipple,(minF)/(samplingFrequency/2),'high');

    lpf = tf(a,b);
    hpf = tf(c,d);

    bpf = series(lpf,hpf);

    [num,den] = tfdata(bpf, 'v');

    [h,frequencies] = freqz(den,num);
    h = db(h);
    frequencies = frequencies * (samplingFrequency) / (2* pi);

    figure;
    hold on
    plot(frequencies, h);
    xline(maxF, 'r');
    xline(minF, 'g');
    axis tight
    xlim([0, F_sharp_4*1.1]);

    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(ttl);
    
    %Find the point closest to the note frequency
    for note_frequency = note_freq
        fq = note_frequency * 1000;
        xline(fq);
        absDiff = abs(frequencies-fq);
        minPt = min(absDiff);
        plot(fq, h(absDiff == minPt), "r*");
    end
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) Create filter bank with optimal parameters and record frequency amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each value in each array corresponds to the optimal parameters for each
% passband filter
filterOrders = [3 5 3 4];
passBandRipples = [2 2 5 8];
deltaFs = [2 0.01 1 1];
names = ["D3","A3","D4","F#4"];

filteredSignals = zeros(signalLength,numberOfNotes);
filteredFrequencyAmplitude = zeros(numberOfNotes,1);

for i = 1:numberOfNotes
    name = names(i);
    freq = note_freq(i) * 1000;
    filterOrder = filterOrders(i);
    passBandRipple = passBandRipples(i);
    deltaF = deltaFs(i);
    
    maxF = freq + (deltaF * freq / 100);
    minF = freq - (deltaF * freq / 100);
    
    [b, a] = cheby1(filterOrder,passBandRipple,maxF/(samplingFrequency/2),'low');
    y = filter(b,a,inputSignal);
    
    [b, a] = cheby1(filterOrder,passBandRipple,minF/(samplingFrequency/2),'high');
    y = filter(b,a,y);
    
    filteredSignals(:,i) = y;
    
    % Plot the Time Domain Signal
    figure,plot(time,y), axis tight
    hold on
    title(sprintf('Filtered %s Time Domain Signal', name))                  	% labels
    xlabel('t (s)')                             % labels
    ylabel('y(t)')                              % labels
    hold off
    
    Y = fft(y);                             % Discrete Fourier transform

    F1 = abs(Y/signalLength);                        % frequency
    F2 = F1(1:signalLength/2+1);                     % half of frequency
    F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
    
    %Plot the frequency domain signal
    figure,plot(f_kHz,F2)                   % plot 
    axis([0  max(f_kHz) 0 max(F2)])         % axis details
    
    % get the amplitude of the peak frequency component
    pk = max(F2);
    filteredFrequencyAmplitude(i) = pk;
    
    hold on
    xline(note_freq);
    plot(note_freq, F2(note_freq_int), "r*");
    
    title(sprintf('Filtered %s Frequency Domain Signal', name))                  	% labels
    xlabel('F (kHz)')                       % labels
    ylabel('Y(F)')                          % labels
    grid("on");
    axis("tight");
    xlim([0 note_freq(end)*1.1])
    ylim([0 max(F2)*1.1]);
    hold off
end

% inverse of the amplitudes
frequencyAmplitudeScales = ones(numberOfNotes,1) ./ filteredFrequencyAmplitude;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (10) Normalize Filtered Signals with Respect to the Passband Frequency Ampitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalizedSignals = zeros(signalLength,4);
normalizedSignalsFFT = zeros(signalLength,4);
for i = 1:numberOfNotes
    name = names(i);
    filteredSignal = filteredSignals(:,i);
    
    normalizedSignal = filteredSignal * frequencyAmplitudeScales(i);
    normalizedSignals(:,i) = normalizedSignal;

    Y = fft(normalizedSignal);              % Discrete Fourier transform
    F1 = abs(Y/signalLength);               % frequency
    F2 = F1(1:signalLength/2+1);            % half of frequency
    F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
    figure,plot(f_kHz,F2)                   % plot 

    normalizedSignalsFFT(:,i) = Y;
    
    hold on
    xline(note_freq);
    plot(note_freq, F2(note_freq_int), "r*");
    
    title(sprintf('Normalized Filtered %s Frequency Domain Signal', name)) % labels
    xlabel('F (kHz)')                       % labels
    ylabel('Y(F)')                          % labels
    grid("on");
    axis("tight");
    xlim([0 note_freq(end)*1.1])
    ylim([0 max(F2)*1.1]);
    hold ("off");
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (11) Add Signals and Scale in the Frequency Domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compressedSignalFFT = sum(normalizedSignalsFFT,2) / min(frequencyAmplitudeScales);
F1 = abs(compressedSignalFFT/signalLength);               % frequency
F2 = F1(1:signalLength/2+1);            % half of frequency
F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
figure,plot(f_kHz,F2)                   % plot 
title('Frequency Domain Representation of Compressed Signal') % labels
xlabel('F (kHz)')                                   % labels
ylabel('Y(F)')                                      % labels
hold on
xline(note_freq);
plot(note_freq, F2(note_freq_int), "r*");
xlim([0 note_freq(end)*1.1])
ylim([0 max(F2)*1.1]);
hold off

compressedSignal = ifft(compressedSignalFFT);
figure,plot(time,compressedSignal), axis tight
title('Time Domain Representation of Compressed Signal')      % labels
xlabel('t (s)')                             % labels
ylabel('y(t)')                              % labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (12) Add Signals and Scale in the Time Domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compressedSignal = sum(normalizedSignals,2) / min(frequencyAmplitudeScales);

figure,plot(time,compressedSignal), axis tight
title('Time Domain Representation of Compressed Signal')      % labels
xlabel('t (s)')                             % labels
ylabel('y(t)')                              % labels

Y = fft(compressedSignal);                          % Discrete Fourier transform
F1 = abs(Y/signalLength);                           % frequency
F2 = F1(1:signalLength/2+1);                        % half of frequency
F2(2:end-1) = 2*F2(2:end-1);                        % Discrete Fourier transform
figure,plot(f_kHz,F2)                               % plot 
xlim([0 note_freq(end)*1.1])
ylim([0 max(F2)*1.1]);
title('Frequency Domain Representation of Compressed Signal') % labels
xlabel('F (kHz)')                                   % labels
ylabel('Y(F)')                                      % labels
hold on
xline(note_freq);
plot(note_freq, F2(note_freq_int), "r*");
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (13) Simulink Model Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simin = [time' inputSignal];

filterOrders = [3 5 3 4];
passBandRipples = [2 2 5 8];
deltaFs = [2 0.01 1 1];

sys = tf(zeros(1,1,numberOfNotes*2));

for i = 1:numberOfNotes
    freq = note_freq(i) * 1000;
    filterOrder = filterOrders(i);
    passBandRipple = passBandRipples(i);
    deltaF = deltaFs(i);
    filteredSignalPk = filteredFrequencyAmplitude(i);
    
    maxF = freq + (deltaF * freq / 100);
    minF = freq - (deltaF * freq / 100);
    
    [b, a] = cheby1(filterOrder,passBandRipple,(freq+deltaF)/(samplingFrequency/2),'low');
    lpf = tf(a,b);
    sys(:,:,i) = lpf;
    
    [d, c] = cheby1(filterOrder,passBandRipple,(freq-deltaF)/(samplingFrequency/2),'high');
    hpf = tf (c,d);

    sys(:,:,i+4) = hpf;
end

D3LPF = sys(:,:,1);
[D3LPFNUM,D3LPFDEN] = tfdata(D3LPF, 'v');

D3HPF = sys(:,:,5);
[D3HPFNUM,D3HPFDEN] = tfdata(D3HPF, 'v');

A3LPF = sys(:,:,2);
[A3LPFNUM,A3LPFDEN] = tfdata(A3LPF, 'v');

A3HPF = sys(:,:,6);
[A3HPFNUM,A3HPFDEN] = tfdata(A3HPF, 'v');

D4LPF = sys(:,:,3);
[D4LPFNUM,D4LPFDEN] = tfdata(D4LPF, 'v');

D4HPF = sys(:,:,7);
[D4HPFNUM,D4HPFDEN] = tfdata(D4HPF, 'v');

FS4LPF = sys(:,:,4);
[FS4LPFNUM,FS4LPFDEN] = tfdata(FS4LPF, 'v');

FS4HPF = sys(:,:,8);
[FS4HPFNUM,FS4HPFDEN] = tfdata(FS4HPF, 'v');

D3Gain = frequencyAmplitudeScales(1);
A3Gain = frequencyAmplitudeScales(2);
D4Gain = frequencyAmplitudeScales(3);
FS4Gain = frequencyAmplitudeScales(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (14) Equivalent Transfer Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filterOrders = [3 5 3 4];
passBandRipples = [2 2 5 8];
deltaFs = [2 0.01 1 1];

sys = tf(zeros(1,1,numberOfNotes));

for i = 1:numberOfNotes
    freq = note_freq(i) * 1000;
    filterOrder = filterOrders(i);
    passBandRipple = passBandRipples(i);
    deltaF = deltaFs(i);
    filteredSignalPk = filteredFrequencyAmplitude(i);
    
    maxF = freq + (deltaF * freq / 100);
    minF = freq - (deltaF * freq / 100);
    
    [b, a] = cheby1(filterOrder,passBandRipple,(freq+deltaF)/(samplingFrequency/2),'low');
    lpf = tf(a,b);
    
    [d, c] = cheby1(filterOrder,passBandRipple,(freq-deltaF)/(samplingFrequency/2),'high');
    hpf = tf (c,d);

    bpf = series(lpf,hpf) / filteredSignalPk;
    sys(:,:,i) = bpf;
end

equivalentTF1 = parallel(sys(:,:,1), sys(:,:,2));
equivalentTF2 = parallel(sys(:,:,3), sys(:,:,4));

equivalentTF = parallel(equivalentTF1, equivalentTF2);
[num,den] = tfdata(equivalentTF, 'v');

[h,frequencies] = freqz(den,num);
h = db(h);
frequencies = frequencies * (samplingFrequency) / (2* pi);

figure;
hold on
plot(frequencies, h);
axis tight
xlim([0, F_sharp_4*1.1]);

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Equivalent Transfer Function Frequency Response');

for note_frequency = note_freq
    fq = note_frequency * 1000;
    xline(fq);
    absDiff = abs(frequencies-fq);
    minPt = min(absDiff);
    plot(fq, h(absDiff == minPt), "r*");
end
hold off

% compressedSignal2 = filter(den,num,inputSignal);
% figure,plot(time,compressedSignal2), axis tight
% title('Time Domain Representation of Compressed Signal')      % labels
% xlabel('t (s)')                             % labels
% ylabel('y(t)')                              % labels
% 
% Y = fft(compressedSignal2);                         % Discrete Fourier transform
% F1 = abs(Y/signalLength);                           % frequency
% F2 = F1(1:signalLength/2+1);                        % half of frequency
% F2(2:end-1) = 2*F2(2:end-1);                        % Discrete Fourier transform
% figure,plot(f_kHz,F2)                               % plot 
% xlim([0 note_freq(end)*1.1]);
% maxY = max(F2)*1.1;
% ylim([0 maxY]);
% title('Frequency Domain Representation of Compressed Signal') % labels
% xlabel('F (kHz)')                                   % labels
% ylabel('Y(F)')                                      % labels
% hold on
% xline(note_freq);
% plot(note_freq, F2(note_freq_int), "r*");
% hold off
