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
Fs = 48000;                             % sampling freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) Play sound of raw audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soundsc(inputSignal,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) Graph in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 1/Fs;                                  % sampling period       
Length_y = length(inputSignal(:,1));      	% length of signal
time = (0:Length_y-1)*Ts;                   % time vector
figure,plot(time,inputSignal), axis tight
title('Time Domain Signal of Recording')	% labels
xlabel('t (s)')                             % labels
ylabel('y(t)')                              % labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) Graph with DFT/FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = fft(inputSignal);                               % Discrete Fourier transform
F1 = abs(Y/Length_y);                               % frequency
F2 = F1(1:Length_y/2+1);                            % half of frequency
F2(2:end-1) = 2*F2(2:end-1);                        % Discrete Fourier transform
f = Fs*(0:(Length_y/2))/Length_y;                   % freq vector [Hz]
f_kHz = f/1000;                                     % freq vector [kHz]
figure,plot(f_kHz,F2)                               % plot 
axis([0  max(f_kHz) 0 max(F2)])                     % axis details
title('Frequency Domain Plot of Recorded Signal') 	% labels
xlabel('F (kHz)')                                   % labels
ylabel('Y(F)')                                      % labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) D major chord frequencies for notes D3, A3, D4, F#4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% (7) Find true maxima and corresponding frequncy points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
note_freq_amp = zeros(1,length(note_freq_int)); % amplitude of each frequency
for i = 1:length(note_freq_int)                 % loop through all frequencies
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
    
    freq = f_kHz(maxAmpIdx);                        %get the frequency, its index, and its amplitude
    
    note_freq(i) = freq;
    note_freq_int(i) = maxAmpIdx;
    note_freq_amp(i) = maxAmp;
    
    plot(freq,maxAmp,'r*')                          %plot it
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8) Experimentally study filter bank parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterOrders = 1:5;
passRippleRipples = [0.1,0.5,1,2,5,10];
deltaFs = [0,1,5,10,20];
names = ["D3","A3","D4","F#4"];

for i = 1:length(note_freq)
    freq = note_freq(i) * 10^3;
    name = names(i);
    v = VideoWriter(sprintf('plots_%s.avi',name));
    open(v);
    for filterOrder = filterOrders
        for passbandRipple = passRippleRipples
            for deltaF = deltaFs
                try
                    
                    maxF = freq + (deltaF * freq / 100);
                    minF = freq - (deltaF * freq / 100);

                    f = figure;
                    figureName = sprintf("n = %d, Rp = %f, minF = %f, maxF = %f",...
                        filterOrder, passbandRipple, minF, maxF);
                    f.Name = figureName;

                    [b,a] = cheby1(filterOrder,passbandRipple,(maxF)/(Fs/2), "low");
                    [d,c] = cheby1(filterOrder,passbandRipple,(minF)/(Fs/2), "high");
                    
                    lpf = tf(a,b);
                    hpf = tf(c,d);

                    bpf = series(lpf,hpf);
                    [num,den] = tfdata(bpf, 'v');

                    [h,freqs] = freqz(den,num,[],Fs);
                    h = db(h);

                    idx = find(freqs < (10 * freqs));

                    plot(freqs(idx),h(idx));

                    ttl = sprintf('Frequency Response(dB) of %s Passband Filter \n n = %d, Rp = %.1f, deltaF = %d, minF = %.2f, maxF = %.2f',...
                        name, filterOrder, passbandRipple, minF, maxF);

                    hold("on");
                    for note_frequency = note_freq
                        note_frequency = note_frequency * 10^3;
                        xline(note_frequency);
                        absDiff = abs(freqs-note_frequency);
                        minPt = min(absDiff);
                        plot(note_frequency, h(absDiff == minPt), "r*");
                    end
                    yline(-3);
                    xlabel('Frequency (Hz)');
                    ylabel('Magnitude (dB)');
                    title(ttl);
                    grid("on");
                    axis("tight");
                    ylim([-10 1]);
                    hold ("off");

                    frame = getframe(f);
                    close(f);
                    writeVideo(v,frame);
                catch ME
                    disp(ME)
                    fprintf("%s: Failed\n", name);
                end
            end
        end
    end
    close(v);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) Create filter bank with optimal parameters and record frequency amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterOrders = [3 3 3 4];
passBandRipples = [0.5 1 2 0.5];
deltaFs = [20 1 20 5];

filteredSignals = zeros(length(inputSignal),4);
filteredFrequencyAmplitude = zeros(4,1);

for i = 1:length(note_freq)
    freq = note_freq(i) * 1000;
    filterOrder = filterOrders(i);
    passBandRipple = passBandRipples(i);
    deltaF = deltaFs(i);
    
    maxF = freq + (deltaF * freq / 100);
    minF = freq - (deltaF * freq / 100);
    
    [b, a] = cheby1(filterOrder,passBandRipple,(freq+deltaF)/(Fs/2),'low');
    y = filter(b,a,inputSignal);
    
    [b, a] = cheby1(filterOrder,passBandRipple,(freq-deltaF)/(Fs/2),'high');
    y = filter(b,a,y);
    
    filteredSignals(:,i) = y;
    
    figure,plot(time,y), axis tight
    
    Y = fft(y);                             % Discrete Fourier transform
    len = length(y(:,1));                   % length of signal

    F1 = abs(Y/len);                        % frequency
    F2 = F1(1:len/2+1);                     % half of frequency
    F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
    figure,plot(f_kHz,F2)                   % plot 
    axis([0  max(f_kHz) 0 max(F2)])         % axis details
    
    pk = max(F2);
    filteredFrequencyAmplitude(i) = pk;
    
    hold on
    for j = 1:length(note_freq)
        note_frequency = note_freq(j);
        idx = note_freq_int(j);
        
        xline(note_frequency);
        plot(note_frequency, F2(idx), "r*");
    end
    
    title('Y(F) vs F')                  	% labels
    xlabel('F (kHz)')                       % labels
    ylabel('Y(F)')                          % labels
    grid("on");
    axis("tight");
    ylim([0 pk*1.5]);
    hold ("off");
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) Normalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalizedSignals = zeros(length(inputSignal),4);
for i = 1:length(note_freq)
    filteredSignal = filteredSignals(:,i);
    filteredSignalPk = filteredFrequencyAmplitude(i);
    
    normalizedSignal = filteredSignal/filteredSignalPk;
    normalizedSignals(:,i) = normalizedSignal;
end
%%
Y = fft(normalizedSignal);                             % Discrete Fourier transform
len = length(normalizedSignal(:,1));                   % length of signal

F1 = abs(Y/len);                        % frequency
F2 = F1(1:len/2+1);                     % half of frequency
F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
figure,plot(f_kHz,F2)                   % plot 
axis([0  max(f_kHz) 0 max(F2)])         % axis details

hold on
for j = 1:length(note_freq)
    note_frequency = note_freq(j);
    idx = note_freq_int(j);

    xline(note_frequency);
    plot(note_frequency, F2(idx), "r*");
end

title('Y(F) vs F')                  	% labels
xlabel('F (kHz)')                       % labels
ylabel('Y(F)')                          % labels
grid("on");
axis("tight");
ylim([0 1.1]);
hold ("off");
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) Add Signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compressedSignal = zeros(length(inputSignal),1);

for j = 1:length(note_freq)
    normalizedSignal = normalizedSignals(:,j);
    compressedSignal = compressedSignal + normalizedSignal;
end

figure,plot(time,compressedSignal), axis tight
title('Compressed Time Domain Signal')      % labels
xlabel('t (s)')                             % labels
ylabel('y(t)')                              % labels

Y = fft(compressedSignal);                               % Discrete Fourier transform
F1 = abs(Y/Length_y);                               % frequency
F2 = F1(1:Length_y/2+1);                            % half of frequency
F2(2:end-1) = 2*F2(2:end-1);                        % Discrete Fourier transform
figure,plot(f_kHz,F2)                               % plot 
axis([0  max(f_kHz) 0 max(F2)])                     % axis details
title('Frequency Domain Plot of Recorded Signal') 	% labels
xlabel('F (kHz)')                                   % labels
ylabel('Y(F)')                                      % labels
hold on
for j = 1:length(note_freq)
    note_frequency = note_freq(j);
    idx = note_freq_int(j);

    xline(note_frequency);
    plot(note_frequency, F2(idx), "r*");
end
hold off
%%
Ts = 1/ Fs;

filteredFrequencyAmplitudeX = zeros(4,1);
normalizedSignals = zeros(length(inputSignal),4);


for i = 1:length(note_freq)
    filteredSignal = filteredSignals(i);
    filteredSignalPk = filteredFrequencyAmplitude(i);
    
    filteredFFT = fft(filteredSignal);         % Discrete Fourier transform
    normalizedFFT = filteredFFT/filteredSignalPk;
    
    normalizedSignal = ifft(normalizedFFT) * length(filteredSignal);
    normalizedSignals(:,i) = normalizedSignal;
end