D3 = 146.83;                            % freq of note D3 [Hz]
A3 = 220.00;                            % freq of note A3 [Hz]
D4 = 293.66;                            % freq of note D4 [Hz]
F_sharp_4 = 369.99;                     % freq of note F#4 [Hz]

Fs = 14000;

frequencies = [D3 A3 D4 F_sharp_4];

filterOrders = 1:5;
passRippleRipples = [0.1,0.5,1,2,5,10];
deltaFs = [0,1,5,10,20];
names = ["D3","A3","D4","F#4"];

% filterOrders = 1:1;
% passRippleRipples = [0.1,0.5,1,2,5,10];
% deltaFs = 0:0;
% names = ["D3","A3","D4","F#4"];

% optimalFilterOrders = zeros(4,1);
% optimalbandpassRipples = zeros(4,1);
% optimalDeltaFs = zeros(4,1);

for i = 1:length(frequencies)
    freq = frequencies(i);
    name = names(i);
    v = VideoWriter(sprintf('plots_%s.avi',name));
    open(v);
    for filterOrder = filterOrders
        for passbandRipple = passRippleRipples
            for deltaF = deltaFs
                try
                    
                    maxF = freq + deltaF;
                    minF = freq - deltaF;

                    f = figure;
                    figureName = sprintf("n = %d, Rp = %f, minF = %f, maxF = %f",...
                        filterOrder, passbandRipple, minF, maxF);
                    f.Name = figureName;

                    [b,a] = cheby1(filterOrder,passbandRipple,(maxF)/(Fs/2), "low");
                    [d,c] = cheby1(filterOrder,passbandRipple,(minF)/(Fs/2), "high");

                    sys = series(tf(a,b),tf(c,d));
                    [num,den] = tfdata(sys, 'v');

                    [h,freqs] = freqz(den,num,[],Fs);
                    h = db(h);

                    idx = find(freqs < (10 * freqs));

                    plot(freqs(idx),h(idx));

                    ttl = sprintf('Magnitude (dB) Response of %s Passband Filter \n n = %d, Rp = %f, minF = %f, maxF = %f',...
                        name, filterOrder, passbandRipple, minF, maxF);

                    hold("on");
                    for note_frequency = frequencies
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
