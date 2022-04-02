numberOfNotes = 4;
samplingFrequency = 48000;
note_freq = [0.146277993675304	0.220588868878874	0.294624007996347	0.369486355372115];

filterOrders = 1:5;
passRippleRipples = [0.1,0.5,1,2,5,10];
deltaFs = [0,1,5,10,20];
names = ["D3","A3","D4","F#4"];

for i = 1:numberOfNotes
    freq = note_freq(i) * 10^3;
    name = names(i);
    v = VideoWriter(sprintf('%s.avi',name));
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

                    [b,a] = cheby1(filterOrder,passbandRipple,(maxF)/(samplingFrequency/2), 'low');
                    [d,c] = cheby1(filterOrder,passbandRipple,(minF)/(samplingFrequency/2), 'high');
                    
                    lpf = tf(a,b);
                    hpf = tf(c,d);

                    bpf = series(lpf,hpf);
                    [num,den] = tfdata(bpf, 'v');

                    [h,freqs] = freqz(den,num,[],samplingFrequency);
                    h = db(h);

                    idx = find(freqs < (10 * freqs));

                    plot(freqs(idx),h(idx));

                    ttl = sprintf('Frequency Response(dB) of %s Passband Filter \n n = %d, Rp = %.1f, deltaF = %d, minF = %.2f, maxF = %.2f',...
                        name, filterOrder, passbandRipple, minF, maxF);

                    hold("on");
                    for note_frequency = note_freq
                        noteF = note_frequency * 10^3;
                        xline(noteF);
                        absDiff = abs(freqs-noteF);
                        minPt = min(absDiff);
                        plot(noteF, h(absDiff == minPt), "r*");
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