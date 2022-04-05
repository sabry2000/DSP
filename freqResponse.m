%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (10) Plot Frequency Response of each filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterOrders = [3 5 3 4];
passBandRipples = [2 2 5 2];
deltaFs = [2 0.01 1 1];

for i = 1:numberOfNotes
    freq = note_freq(i) * 1000;
    filterOrder = filterOrders(i);
    passBandRipple = passBandRipples(i);
    deltaF = deltaFs(i);

    maxF = freq + (deltaF * freq / 100);
    minF = freq - (deltaF * freq / 100);

    [b, a] = cheby1(filterOrder,passBandRipple,(maxF)/(samplingFrequency/2),'low');
    [d, c] = cheby1(filterOrder,passBandRipple,(minF)/(samplingFrequency/2),'high');

    lpf = tf(a,b);
    hpf = tf(c,d);

    bpf = series(lpf,hpf);

    [num,den] = tfdata(bpf, 'v');

    figure;
    freqz(den,num);
    xline((maxF)/(samplingFrequency/2));
    xline((minF)/(samplingFrequency/2));
end