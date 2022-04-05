%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8) Compute ideal attenuation for each frequency component associated with each filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

attenuationMatrix = zeros(numberOfNotes);
for i = 1:numberOfNotes
    attenuationMatrix(i,:) = note_freq_amp / (note_freq_amp(i)/2);
end

attenuationMatrix = 10 * log10(ones(numberOfNotes) ./ attenuationMatrix);
attenuationMatrix(attenuationMatrix > -3) = -3;
z = zeros(numberOfNotes,1);
attenuationMatrix(1:(numberOfNotes+1):end) = z;