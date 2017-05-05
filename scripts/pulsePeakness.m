function PPaltika = pulsePeakness(waveform, Nright)
%PULSEPEAKNESS Calculate the pulse peakness of a waveform
%   PP = pulsePeakness(WAVEFORM, Nright) returns the pulse peakness
%   of a waveform given Nright.
%
%   Source: "Sea Ice Leads Detection Using SARAL/AltiKa Alimeter"
%           by E. Zakharova, S. Fleury, K. Guerreiro, S. Willmes,
%              F. Rémy, A. Kouraev, and G. Heinemann

PPaltika = max(waveform) * Nright / sum(waveform);
end