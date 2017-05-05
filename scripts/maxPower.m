function value = maxPower(waveform, agc)
%MAXPOWER Calculate the maximum power of a waveform
%   MP = maxPower(WAVEFORM, AGC) returns max power in dB given
%   a waveform and the automatic gain control in dB (agc).
%
%   Source: "Sea Ice Leads Detection Using SARAL/AltiKa Alimeter"
%           by E. Zakharova, S. Fleury, K. Guerreiro, S. Willmes,
%              F. Rémy, A. Kouraev, and G. Heinemann

   value = 10 * log10(max(waveform) * 10^(agc/10));
end