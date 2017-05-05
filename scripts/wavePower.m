function powerCurve = wavePower(waveform, agc)
%WAVEPOWER Calculate the power of a waveform
%   MP = wavePower(WAVEFORM, AGC) returns the power curve in dB given
%   a waveform and the automatic gain control in dB (agc).
%
%   Source: "Sea Ice Leads Detection Using SARAL/AltiKa Alimeter"
%           by E. Zakharova, S. Fleury, K. Guerreiro, S. Willmes,
%              F. Rémy, A. Kouraev, and G. Heinemann

   powerCurve = 10 * log10(waveform * 10^(agc/10));
end