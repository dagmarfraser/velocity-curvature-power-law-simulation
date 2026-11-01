function [nMaj, nMin] = shapedCommonSurrogate_v001(rm, rn, FS, aMaj, aMin)
%SHAPEDCOMMONSURROGATE_V001  Xa noise (spec 10.55 as amended 10.56): colour AND cross-axis relation kept, template-locked phase replaced.
%   Along axis = shaped_xu(rm) (empirical amplitude spectrum on Xu fractional phase at alpha aMaj; alpha survives, unlike random phase).
%   Across axis = |FFT(rn)| with phase angle(FFT(nMaj)) + angle(FFT(rn)) - angle(FFT(rm)), so conj(FFT(nMaj)).*FFT(nMin) equals conj(FFT(rm)).*FFT(rn)
%   up to the final SD rescale of nMin to std(rn). Per-axis SD equals the input's. Whether the across-axis IRASA alpha survives is tested in checkR30Arms_v001.
%   Fraser, D.S. (2026)  v001
rm = rm(:); rn = rn(:);
nMaj = generateLoopClosureNoise_v003('shaped_xu', rm, rn, FS, aMaj, aMin);
nMaj = nMaj(:);
if numel(nMaj) ~= numel(rm) || any(~isfinite(nMaj)), error('shapedCommonSurrogate:noise','shaped_xu returned length %d (need %d) or non-finite values (alpha %.2f out of range?)',numel(nMaj),numel(rm),aMaj); end
Zm = fft(rm); Zn = fft(rn); Zs = fft(nMaj);
nMin = real(ifft(abs(Zn).*exp(1i*(angle(Zs) + angle(Zn) - angle(Zm)))));
nMin = nMin*std(rn)/std(nMin);
end
