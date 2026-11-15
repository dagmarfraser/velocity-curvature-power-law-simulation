function [ptPerMm, canvasPt] = ipadScaleM4_v001()
%IPADSCALEM4_V001  The one place the iPad Pro 13" M4 scale is defined (docs/PLAN_PencilScaleCorrection_v001.md, section 1.1).
%   StimuliApp -enhanced.csv x,y (Pilot and Fraser) are UIKit POINTS, not device pixels.
%   iPad Pro 13" M4: 2752 x 2064 px, 264.615 ppi, UIKit scale 2, so
%     ptPerMm  = 264.615 / 25.4 / 2 = 5.208965 points per mm
%     canvasPt = [1376 1032]  (landscape point canvas)
%   Evidence: Dagmar's ruler on the on-screen ellipse (about 116 x 55 mm), checkPencilScale_v001 (median trace 629 x 292.5
%   units, i.e. 121 mm if points, 60 mm if pixels), checkPencilFraction_v001 and checkFraserN48Canvas_v001 (every value of
%   all 51 Fraser participants inside the point canvas, filling it). The retired constants 10.41793 (Fraser) and 9.73
%   (Pilot) are pixel densities, about 2x too large as divisors.
%   Callers use 1/ptPerMm as SigmaToMM; nothing else should hard-code an iPad constant.
%   Note: importDB_fraser_v002 and runLoopClosureFftnoise_v020 (both 2026-11-03) carry the rounded literal 5.209, which differs
%   from this value by 7 parts per million; they are left as written (the N=48 HPC run uses them) and later versions call this.
%   Fraser, D.S. (2026)  v001
ptPerMm  = 264.615 / 25.4 / 2;
canvasPt = [1376 1032];
end
