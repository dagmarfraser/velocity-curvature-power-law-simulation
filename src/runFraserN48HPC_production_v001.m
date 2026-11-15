function runFraserN48HPC_production_v001()
%RUNFRASERN48HPC_PRODUCTION_V001  One press: FraserN48 under the production template only (2026-11-03).
%   Open in the MATLAB GUI on BlueBEAR (sinteractive, 72 cores) and press Run. About 6.5 h. If it PAUSES, press Run again
%   in the next session: finished chunks are kept. The other template runs from its own wrapper in a parallel session.
%   Fraser, D.S. (2026)  v001
runFraserN48HPC_v002("production");
end
