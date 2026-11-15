function p = ensurePool_v001(nWorkers)
% ENSUREPOOL_V001  Return the current parallel pool, opening one only if none exists. NEVER deletes or resizes an existing pool.
%   Rule (Dagmar, 2026-11-07): called functions must not restart the caller's pool. A restart inside a cell (checkAdditiveVsRetimed_v003:
%   "leave one core free", delete + parpool) killed a BlueBEAR worker for memory on 2026-11-07. Size is the caller's business: the launcher
%   opens the pool once, at full size on HPC (all allocated cores; no desktop "one core free" rule), and everything below reuses it.
%   p = ensurePool_v001()           existing pool, or a new "Processes" pool at the profile's default size (all cores)
%   p = ensurePool_v001(nWorkers)   existing pool (any size, reported, not changed), or a new pool of nWorkers
%   Fraser, D.S. (2026)  v001
p = gcp("nocreate");
if isempty(p)
    if nargin < 1 || isempty(nWorkers), p = parpool("Processes"); else, p = parpool("Processes", nWorkers); end
end
fprintf("Pool: %d workers on %d cores (reused if it existed)\n", p.NumWorkers, feature("numcores"));
end
