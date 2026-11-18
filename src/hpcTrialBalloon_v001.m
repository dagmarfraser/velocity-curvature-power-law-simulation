function ok = hpcTrialBalloon_v001(deadline)
%HPCTRIALBALLOON_V001  Ten-minute end-to-end test of the BlueBEAR batch route before any real v022 job (2026-11-15).
%   Run by hpcTrialBalloon_v001.sh (sbatch, bbshort QOS) or by hand: ok = hpcTrialBalloon_v001("2026-11-15T09:00:00").
%   Checks, each printed as PASS or FAIL with the reason; nothing stops at the first failure:
%     1. environment: host, MATLAB version, Slurm job id, cores, $TMPDIR, Parallel Computing Toolbox licence
%     2. parallel pool opened the BEAR way (parcluster('local'), JobStorageLocation = $TMPDIR, size = cores Slurm gave),
%        and a parfor over all workers returns the right answer
%     3. every function and data file the v022 endMatched run needs is on the path / in place (src/ and src/functions/)
%     4. Cook CTRL imports (94 trials expected) and its noise file loads
%     5. the end-matched pieces run on one real trial (first Cook CTRL trial, production template): segment indices and
%        one noise draw, compared with the reference values computed on the M5 (R2026a) on 2026-11-15; a mismatch in the
%        segment is reported as FAIL (deterministic code), in the noise draw as WARN (alpha selection can flip on rounding)
%     6. the deadline string from Slurm parses (the launcher's pause-before-time-limit logic)
%     7. write access: a small result file in results/hpcTrialBalloon/ (kept, named by job id)
%   ok is true when nothing FAILed. In batch, a false ok makes MATLAB exit non-zero, so sacct shows FAILED.
arguments
    deadline (1,1) string = ""
end
here = fileparts(mfilename('fullpath'));  addpath(here, genpath(fullfile(here, 'functions')));
R = strings(0, 3);   % check, status, detail
    function rec(name, status, detail), R(end + 1, :) = [string(name), string(status), string(detail)]; fprintf('%-34s %-4s %s\n', name, status, detail); end

% 1. environment
[~, host] = system('hostname');  host = strtrim(host);
jid = getenv('SLURM_JOB_ID');  ncpu = getenv('SLURM_CPUS_PER_TASK');  ntask = getenv('SLURM_TASKS_PER_NODE');  tmp = getenv('TMPDIR');
rec('host / MATLAB', "INFO", sprintf('%s, MATLAB %s', host, version));
rec('Slurm job / cores', "INFO", sprintf('job %s, cpus-per-task %s, tasks-per-node %s, numcores %d', jid, ncpu, ntask, feature('numcores')));
if ~isempty(tmp) && isfolder(tmp), rec('TMPDIR', "PASS", tmp); else, rec('TMPDIR', "FAIL", sprintf('"%s" is not a folder', tmp)); end
if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel')), rec('Parallel Computing Toolbox', "PASS", 'licensed and installed');
else, rec('Parallel Computing Toolbox', "FAIL", 'not licensed or not installed'); end

% 2. pool, BEAR's recipe
try
    n = str2double(ncpu);  if ~isfinite(n), n = str2double(ntask); end;  if ~isfinite(n), n = feature('numcores'); end
    p = gcp('nocreate');
    if isempty(p)
        pc = parcluster('local');  if ~isempty(tmp), pc.JobStorageLocation = tmp; end
        pc.NumWorkers = n;  p = parpool(pc, n);
    end
    s = zeros(1, 4 * p.NumWorkers);  parfor k = 1:numel(s), s(k) = k^2; end
    if p.NumWorkers == n && isequal(s, (1:numel(s)).^2), rec('parallel pool', "PASS", sprintf('%d workers, storage %s', p.NumWorkers, p.Cluster.JobStorageLocation));
    else, rec('parallel pool', "FAIL", sprintf('%d workers (expected %d) or parfor result wrong', p.NumWorkers, n)); end
catch ME
    rec('parallel pool', "FAIL", ME.message);
end

% 3. files the v022 endMatched run needs
need = ["runLoopClosureFftnoise_v022", "endMatchedNoise_v001", "endToEndSegment_v004", "alphaPWLW_v001", "dfaAlpha_v001", ...
    "ensurePool_v001", "differentiateKinematicsEBR", "curvatureKinematicEBR", "regressDataEBR", "iraAlphaSigma_v001", ...
    "generateLoopClosureNoise_v003", "findBothBranches_v008", "importDB_cook_v002", "importDB_hickman_v003", "importDB_fraser_v002", ...
    "checkV022Acceptance_v001", "runLoopClosureFftnoise_v021", "runV022HPC_v001"];
miss = need(arrayfun(@(f) isempty(which(f)), need));
if isempty(miss), rec('functions on path', "PASS", sprintf('%d of %d found', numel(need), numel(need)));
else, rec('functions on path', "FAIL", "missing: " + strjoin(miss, ", ")); end
nf = ["noiseCharacterisation_cook.mat", "noiseCharacterisation_cookASD.mat"];
nm = nf(~arrayfun(@(f) isfile(fullfile(here, f)), nf));
if isempty(nm), rec('noise files (Cook)', "PASS", strjoin(nf, ", ")); else, rec('noise files (Cook)', "FAIL", "missing in src/: " + strjoin(nm, ", ")); end

% 4. data import and noise file
tr = [];
try
    tr = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
    if numel(tr) == 94, rec('Cook CTRL import', "PASS", '94 trials'); else, rec('Cook CTRL import', "FAIL", sprintf('%d trials (expected 94)', numel(tr))); end
    b = load(fullfile(here, 'noiseCharacterisation_cook.mat'), 'bioResults');
    rec('Cook noise file loads', "PASS", sprintf('%d rows', height(b.bioResults)));
catch ME
    rec('Cook CTRL import / noise file', "FAIL", ME.message);
end

% 5. end-matched pieces on one real trial, against the M5 reference
REF_SEG = [1987 4001; 1962 3976];  REF_NOISE = 97.66580567;   % M5 (UOB260130A, R2026a Update 5), 2026-11-15 08:21; L 2015
try
    if isempty(tr), error('balloon:noTrial', 'no trials imported'); end
    [segL, segInfo, n1, info1] = endMatchedOnTrial(tr(1));
    rec('end-matched segment (trial 1)', "INFO", sprintf('L %d (%.3f of residual), maj [%d %d], min [%d %d]', segInfo.L, segInfo.Lfrac, segL(1, 1), segL(1, 2), segL(2, 1), segL(2, 2)));
    rec('end-matched noise draw (trial 1)', "INFO", sprintf('accepted maj %d min %d, draws %d/%d, first sample %.10g', info1.acceptedMaj, info1.acceptedMin, info1.nDrawMaj, info1.nDrawMin, n1(1)));
    if all(isfinite(REF_SEG(:)))
        if isequal(segL, REF_SEG), rec('segment matches M5', "PASS", ''); else, rec('segment matches M5', "FAIL", sprintf('M5 %s', mat2str(REF_SEG))); end
    end
    if isfinite(REF_NOISE)
        if abs(n1(1) - REF_NOISE) <= 1e-9 * max(1, abs(REF_NOISE)), rec('noise draw matches M5', "PASS", ''); else, rec('noise draw matches M5', "WARN", sprintf('M5 %.10g', REF_NOISE)); end
    end
catch ME
    rec('end-matched pieces', "FAIL", ME.message);
end

% 6. deadline parsing (as runV022HPC_v001)
try
    if deadline == "", rec('deadline parse', "INFO", 'no deadline given');
    else
        dl = datetime(deadline, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss");
        rec('deadline parse', "PASS", sprintf('%s (job ends in %.0f min)', char(dl), minutes(dl - datetime('now'))));
    end
catch ME
    rec('deadline parse', "FAIL", ME.message);
end

% 7. write access
try
    od = fullfile(here, '..', 'results', 'hpcTrialBalloon');  if ~isfolder(od), mkdir(od); end
    tag = jid;  if isempty(tag), tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')); end
    f = fullfile(od, "hpcTrialBalloon_v001_" + tag + ".txt");
    writematrix(R, f, 'Delimiter', '\t', 'FileType', 'text');
    rec('write access', "PASS", f);
catch ME
    rec('write access', "FAIL", ME.message);
end

ok = ~any(R(:, 2) == "FAIL");
fprintf('\nTRIAL BALLOON %s: %d PASS, %d FAIL, %d WARN\n', string(ok).replace("true", "PASSED").replace("false", "FAILED"), nnz(R(:, 2) == "PASS"), nnz(R(:, 2) == "FAIL"), nnz(R(:, 2) == "WARN"));
if ~ok && batchStartupOptionUsed, error('balloon:failed', 'trial balloon failed: see the FAIL lines above'); end
end

function [segL, info, n1, info1] = endMatchedOnTrial(t)
% the v022 steps for one trial: residual (production template, the runner's own local functions are not reachable, so the
% same template code as endMismatchSensitivity_v006 is used), common-length segments, alpha, one noise draw on a fixed stream
x = double(t.x(:));  y = double(t.y(:));  FS = 133;
f0 = estimateF0_b(x, y, FS);  [xR, yR] = templateSubtract_b(x, y, FS, f0, 4, 4);
xc = x - mean(x);  yc = y - mean(y);  [V, D] = eig(cov(xc, yc));  [~, ord] = sort(diag(D), 'descend');  V = V(:, ord);  th = atan2(V(2,1), V(1,1));
rMaj = xR*cos(th) + yR*sin(th);  rMin = -xR*sin(th) + yR*cos(th);
[i0, i1, info] = endToEndSegment_v004([rMaj, rMin], Tol=1, Mode="commonLength");  segL = [i0(:), i1(:)];
sMaj = rMaj(i0(1):i1(1));  sMin = rMin(i0(2):i1(2));
a = [alphaPWLW_v001(sMaj, FS, [1 20]), alphaPWLW_v001(sMin, FS, [1 20])];
rs = RandStream("threefry4x64_20", "Seed", 1729);  rs.Substream = 1;  old = RandStream.setGlobalStream(rs);
[n1, ~, info1] = endMatchedNoise_v001(sMaj, sMin, FS, a(1), a(2));
RandStream.setGlobalStream(old);
end

% ---- copied verbatim from runLoopClosureFftnoise_v021.m (production template path) ----
function f0 = estimateF0_b(x, y, fs)
    N    = numel(x); nfft = 2^nextpow2(4*N);
    Xf   = abs(fft(detrend(x(:), 'linear'), nfft));
    Yf   = abs(fft(detrend(y(:), 'linear'), nfft));
    fAx  = (0:nfft-1)' * fs / nfft;
    band = fAx > 0.1 & fAx < 5; idx = find(band);
    [~, pkX] = max(Xf(band)); [~, pkY] = max(Yf(band));
    f0x = fAx(idx(pkX)); f0y = fAx(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2
        f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band))
        f0 = f0x;
    else
        f0 = f0y;
    end
end

function [resX, resY] = templateSubtract_b(x, y, fs, f0, nH, nCW)
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);
    if win > N, starts = []; end
    for s = starts
        e   = s + win - 1;
        tw  = (0:win-1)'/fs;
        han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
        D   = [ones(win,1), tw];
        for h = 1:nH
            D(:,end+1) = cos(2*pi*h*f0*tw); %#ok<AGROW>
            D(:,end+1) = sin(2*pi*h*f0*tw); %#ok<AGROW>
        end
        Dw        = D .* han;
        bX        = Dw \ (x(s:e) .* han);
        bY        = Dw \ (y(s:e) .* han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;
        resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;
        ww(s:e)   = ww(s:e) + han;
    end
    ok = ww > 0;
    resX(ok) = resX(ok) ./ ww(ok);
    resY(ok) = resY(ok) ./ ww(ok);
    cl   = max(round(win/4), 1);
    cs   = cl;
    ce   = min(max(N - cl, cs + 100), N);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
