%% checkSpecSetup_v001  Deep Review GD17 items 4-5: geometry extraction, VGF grid, IRASA settings and harmonic-template construction, read from the code.
% Parses runLoopClosureFftnoise_v019.m / _v015.m and ModelAdequacy_Stage1_KitchenSink_v2_001.m (a missing pattern is an error) and records the literals the
% manuscript quotes: ellipse axes from the covariance eigenvalues, the VGF grid, IRASA band and dilation set, template order/window/hop/taper/design, trimming,
% the v019 end-anchored final window (absent from v015), and the length floor (applied to the trimmed residual length M).
% Run from src/. Writes results/checkSpecSetup_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
r19 = fileread(fullfile(ROOT, "src", "runLoopClosureFftnoise_v019.m"));
r15 = fileread(fullfile(ROOT, "src", "runLoopClosureFftnoise_v015.m"));
ks  = fileread(fullfile(ROOT, "src", "ModelAdequacy_Stage1_KitchenSink_v2_001.m"));
need = @(txt, pat, what) assert(~isempty(regexp(txt, pat, "once")), "specSetup:parse", "Pattern not found: %s", what);
need(r19, "\[V, D\]\s+= eig\(cov\(xMM, yMM\)\);", "eig of the 2x2 covariance of mean-centred mm coordinates");
need(r19, "a_mm\s+= sqrt\(2\*max\(lams\(1\), 0\)\);", "a = sqrt(2*lambda1)");
need(r19, "b_mm\s+= sqrt\(2\*max\(lams\(2\), 0\)\);", "b = sqrt(2*lambda2)");
need(r19, "theta\s+= atan2\(V\(2,1\), V\(1,1\)\);", "theta from leading eigenvector");
need(r19, "IRA_FLOW\s+= 1\.0;", "IRASA fLow 1.0");
need(r19, "IRA_FHIGH = min\(20\.0, FS/2 - 1\);", "IRASA fHigh");
need(r19, "IRA_HSET\s+= 1\.1:0\.05:1\.9;", "IRASA h set");
need(r15, "IRA_HSET\s+= 1\.1:0\.05:1\.9;", "IRASA h set (v015)");
need(r19, "\[xResid, yResid, xFit, yFit\] = templateSubtract_local\( \.\.\.\s*\n\s*double\(tr\.x\), double\(tr\.y\), FS, res\.f0, 4, 4\);", "template nH = 4, nCW = 4");
need(r19, "win = round\(nCW/f0\*fs\);", "window = nCW cycles");
need(r19, "hop = round\(win/2\);", "50% hop");
need(r19, "win = max\(win, round\(2/f0\*fs\)\);", "window at least 2 cycles");
need(r19, "han = 0\.5\*\(1 - cos\(2\*pi\*\(0:win-1\)'/\(win-1\)\)\);", "Hann taper");
need(r19, "D\s+= \[ones\(win,1\), tw\];", "design: constant and linear trend");
need(r19, "cos\(2\*pi\*h\*f0\*tw\)", "cosine harmonics"); need(r19, "sin\(2\*pi\*h\*f0\*tw\)", "sine harmonics");
need(r19, "bX\s+= Dw \\ \(x\(s:e\) \.\* han\);", "taper-weighted least squares");
need(r19, "resX\(ok\) = resX\(ok\) \./ ww\(ok\);", "overlap-add normalised by summed taper");
need(r19, "cl\s+= max\(round\(win/4\), 1\);", "trim round(win/4) each end");
need(r19, "ce\s+= min\(max\(N - cl, cs \+ 100\), N\);", "trim end index");
need(r19, "starts\(end\) \+ win - 1 < N, starts\(end\+1\) = N - win \+ 1;", "v019 end-anchored final window");
assert(isempty(regexp(r15, "starts\(end\+1\)", "once")), "specSetup:v015", "v015 unexpectedly has an end-anchored window");
need(r15, "while s \+ win - 1 <= N", "v015 windows by while loop (no end anchor)");
need(r19, "if M < 4\*\(EDGE_CLIP \+ 50\)", "length floor on residual length M");
need(r19, "addParameter\(p, 'EdgeClip',\s+20,", "EdgeClip default 20");
need(ks,  "vgf_full = exp\(linspace\(4\.5, 5\.8, 14\)\);", "VGF grid");
c57 = fileread(fullfile(ROOT, "src", "Toolchain_caller_v057.m"));
c58 = fileread(fullfile(ROOT, "src", "Toolchain_caller_v058.m"));
tf  = fileread(fullfile(ROOT, "src", "Toolchain_func_v032.m"));
need(c57, "cfg\.edgeClip\s*=\s*50;", "registered caller v057: edgeClip 50");
need(c58, "cfg\.edgeClip\s*=\s*50;", "caller v058: edgeClip 50");
need(tf, "edgeClip = double\(cfg\.edgeClip\);", "func reads cfg.edgeClip");
need(tf, "velocityX = dx\(edgeClip:end-edgeClip,2\);", "func clips dx(edgeClip:end-edgeClip)");
need(tf, "if size\(dx, 1\) <= 2\*edgeClip", "func errors if rows <= 2*edgeClip");
need(tf, "if sum\(validIdx\) < 10", "func errors below 10 valid points");
S = struct();
S.harnessEdgeClip = 50; S.harnessMinValid = 10;
S.harnessFailNmax = S.harnessMinValid + 2*S.harnessEdgeClip - 1 - 1;   % rows retained N-2*clip+1 < 10  <=>  N <= 108
assert(S.harnessFailNmax == 108, "specSetup:failN", "Expected failing N of 108 per auditV058Failures_v004");
S.vgfGrid = exp(linspace(4.5, 5.8, 14)); S.vgfLo = S.vgfGrid(1); S.vgfHi = S.vgfGrid(end); S.vgfLnStep = (5.8-4.5)/13;
S.iraFLow = 1.0; S.iraFHighMax = 20; S.iraH = 1.1:0.05:1.9; S.nH = 4; S.nCW = 4; S.edgeClip = 20; S.mFloor = 4*(S.edgeClip + 50);
S.nTemplateCols = 2 + 2*S.nH;
disp(S); save(fullfile(ROOT, "results", "checkSpecSetup_v001.mat"), "-struct", "S");
