function T = designRecsTempoNumbers_v001()
% DESIGNRECSTEMPONUMBERS_V001  Rebuild the numbers of design recommendation 1 (tempo gaps; Findings #255, #258) on N=45 and N=48.
%   R31b (no true difference): spurious beta_obs difference near 0.4 and 0.7 doublings, the gap at which it reaches the MDC
%   (linear interpolation between grid cells, Fraser W and B), worst adjusted error up to and beyond 0.75 doublings, share
%   of trials matching discards. R31c S1 (real difference): passes per outcome and estimator, worst unadjusted error, the
%   Fraser W k4 realised gap. N=48 R31c: Fraser cells from results/N48/tempoSim/analyseR31cS1_v002_N48.mat, Hickman PLAC
%   cells from the N=45 file on the RDS (PLAC unchanged by N=48). Saves results/N48/results/designRecsTempoNumbers_v001.mat.
%   Run from src/.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename('fullpath')));  MDC = 0.03;
rds = '/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/results/tempoSim';
b45 = load(fullfile(ROOT, 'results', 'tempoSim', 'analyseR31b_v001.mat')).rows;
b48 = load(fullfile(ROOT, 'results', 'N48', 'tempoSim', 'analyseR31b_v001_N48.mat')).rows;
c45 = load(fullfile(rds, 'analyseR31cS1_v002.mat')).out.perCell;
f48 = load(fullfile(ROOT, 'results', 'N48', 'tempoSim', 'analyseR31cS1_v002_N48.mat')).out.perCell;
if any(f48.dataset ~= "Fraser"), error('drt:n48', '%s', 'N=48 S1 file should hold Fraser cells only'), end
c48 = [f48; c45(c45.dataset ~= "Fraser", :)];
T = struct();
for base = ["N45" "N48"]
    if base == "N45", B = b45; C = c45; else, B = b48; C = c48; end
    o = B(B.outcome == "betaObs", :);  F = o(o.dataset == "Fraser", :);
    r.obsNear04 = [F.biasUnadj(F.design == "W" & F.k == 4), F.biasUnadj(F.design == "B" & F.k == 1)];
    r.gapNear04 = [F.deltaLog2f0(F.design == "W" & F.k == 4), F.deltaLog2f0(F.design == "B" & F.k == 1)];
    r.obsNear07 = F.biasUnadj(F.design == "B" & F.k == 4);  r.gapNear07 = F.deltaLog2f0(F.design == "B" & F.k == 4);
    for dz = ["W" "B"]
        G = sortrows(F(F.design == dz, :), 'deltaLog2f0');  r.("mdcGap" + dz) = interp1(G.biasUnadj, G.deltaLog2f0, MDC);
    end
    r.maxAdjUpTo075 = max(abs(o.biasAdj(o.deltaLog2f0 <= 0.75)));
    r.adjBeyond1 = o.biasAdj(o.deltaLog2f0 > 1);
    r.discardPct = 100 * (1 - [min(o.fracKeptMatched(o.k > 0)), max(o.fracKeptMatched)]);  r.discardPct = sort(100 * (1 - [max(o.fracKeptMatched), min(o.fracKeptMatched)]));
    g = C(C.outcome == "betaGen", :);  ob = C(C.outcome == "betaObs", :);
    r.genPass = [nnz(g.passAdj), nnz(g.passMatched), height(g)];  r.obsPass = [nnz(ob.passAdj), nnz(ob.passMatched), height(ob)];
    r.genMaxUnadj = max(g.errUnadj);  r.obsMaxUnadj = max(ob.errUnadj);
    r.fraserWgap = unique(round(C.deltaLog2f0(C.dataset == "Fraser" & C.design == "W"), 2))';
    r.missesObsAdj = C(C.outcome == "betaObs" & ~C.passAdj, {'dataset', 'design', 'k', 'd', 'errAdj'});
    r.missesObsMat = C(C.outcome == "betaObs" & ~C.passMatched, {'dataset', 'design', 'k', 'd', 'errMatched'});
    r.missesGen = C(C.outcome == "betaGen" & (~C.passAdj | ~C.passMatched), {'dataset', 'design', 'k', 'd', 'errAdj', 'errMatched'});
    T.(base) = r;
    fprintf('\n=== %s ===\n', base);
    fprintf('R31b beta_obs spurious: %+.3f / %+.3f at gaps %.2f / %.2f; %+.3f at %.2f\n', r.obsNear04, r.gapNear04, r.obsNear07, r.gapNear07);
    fprintf('R31b MDC reached at gap W %.2f, B %.2f doublings\n', r.mdcGapW, r.mdcGapB);
    fprintf('R31b adjusted: max |bias| %.4f up to 0.75 doublings; beyond 1 doubling %s\n', r.maxAdjUpTo075, mat2str(round(r.adjBeyond1', 4)));
    fprintf('R31b matching discards %.0f to %.0f%% of trials\n', r.discardPct);
    fprintf('R31c beta_gen*: adjusted %d, matched %d of %d; beta_obs: adjusted %d, matched %d of %d\n', r.genPass, r.obsPass);
    fprintf('R31c unadjusted max error beta_gen* %+.4f, beta_obs %+.4f; Fraser W realised gap %s\n', r.genMaxUnadj, r.obsMaxUnadj, mat2str(r.fraserWgap));
    disp(r.missesGen); disp(r.missesObsAdj); disp(r.missesObsMat);
end
save(fullfile(ROOT, 'results', 'N48', 'results', 'designRecsTempoNumbers_v001.mat'), 'T');  fprintf('Saved results/N48/results/designRecsTempoNumbers_v001.mat\n');
end
