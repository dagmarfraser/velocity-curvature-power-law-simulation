function D = datasetOwnTempoExponent_v002(name)
% datasetOwnTempoExponent_v002  Single source for a dataset's own generator exponent, tempo and
% validity-domain flag, as used to place it on the grid's VGF axis (VGF_own = f0 * K(beta)).
%   D = datasetOwnTempoExponent_v002("Cook_CTRL")
%
%   name   one of Fraser, Cook_CTRL, Cook_ASD, Hickman_PLAC, Hickman_HALO, Zarandi, Dhieb.
%   D      struct: name, inDomain, nTrials, nNoF0, nExp, f0Median, f0Mean (ALL trials with a
%          usable f0), betaGenMed (median betaGenStarMed over trials with one; NaN if none),
%          betaDs (the exponent to use), betaRule (a string saying which rule gave it).
%
%   THE RULE (Dagmar, 2026-09-30, supersedes v001). Every dataset, in or out of the validity
%   domain, uses its own generator estimate: betaDs = median betaGenStarMed over trials that have
%   one. Placement on the VGF axis needs an exponent, not a valid generator, so the out-of-domain
%   sets (Zarandi, Dhieb; Finding #242) are treated identically but flagged: inDomain = false and
%   betaRule says its value is a model response, not an estimate of the participants' exponent.
%   v001 used a 1/3 reference there. Trials without an exponent are not imputed; f0 statistics use
%   all trials, so the tempo is not biased towards the fast, invertible subset (Finding #227).
%   Coverage is reported (nExp of nTrials): Dhieb has an exponent on only about a third of trials.
%
%   Fail Loud: errors if the corpus is missing, if an in-domain trial count differs from
%   Table 6a (2829, 94, 102, 359, 338), or if ANY dataset has no usable exponent.
% Reads: src/loopClosureResults_<name>_all_shaped_xu_v015.mat
% Fraser, D.S. (2026)  v002
    arguments
        name (1,1) string {mustBeMember(name, ["Fraser","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Zarandi","Dhieb"])}
    end
    ALL      = ["Fraser","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Zarandi","Dhieb"];
    IN_DOM   = [true true true true true false false];               % Finding #242
    N_6A     = [2829 94 102 359 338 NaN NaN];                        % Table 6a all-trial Ns
    k    = find(ALL == name);
    root = fileparts(fileparts(mfilename("fullpath")));
    f    = fullfile(root, "src", "loopClosureResults_" + name + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("datasetOwnTempo:NoCorpus", "%s", "FAILED PATH: " + f); end
    R = load(f, "results").results;
    n = numel(R);
    if ~isnan(N_6A(k)) && n ~= N_6A(k)
        error("datasetOwnTempo:Count", "%s", sprintf("%s: %d trials, Table 6a says %d", name, n, N_6A(k)));
    end
    f0 = arrayfun(@(r) double(r.f0), R(:));
    g  = arrayfun(@(r) double(r.betaGenStarMed), R(:));
    okF = isfinite(f0) & f0 > 0;
    D = struct("name", name, "inDomain", IN_DOM(k), "nTrials", n, "nNoF0", nnz(~okF), ...
        "nExp", nnz(okF & isfinite(g)), "f0Median", median(f0(okF)), "f0Mean", mean(f0(okF)), ...
        "betaGenMed", median(g(okF & isfinite(g))));
    if D.nNoF0 > 0
        warning("datasetOwnTempo:NoF0", "%s", sprintf("%s: %d trials without a usable f0 excluded", name, D.nNoF0));
    end
    if ~isfinite(D.betaGenMed) || D.betaGenMed < 0
        error("datasetOwnTempo:Beta", "%s", name + ": no usable dataset exponent");
    end
    D.betaDs = D.betaGenMed;
    if IN_DOM(k)
        D.betaRule = "own: median betaGenStarMed";
    else
        D.betaRule = "own: median betaGenStarMed (outside validity domain: model response, not a generator estimate)";
    end
end
