function sets = dagAdjustmentSets_v002(G, exposure, outcome, candidates)
%DAGADJUSTMENTSETS_V002  Minimal backdoor adjustment sets by exhaustive search.
%   As v001, with one fix (2026-10-24): v001 called nchoosek(1:n, k), which for n = 1 reads 1:1 as
%   the scalar 1 and returns the count, not the combinations, so with exactly one candidate the
%   empty set was never tested (found by tempoDAG_v002: the mediator reading returned {U} where no
%   adjustment is needed). Readings with two or more candidates are unaffected.
%   Returns a cell array of string column vectors. An empty cell array means no admissible set
%   exists among the candidates; a cell holding an empty string array means the empty set suffices.
%   Candidates that are descendants of the exposure are excluded (backdoor criterion). Run from src/.
arguments
    G (1,1) digraph
    exposure (1,1) string
    outcome (1,1) string
    candidates (:,1) string
end
MAX_CANDIDATES = 16;                   % 2^16 subsets is the practical ceiling
if numel(candidates) > MAX_CANDIDATES
    error('dagAdjustmentSets_v002:tooMany', 'At most %d candidates (got %d).', MAX_CANDIDATES, numel(candidates));
end
names = string(G.Nodes.Name);
iX = find(names == exposure);
if isempty(iX), error('dagAdjustmentSets_v002:unknownNode', 'Unknown exposure "%s".', exposure); end

desc = names(dfsearch(G, iX));         % exposure and all its descendants
cand = candidates(~ismember(candidates, [desc; outcome]));
Gb = rmedge(G, outedges(G, iX));       % backdoor graph: drop edges leaving the exposure

sets = {};
n = numel(cand);
for k = 0:n
    C = combos_local(n, k);
    for r = 1:size(C, 1)
        Z = cand(C(r, :));  Z = Z(:);
        if any(cellfun(@(s) all(ismember(s, Z)), sets)), continue, end   % superset of a smaller valid set
        if dagDSeparated_v001(Gb, exposure, outcome, Z), sets{end+1} = Z; end %#ok<AGROW>
    end
end
end

function C = combos_local(n, k)
% All k-subsets of 1:n as rows; one empty row for k = 0 (safe for n = 1, unlike nchoosek(1:n, k)).
if k == 0, C = zeros(1, 0);
elseif k == n, C = 1:n;
else, C = nchoosek(1:n, k);
end
end
