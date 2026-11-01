function sets = dagAdjustmentSets_v001(G, exposure, outcome, candidates)
%DAGADJUSTMENTSETS_V001  Minimal backdoor adjustment sets by exhaustive search.
%   Returns a cell array of string column vectors. An empty cell array means no
%   admissible set exists among the candidates; a cell holding an empty string
%   array means the empty set suffices. Candidates that are descendants of the
%   exposure are excluded (backdoor criterion). Run from src/.
arguments
    G (1,1) digraph
    exposure (1,1) string
    outcome (1,1) string
    candidates (:,1) string
end
MAX_CANDIDATES = 16;                   % 2^16 subsets is the practical ceiling
if numel(candidates) > MAX_CANDIDATES
    error('dagAdjustmentSets_v001:tooMany', 'At most %d candidates (got %d).', ...
        MAX_CANDIDATES, numel(candidates));
end
names = string(G.Nodes.Name);
iX = find(names == exposure);
if isempty(iX)
    error('dagAdjustmentSets_v001:unknownNode', 'Unknown exposure "%s".', exposure);
end

desc = names(dfsearch(G, iX));         % exposure and all its descendants
cand = candidates(~ismember(candidates, [desc; outcome]));
Gb = rmedge(G, outedges(G, iX));       % backdoor graph: drop edges leaving the exposure

sets = {};
n = numel(cand);
for k = 0:n
    C = nchoosek(1:n, k);
    for r = 1:size(C, 1)
        Z = cand(C(r, :));
        Z = Z(:);
        if any(cellfun(@(s) all(ismember(s, Z)), sets))
            continue                   % superset of a smaller valid set
        end
        if dagDSeparated_v001(Gb, exposure, outcome, Z)
            sets{end+1} = Z;           %#ok<AGROW>
        end
    end
end
end
