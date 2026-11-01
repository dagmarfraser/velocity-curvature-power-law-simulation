function tf = dagDSeparated_v001(G, X, Y, Z)
%DAGDSEPARATED_V001  True if node sets X and Y are d-separated given Z in digraph G.
%   Moral ancestral graph method. Fails loudly on unknown, overlapping or cyclic input.
%   Used by tempoDAG_v001.m and dagAdjustmentSets_v001.m (run from src/).
arguments
    G (1,1) digraph
    X (:,1) string
    Y (:,1) string
    Z (:,1) string = strings(0,1)
end
if ~isdag(G)
    error('dagDSeparated_v001:cyclic', 'Graph contains a cycle.');
end
names = string(G.Nodes.Name);
ix = nodeIdx(X);  iy = nodeIdx(Y);  iz = nodeIdx(Z);
if ~isempty(intersect(ix, iy)) || ~isempty(intersect(ix, iz)) || ~isempty(intersect(iy, iz))
    error('dagDSeparated_v001:overlap', 'X, Y and Z must be disjoint.');
end

Grev = flipedge(G);                    % ancestors of s = nodes reachable from s in Grev
seed = unique([ix; iy; iz]);
keep = seed;
for k = 1:numel(seed)
    keep = [keep; dfsearch(Grev, seed(k))];  %#ok<AGROW>
end
H = subgraph(G, unique(keep));

A = full(adjacency(H)) > 0;            % A(i,j): edge i -> j
M = A | A';
for j = 1:size(A, 2)                   % marry the parents of every node
    p = find(A(:, j));
    M(p, p) = true;
end
M(logical(eye(size(M)))) = false;

hn = string(H.Nodes.Name);
isZ = ismember(hn, Z);
M(isZ, :) = false;
M(:, isZ) = false;

cc = conncomp(graph(double(M)));
tf = isempty(intersect(cc(ismember(hn, X)), cc(ismember(hn, Y))));

    function idx = nodeIdx(S)
        idx = zeros(numel(S), 1);
        for q = 1:numel(S)
            f = find(names == S(q));
            if isempty(f)
                error('dagDSeparated_v001:unknownNode', 'Unknown node "%s".', S(q));
            end
            idx(q) = f;
        end
    end
end
