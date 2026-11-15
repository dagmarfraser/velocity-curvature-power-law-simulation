function S = applyLedgerSwap_v001(opts)
% APPLYLEDGERSWAP_V001  Rewrite a manuscript copy so each claims-ledger literal takes its value on a second result set.
%   For every row of a checkClaimsV006 ledger that is "ok" on the OLD results and not ok on the NEW results, find the
%   row's anchor in the target manuscript, find the OLD literal within the ledger's 2500-character window after it
%   (whitespace-tolerant, U+2212 read as '-'), and replace that one occurrence with the NEW literal. Rows whose OLD
%   literal cannot be found in the window are listed and left alone (Fail Loud: nothing is guessed).
%   S = applyLedgerSwap_v001(Old=<ledger .mat on N=45>, New=<ledger .mat on N=48>, Target=<manuscript .md>)
%   Writes Target in place (the caller makes the copy first) and returns the per-row status.   Fraser, D.S. (2026)  v001
arguments
    opts.Old (1,1) string
    opts.New (1,1) string
    opts.Target (1,1) string
end
O = load(opts.Old);  N = load(opts.New);
if size(O.C, 1) ~= size(N.C, 1) || ~isequal(string(O.C(:, 1)), string(N.C(:, 1))), error('swap:rows', '%s', 'ledger rows differ between Old and New'), end
raw = string(fileread(opts.Target));  nR = size(O.C, 1);
S = struct('row', {}, 'status', {}, 'old', {}, 'new', {});
for k = 1:nR
    if N.st(k) == "ok", continue, end
    r = struct('row', string(O.C{k, 1}), 'status', "", 'old', O.lit(k), 'new', N.lit(k));
    if O.st(k) ~= "ok", r.status = "OLD NOT OK"; S(end + 1) = r; continue, end %#ok<AGROW>
    anc = regexptranslate('escape', string(O.C{k, 3}));  anc = regexprep(anc, "\s+", "\\s+");
    a = regexp(raw, anc, 'once');
    if isempty(a), r.status = "ANCHOR LOST"; S(end + 1) = r; continue, end %#ok<AGROW>
    pat = regexptranslate('escape', O.lit(k));  pat = regexprep(pat, "(\\ )+|\s+", "\\s+");  pat = strrep(pat, "\-", "[-" + char(8722) + "]");   % regexptranslate escapes '-' as '\-'
    win = extractBetween(raw, a, min(strlength(raw), a + 3500));
    [s0, e0] = regexp(win, pat, 'once');
    if isempty(s0), r.status = "OLD NOT IN WINDOW"; S(end + 1) = r; continue, end %#ok<AGROW>
    p0 = a + s0 - 1;  p1 = a + e0 - 1;  newLit = N.lit(k);
    if contains(extractBetween(raw, p0, p1), char(8722))   % keep the manuscript's typographic minus: a sign '-' (after start, space, ( [ | or ; and before a digit or '.') becomes U+2212; range hyphens stay
        newLit = regexprep(newLit, '(^|\s|\(|\[|\||;)-(?=[0-9.])', ['$1' char(8722)]);   % alternation, not a class: MATLAB mis-reads '(' inside [...]
    end
    raw = extractBefore(raw, p0) + newLit + extractAfter(raw, p1);
    r.status = "SWAPPED";  S(end + 1) = r; %#ok<AGROW>
end
fid = fopen(opts.Target, 'w', 'n', 'UTF-8');  fwrite(fid, unicode2native(char(raw), 'UTF-8'));  fclose(fid);
st = string({S.status});
fprintf('%d rows to change: %d swapped, %d not swapped\n', numel(S), nnz(st == "SWAPPED"), nnz(st ~= "SWAPPED"));
for r = S(st ~= "SWAPPED"), fprintf('  %-22s %s\n', r.row, r.status); end
end
