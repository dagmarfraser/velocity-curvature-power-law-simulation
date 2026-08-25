function test_db_roundtrip_v001()
% test_db_roundtrip_v001  Beta values survive SQLite write-read with CAST workaround.
%
% Writes known float values (including non-binary-representable fractions
% like 1/3, 0.1, 0.285) to a temp SQLite DB, reads back both without and
% with CAST(col AS REAL) + double() workaround, and verifies:
%   (a) CAST workaround restores values within ROUNDTRIP_TOL
%   (b) Raw class of fetched values is reported (int64 truncation diagnostic)
%
% Temp DB is deleted via onCleanup regardless of test outcome.
%
% Requires: Database Toolbox
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

ROUNDTRIP_TOL = 1e-9;

% Production-relevant stress values: attractor ~0.285, 1/3, Zarandi alpha ~2.78
testVals = [0, 0.1, 1/3, 0.285, 0.5, 2.78, pi, -1/6, 0.001, 3.0];

fprintf('=== test_db_roundtrip_v001 ===\n');
fprintf('Testing %d float values, roundtrip tol=%.0e\n\n', numel(testVals), ROUNDTRIP_TOL);

%% ---- Create temp DB ----
tmpDB = fullfile(tempdir, sprintf('test_roundtrip_%s.db', ...
    strrep(datestr(now, 'yyyymmddHHMMSSFFF'), '.', '')));
cleanup = onCleanup(@() deleteTempDB(tmpDB));  %#ok<NASGU>

try
    conn = sqlite(tmpDB, 'create');
catch ME
    fprintf('FAIL: could not create temp SQLite DB: %s\n', ME.message);
    return;
end

execute(conn, ['CREATE TABLE IF NOT EXISTS rt (' ...
    'id INTEGER PRIMARY KEY, beta REAL, alpha_noise REAL)']);

for i = 1:numel(testVals)
    execute(conn, sprintf( ...
        'INSERT INTO rt (id, beta, alpha_noise) VALUES (%d, %.17g, %.17g)', ...
        i, testVals(i), testVals(i)));
end

%% ---- Read back without and with CAST ----
rawT  = fetch(conn, 'SELECT id, beta FROM rt ORDER BY id');
castT = fetch(conn, 'SELECT id, CAST(beta AS REAL) AS beta FROM rt ORDER BY id');
close(conn);

%% ---- Compare ----
fprintf('  %-6s  %-20s  %-20s  %-12s  %-10s  %s\n', ...
    'id', 'written', 'cast_read', 'abs_err', 'raw_class', 'verdict');
fprintf('  %s\n', repmat('-', 1, 82));

allPass    = true;
nTruncated = 0;

for i = 1:numel(testVals)
    written  = testVals(i);
    rawVal   = extractVal(rawT,  i, 'beta');
    castVal  = double(extractVal(castT, i, 'beta'));
    rawClass = class(rawVal);
    if isa(rawVal, 'int64'), nTruncated = nTruncated + 1; end

    err  = abs(castVal - written);
    pass = (err < ROUNDTRIP_TOL);
    if ~pass, allPass = false; end

    fprintf('  %-6d  %-20.17g  %-20.17g  %-12.2e  %-10s  %s\n', ...
        i, written, castVal, err, rawClass, verdict(pass));
end

fprintf('\n  CAST roundtrip all < %.0e: %s\n', ROUNDTRIP_TOL, verdict(allPass));
if nTruncated > 0
    fprintf('  Info: raw fetch returned int64 for %d/%d values (JDBC truncation confirmed)\n', ...
        nTruncated, numel(testVals));
end

fprintf('\n');
if allPass
    fprintf('PASS: all beta values survive SQLite round-trip with CAST workaround\n');
else
    fprintf('FAIL: one or more values exceeded roundtrip tolerance %.0e\n', ROUNDTRIP_TOL);
end

end

% =========================================================================
function val = extractVal(result, rowIdx, colName)
if istable(result)
    val = result.(colName)(rowIdx);
    if iscell(val), val = val{1}; end
else
    colNames = {'id','beta'};
    col = find(strcmp(colNames, colName), 1);
    val = result{rowIdx, col};
end
end

function deleteTempDB(path)
if isfile(path), try, delete(path); catch, end; end
end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
