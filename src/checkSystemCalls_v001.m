function checkSystemCalls_v001()
%CHECKSYSTEMCALLS_V001  Can this MATLAB session make the system() calls buildN48ShadowTree_RDS_v001 needs? (2026-11-04)
%   Press Run. Touches only a temporary folder, which it removes. Prints PASS/FAIL per call and an overall verdict.
%   Needed: system() itself, ln -s (the data/ link), and a link that MATLAB can read through. md5sum is optional (logging).
d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
r = struct('name', {}, 'ok', {}, 'out', {});
[st, out] = system('hostname');                                   r(end+1) = struct('name', "system('hostname')", 'ok', st == 0, 'out', strtrim(out));
mkdir(fullfile(d, 'real'));  writelines("hello", fullfile(d, 'real', 'f.txt'));
[st, out] = system(sprintf('ln -s "%s" "%s"', fullfile(d, 'real'), fullfile(d, 'link')));
r(end+1) = struct('name', "ln -s", 'ok', st == 0, 'out', strtrim(out));
okRead = false;  try, okRead = strtrim(string(fileread(fullfile(d, 'link', 'f.txt')))) == "hello"; catch, end
r(end+1) = struct('name', "read through the link", 'ok', okRead, 'out', "");
[st, out] = system(sprintf('md5sum "%s"', fullfile(d, 'real', 'f.txt')));
r(end+1) = struct('name', "md5sum (optional)", 'ok', st == 0, 'out', strtrim(out));
for k = 1:numel(r), fprintf('%-24s %s  %s\n', r(k).name, string(ifelse_local(r(k).ok, "PASS", "FAIL")), r(k).out); end
if all([r(1:3).ok]), fprintf('=== READY: buildN48ShadowTree_RDS_v001 can run here ===\n');
else, fprintf('=== NOT READY: tell Claude which line failed ===\n'); end
end
function v = ifelse_local(c, a, b), if c, v = a; else, v = b; end, end
