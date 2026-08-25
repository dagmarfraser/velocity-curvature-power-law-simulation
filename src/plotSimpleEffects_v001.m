function plotSimpleEffects_v001()
% plotSimpleEffects_v001  Fig 5: Simple effects forest plot with noise-conditional breakdown
%
% Parses simpleEffects_log_20260525_093605.txt from src/ or ../results/.
% Panel A: Grand-mean simple effects -- six S8.1 contrasts as horizontal
%          forest plot with 95% CIs, coloured by EXPECTED/UNEXPECTED/NO_PRED.
% Panel B: Noise-conditional heatmap -- delta per (dataset x contrast),
%          annotated with E/U/NP flag.
%
% Output: ../figures/plotSimpleEffects_v001_forest.png
%         ../figures/plotSimpleEffects_v001_heatmap.png
%
% Run from src/. Log file searched in src/ then ../results/.
% Fraser, Di Luca, Cook (2026) -- PreReg simulation study.

    logFile = locateLog();

    [gmContrasts, gmDelta, gmLo, gmHi, gmD, gmFlag] = parseS3(logFile);
    [datasets, s4Contrasts, s4Delta, ~, s4Flag]      = parseS4(logFile);

    %% Panel A: Grand-mean forest plot
    colExp   = [0.13 0.55 0.13];
    colUnexp = [0.80 0.20 0.20];
    colNP    = [0.50 0.50 0.50];

    nC = numel(gmContrasts);
    colours = zeros(nC, 3);
    for k = 1:nC
        switch gmFlag{k}
            case 'EXPECTED',   colours(k,:) = colExp;
            case 'UNEXPECTED', colours(k,:) = colUnexp;
            otherwise,         colours(k,:) = colNP;
        end
    end

    figA = figure('Units', 'centimeters', 'Position', [2 2 22 10], 'Color', 'w');
    axA  = axes(figA);
    hold(axA, 'on');

    % Draw ascending y (k=1 at bottom); contrast 1 appears at top via flipped labels
    for k = 1:nC
        y = nC - k + 1;
        line(axA, [gmLo(k) gmHi(k)], [y y], 'Color', colours(k,:), 'LineWidth', 2.5);
        scatter(axA, gmDelta(k), y, 50, colours(k,:), 'filled');
        text(axA, gmHi(k) + 0.0003, y, ...
             sprintf(' d=%.2f  %s', gmD(k), gmFlag{k}), ...
             'FontSize', 7, 'Color', colours(k,:), 'VerticalAlignment', 'middle', ...
             'Interpreter', 'none');
    end
    xline(axA, 0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.4);

    yticks(axA, 1:nC);
    yticklabels(axA, flip(gmContrasts));
    axA.TickLabelInterpreter = 'none';
    axA.FontSize = 9;
    xlabel(axA, '\Delta\beta_{bias}  (positive = first pipeline has more bias)', 'FontSize', 9);
    title(axA, 'Grand-mean simple effects (S8.1)  --  N=17,458,535', ...
          'FontSize', 10, 'FontWeight', 'bold');

    hE = line(axA, nan, nan, 'Color', colExp,   'LineWidth', 2);
    hU = line(axA, nan, nan, 'Color', colUnexp, 'LineWidth', 2);
    hP = line(axA, nan, nan, 'Color', colNP,    'LineWidth', 2);
    legend(axA, [hE hU hP], {'Expected', 'Unexpected', 'No prediction'}, ...
           'Location', 'southeast', 'FontSize', 8);
    box(axA, 'off');
    xlim(axA, [min(gmLo) - 0.003, max(gmHi) + 0.010]);

    outA = fullfile(fileparts(mfilename('fullpath')), '..', 'figures', ...
                    'plotSimpleEffects_v001_forest.png');
    exportgraphics(figA, outA, 'Resolution', 150);
    fprintf('Saved: %s\n', outA);
    close(figA);

    %% Panel B: Noise-conditional heatmap
    nD = numel(datasets);
    nK = numel(s4Contrasts);

    deltaM = nan(nD, nK);
    flagM  = repmat({''}, nD, nK);
    for di = 1:nD
        for ki = 1:nK
            idx = strcmp(s4Delta.dataset, datasets{di}) & ...
                  strcmp(s4Delta.contrast, s4Contrasts{ki});
            if any(idx)
                i1 = find(idx, 1);
                deltaM(di, ki) = s4Delta.delta(i1);
                flagM{di, ki}  = s4Flag{i1};
            end
        end
    end

    figB = figure('Units', 'centimeters', 'Position', [2 2 24 10], 'Color', 'w');
    axB  = axes(figB);
    imagesc(axB, deltaM);
    colormap(axB, redblue(256));
    clLim = max(abs(deltaM(:)), [], 'omitnan');
    if clLim > 0, clim(axB, [-clLim clLim]); end
    cb = colorbar(axB);
    cb.Label.String = '\Delta\beta_{bias}';
    cb.FontSize = 8;

    for di = 1:nD
        for ki = 1:nK
            if ~isempty(flagM{di,ki})
                sf = strrep(strrep(strrep(flagM{di,ki}, 'UNEXPECTED', 'U'), ...
                                   'EXPECTED', 'E'), 'NO_PRED', 'NP');
                text(axB, ki, di, sf, 'HorizontalAlignment', 'center', ...
                     'FontSize', 7, 'Color', 'k', 'FontWeight', 'bold', ...
                     'Interpreter', 'none');
            end
        end
    end

    xticks(axB, 1:nK); xticklabels(axB, s4Contrasts);
    yticks(axB, 1:nD); yticklabels(axB, datasets);
    axB.XTickLabelRotation = 25;
    axB.TickLabelInterpreter = 'none';
    axB.FontSize = 8;
    title(axB, 'Noise-conditional simple effects  (E=Expected  U=Unexpected  NP=No prediction)', ...
          'FontSize', 9, 'FontWeight', 'bold');
    box(axB, 'off');

    outB = fullfile(fileparts(mfilename('fullpath')), '..', 'figures', ...
                    'plotSimpleEffects_v001_heatmap.png');
    exportgraphics(figB, outB, 'Resolution', 150);
    fprintf('Saved: %s\n', outB);
    close(figB);
end

%% -- Helpers
function logFile = locateLog()
    base = fileparts(mfilename('fullpath'));
    candidates = { ...
        fullfile(base, 'simpleEffects_log_20260525_093605.txt'), ...
        fullfile(base, '..', 'results', 'simpleEffects_L9_summary.txt') };
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            logFile = candidates{k};
            return;
        end
    end
    error('plotSimpleEffects_v001:missingLog', '%s', ...
        ['Could not locate simpleEffects log. ' ...
         'Expected in src/ or ../results/. Copy from RDS if needed.']);
end

function [contrasts, delta, lo, hi, d, flag] = parseS3(logFile)
    lines     = readlines(logFile);
    contrasts = {}; delta = []; lo = []; hi = []; d = []; flag = {};
    inS3 = false;
    for k = 1:numel(lines)
        ln = char(lines(k));
        if contains(ln, 'S3:'), inS3 = true; end
        if inS3 && contains(ln, 'S4:'), break; end
        if ~inS3, continue; end
        tok = regexp(ln, ...
            '^\s+(.+?)\s{2,}([+-]?\d+\.\d+)\s+\[([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\]\s+([+-]?\d+\.\d+)\s+(\S+)', ...
            'tokens', 'once');
        if isempty(tok), continue; end
        contrasts{end+1} = strtrim(tok{1}); %#ok<AGROW>
        delta(end+1)     = str2double(tok{2}); %#ok<AGROW>
        lo(end+1)        = str2double(tok{3}); %#ok<AGROW>
        hi(end+1)        = str2double(tok{4}); %#ok<AGROW>
        d(end+1)         = str2double(tok{5}); %#ok<AGROW>
        flag{end+1}      = tok{6};             %#ok<AGROW>
    end
end

function [datasets, contrasts, deltaStruct, d, flag] = parseS4(logFile)
    lines     = readlines(logFile);
    datasets  = {}; contrasts = {};
    dsVec = {}; ctVec = {}; dlVec = []; dVec = []; flVec = {};
    inS4 = false; curDs = '';
    for k = 1:numel(lines)
        ln = char(lines(k));
        if contains(ln, 'S4:'), inS4 = true; end
        if inS4 && contains(ln, 'S5:'), break; end
        if ~inS4, continue; end
        dsHdr = regexp(ln, '^\s{2}(\S+.*?)\s+\(alpha', 'tokens', 'once');
        if ~isempty(dsHdr)
            curDs = strtrim(dsHdr{1});
            if ~ismember(curDs, datasets), datasets{end+1} = curDs; end %#ok<AGROW>
            continue;
        end
        tok = regexp(ln, ...
            '^\s+\S.*?\s{2,}(.+?)\s{2,}([+-]?\d+\.\d+)\s+\[([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\]\s+([+-]?\d+\.\d+)\s+(\S+)', ...
            'tokens', 'once');
        if isempty(tok) || isempty(curDs), continue; end
        ct = strtrim(tok{1});
        if ~ismember(ct, contrasts), contrasts{end+1} = ct; end %#ok<AGROW>
        dsVec{end+1} = curDs;              %#ok<AGROW>
        ctVec{end+1} = ct;                 %#ok<AGROW>
        dlVec(end+1) = str2double(tok{2}); %#ok<AGROW>
        dVec(end+1)  = str2double(tok{5}); %#ok<AGROW>
        flVec{end+1} = tok{6};             %#ok<AGROW>
    end
    deltaStruct.dataset  = dsVec;
    deltaStruct.contrast = ctVec;
    deltaStruct.delta    = dlVec;
    d    = dVec;
    flag = flVec;
end

function C = redblue(n)
    if nargin < 1, n = 256; end
    half = ceil(n/2);
    r1 = linspace(0.17, 1, half)';   g1 = linspace(0.51, 1, half)';   b1 = linspace(0.73, 1, half)';
    r2 = linspace(1, 0.84, n-half)'; g2 = linspace(1, 0.19, n-half)'; b2 = linspace(1, 0.15, n-half)';
    C  = [[r1; r2], [g1; g2], [b1; b2]];
end
