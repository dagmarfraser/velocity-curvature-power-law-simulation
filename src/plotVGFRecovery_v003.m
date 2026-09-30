% plotVGFRecovery_v003.m
%
% v003 (2026-09-30): the tempo conversion is the grid's exact K(1/3) = 185.0242 (gridKConv_v001, the
% path integral of kappa^(1/3) on the grid's nu=2 shape), replacing K_conv = 175.2636 =
% perimeter / mean(kappa^(-1/3)), which is 5.28% low (testGridKConv_v001). Every condition here is
% generated at beta_gen = 1/3, so no exponent choice enters: f0 = VGF / K(1/3) is exact for the
% generator, and f0_rec applies the same conversion to the recovered VGF (its 'f0 at beta = 1/3').
% The grid's VGF range 90.017-330.30 is 0.487-1.785 Hz (the old '~0.5-2 Hz' check was the K_conv
% error). VGF_gen and VGF_rec panels are unchanged; only f0 axes, the 1 Hz anchor (VGF 185.0, was
% 175.3) and the numbers file change. REGRESSION ANCHOR: with K = 175.2636 the numbers table must
% reproduce results/VGFRecovery_v002_numbers.txt (VGF and f0_rec to print precision).
% Two-panel figure per noise condition: f0_gen vs f0_rec per pipeline.
%
% v002 (2026-07-02): re-snapped Zarandi and Hickman PLAC conditions to
% characteriseBiologicalNoise_v004 canonical centroids (EMPIRICAL_DATASETS.md,
% Finding #6). v001's conditions were snapped to the superseded v003
% centroids (Zarandi alpha=2.53, Hickman PLAC alpha=4.60) -- see Finding #97
% for why v003 alpha was ~1-2 units too low. v001's numbers (f0_rec 1.63-2.84
% Hz for Zarandi, 0.993-1.023 Hz for Hickman PLAC) are stale and should not
% be cited; this version supersedes them. No other logic changed.
%
% Three conditions plotted (3 sub-figures, 2 panels each):
%   Cond 1 -- beta_gen=1/3, fs=120Hz, noiseless (alpha=0, sigma=0mm)
%   Cond 2 -- beta_gen=1/3, fs=120Hz, Zarandi noise (alpha=3.20, sigma=4.0mm)
%   Cond 3 -- beta_gen=1/3, fs=120Hz, Hickman PLAC noise (alpha=5.40, sigma=8.0mm)
%
% Per condition, two panels:
%   Panel A: VGF_gen (mm/s) vs VGF_rec -- shows pipeline ordering directly
%   Panel B: f0_gen (Hz) vs f0_rec (Hz) -- intuitive speed narrative
%            Identity line = perfect recovery; overrecovery = curve above diagonal
%
% Conversion: f0 = VGF / K(1/3), K(1/3) = gridKConv_v001(1/3) = 185.0242
%   (path integral of kappa^(1/3) ds on the grid's nu=2 shape; v002 used perimeter/mean(kappa^-1/3))
%   Grid VGF exp(4.5)-exp(5.8) = 90.017-330.30 -> f0 0.487-1.785 Hz at beta = 1/3 (checkGridTempoRange_v001)
%
% DB schema (powerlaw_debug_v058.db):
%   param_configs: generated_beta REAL, sampling_rate REAL, noise_type REAL (=alpha),
%                  noise_magnitude REAL (=sigma mm), vgf_value REAL,
%                  filter_type INT, regress_type INT
%   results:       vgf REAL (recovered), success INT
% Noise snapping (v004 centroids, alpha step=0.2, sigma grid non-uniform --
% see perCoordinateSEM_v2_001.mat coordTable for the exact grid):
%   Zarandi:      alpha=3.184->3.20, sigma=4.77->4.0mm  (confirmed via live
%                 nearest-neighbour snap against coordTable, 2026-07-02)
%   Hickman PLAC: alpha=5.337->5.40, sigma=7.17->8.0mm  (same check)
%
% Outputs:
%   figures/plotVGFRecovery_v003.png  (300 dpi, 3-column layout)
%   figures/plotVGFRecovery_v003.fig
%   results/VGFRecovery_v003_numbers.txt
%
% Fraser (2026) -- R9 figure / D3 desiderata; session 22 (v001); v004-centroid
% re-snap 2026-07-02 (v002)

function plotVGFRecovery_v003()

clc;
srcDir     = fileparts(mfilename('fullpath'));
projectDir = fileparts(srcDir);
figDir     = fullfile(projectDir, 'figures');
resDir     = fullfile(projectDir, 'results');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(resDir, 'dir'), mkdir(resDir); end

BETA_REF  = 1/3;
FS_TARGET = 120;

% Three noise conditions to plot
conditions = { ...
    'Noiseless',     0.0, 0.0; ...
    'Zarandi',       3.20, 4.0; ...
    'Hickman PLAC',  5.40, 8.0  ...
};

fprintf('==============================================\n');
fprintf('  plotVGFRecovery_v003  (v004 centroids)\n');
fprintf('==============================================\n\n');

%% 1. f0 CONVERSION FROM THE GRID'S OWN TEMPO CONSTANT
addpath(projectDir); addpath(srcDir); addpath(genpath(fullfile(srcDir, 'functions')));
K_RETIRED = 175.2636;  K13_PUB = 185.0242;
K_conv = gridKConv_v001(BETA_REF);          % [px^(2/3)]: VGF / K(1/3) = f0 [Hz]
if abs(K_conv - K13_PUB) > 1e-3
    error('plotVGF:K13', '%s', sprintf('K(1/3) %.4f does not reproduce %.4f', K_conv, K13_PUB));
end
vgf_check = exp([4.5, 5.8]);
f0_check  = vgf_check / K_conv;
fprintf('K(1/3) = %.4f (retired K_conv %.4f is %+.2f%% from it)\n', K_conv, K_RETIRED, 100 * (K_RETIRED / K_conv - 1));
fprintf('Grid VGF %.3f-%.3f -> f0 %.3f-%.3f Hz\n', vgf_check(1), vgf_check(2), f0_check(1), f0_check(2));

%% 2. LOCATE DATABASE
dbLocal = fullfile(resDir, 'powerlaw_debug_v058.db');
dbHPC   = ['/rds/projects/f/fraserds-mpo-evaluation/2026_prereg/' ...
           'velocity-curvature-power-law-simulation-main/' ...
           'velocity-curvature-power-law-simulation-main/results/' ...
           'powerlaw_debug_v058.db'];

dbMount = ['/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/' ...
           'velocity-curvature-power-law-simulation-main/' ...
           'velocity-curvature-power-law-simulation-main/results/powerlaw_debug_v058.db'];

if isfile(dbLocal)
    dbFile = dbLocal;
    fprintf('Using local database: %s\n', dbLocal);
elseif isfile(dbMount)
    dbFile = dbMount;
    fprintf('Using mounted RDS database: %s\n', dbMount);
elseif isfile(dbHPC)
    dbFile = dbHPC;
    fprintf('Using HPC database: %s\n', dbHPC);
else
    error('plotVGF:NoDB', '%s', ...
        'powerlaw_debug_v058.db not found. Ensure RDS is mounted or results/ is local.');
end

%% 3. QUERY PER CONDITION
% beta_gen = 1/3 stored as 10*(2/3)/20 = 0.33333...  use ABS tolerance
% noise_type = alpha value (stored as REAL); noise_magnitude = sigma in mm
BETA_TOL  = 0.005;
SIGMA_TOL = 0.5;
ALPHA_TOL = 0.15;

nCond = size(conditions, 1);
rawData = cell(nCond, 1);

for ci = 1:nCond
    condLabel = conditions{ci, 1};
    alpha_c   = conditions{ci, 2};
    sigma_c   = conditions{ci, 3};

    fprintf('Querying condition %d/%d: %s (alpha=%.1f, sigma=%.1f mm)...\n', ...
        ci, nCond, condLabel, alpha_c, sigma_c);

    sql = sprintf([ ...
        'SELECT CAST(pc.vgf_value AS REAL) AS vgf_gen, ' ...
        'pc.filter_type, pc.regress_type, ' ...
        'AVG(CAST(r.vgf AS REAL)) AS mean_vgf_rec, ' ...
        'COUNT(*) AS n ' ...
        'FROM param_configs pc ' ...
        'JOIN results r ON pc.config_id = r.config_id ' ...
        'WHERE r.success = 1 ' ...
        'AND r.vgf IS NOT NULL AND CAST(r.vgf AS REAL) > 0 ' ...
        'AND ABS(CAST(pc.generated_beta AS REAL) - %.6f) < %.4f ' ...
        'AND ABS(CAST(pc.sampling_rate   AS REAL) - %d)     < 1 ' ...
        'AND ABS(CAST(pc.noise_type      AS REAL) - %.2f)   < %.3f ' ...
        'AND ABS(CAST(pc.noise_magnitude AS REAL) - %.2f)   < %.3f ' ...
        'GROUP BY pc.vgf_value, pc.filter_type, pc.regress_type ' ...
        'ORDER BY vgf_gen, pc.filter_type, pc.regress_type'], ...
        BETA_REF, BETA_TOL, FS_TARGET, alpha_c, ALPHA_TOL, sigma_c, SIGMA_TOL);

    conn = sqlite(dbFile);
    raw  = fetch(conn, sql);
    close(conn);

    if ~istable(raw)
        raw = cell2table(raw, 'VariableNames', ...
            {'vgf_gen','filter_type','regress_type','mean_vgf_rec','n'});
    end
    raw.vgf_gen      = double(raw.vgf_gen);
    raw.mean_vgf_rec = double(raw.mean_vgf_rec);
    raw.filter_type  = double(raw.filter_type);
    raw.regress_type = double(raw.regress_type);

    nRows  = height(raw);
    nPipes = numel(unique(raw.filter_type .* 10 + raw.regress_type));
    fprintf('  %d rows (%d VGF levels x %d pipelines)\n', nRows, nRows/max(nPipes,1), nPipes);

    rawData{ci} = raw;
end

%% 4. PIPELINE SPEC
% filterType 2=BWFD(ref), 6=SG; regressType 3=OLS(ref), 4=LMLS, 5=IRLS
pipelines = { ...
    'BWFD-OLS',  2, 3, [0.10, 0.40, 0.80], '-',  1.8; ...
    'BWFD-LMLS', 2, 4, [0.10, 0.40, 0.80], '--', 1.5; ...
    'BWFD-IRLS', 2, 5, [0.10, 0.40, 0.80], ':',  1.5; ...
    'SG-OLS',    6, 3, [0.80, 0.20, 0.20], '-',  1.8; ...
    'SG-LMLS',   6, 4, [0.80, 0.20, 0.20], '--', 1.5; ...
    'SG-IRLS',   6, 5, [0.80, 0.20, 0.20], ':',  1.5  ...
};
nPipe = size(pipelines, 1);

% VGF grid from first condition with data
vgf_gen_vals = [];
for ci = 1:nCond
    if ~isempty(rawData{ci})
        vgf_gen_vals = sort(unique(rawData{ci}.vgf_gen));
        break;
    end
end
f0_gen_vals = vgf_gen_vals / K_conv;
nVGF = numel(vgf_gen_vals);

fprintf('\nVGF_gen: %.1f-%.1f mm/s  ->  f0: %.3f-%.3f Hz\n', ...
    min(vgf_gen_vals), max(vgf_gen_vals), min(f0_gen_vals), max(f0_gen_vals));

%% 5. EXTRACT RECOVERY CURVES PER PIPELINE x CONDITION
mean_vgf_rec = nan(nVGF, nPipe, nCond);
f0_rec_mat   = nan(nVGF, nPipe, nCond);

for ci = 1:nCond
    raw = rawData{ci};
    if isempty(raw), continue; end
    for p = 1:nPipe
        ft   = pipelines{p, 2};
        rt   = pipelines{p, 3};
        mask = raw.filter_type == ft & raw.regress_type == rt;
        sub  = raw(mask, :);
        [~, ia, ib] = intersect(vgf_gen_vals, sub.vgf_gen);
        mean_vgf_rec(ia, p, ci) = sub.mean_vgf_rec(ib);
        f0_rec_mat(ia, p, ci)   = sub.mean_vgf_rec(ib) / K_conv;
    end
end

%% 6. KEY NUMBERS FOR PAPER
F0_ANCHOR = 1.0;   % Hz -- 1 Hz as the f0 reference point (mid-grid)
VGF_ANCHOR = F0_ANCHOR * K_conv;

fprintf('\n=== KEY NUMBERS (f0_gen = 1 Hz, VGF_gen = %.1f mm/s) ===\n', VGF_ANCHOR);
numFile = fullfile(resDir, 'VGFRecovery_v003_numbers.txt');
fid = fopen(numFile, 'w');
fprintf(fid, 'plotVGFRecovery_v003 -- key numbers for R9 (v004 centroids)\n');
fprintf(fid, 'Generated: %s\n', datestr(now)); %#ok<DATST>
fprintf(fid, 'beta_gen=1/3, fs=120Hz\n');
fprintf(fid, 'K(1/3) = %.4f  (VGF / f0[Hz]; gridKConv_v001)\n\n', K_conv);

for ci = 1:nCond
    condLabel = conditions{ci, 1};
    fprintf(fid, '\n--- %s (alpha=%.1f, sigma=%.1f mm) ---\n', ...
        condLabel, conditions{ci,2}, conditions{ci,3});
    fprintf('\n%s:\n', condLabel);
    fprintf('%-14s  %12s  %12s  %12s\n', 'Pipeline','VGF_gen','VGF_rec','f0_rec(Hz)');
    fprintf(fid,'%-14s  %12s  %12s  %12s\n','Pipeline','VGF_gen','VGF_rec','f0_rec(Hz)');
    for p = 1:nPipe
        vr = interp1(vgf_gen_vals, mean_vgf_rec(:,p,ci), VGF_ANCHOR, 'linear', 'extrap');
        fr = vr / K_conv;
        fprintf('%-14s  %12.1f  %12.1f  %12.3f\n', pipelines{p,1}, VGF_ANCHOR, vr, fr);
        fprintf(fid,'%-14s  %12.1f  %12.1f  %12.3f\n', pipelines{p,1}, VGF_ANCHOR, vr, fr);
    end
end
fclose(fid);
fprintf('\nNumbers -> %s\n', numFile);

%% 6b. REGRESSION ANCHOR: retired K must reproduce v002's numbers file
oldFile = fullfile(resDir, 'VGFRecovery_v002_numbers.txt');
if ~isfile(oldFile), error('plotVGF:NoAnchor', '%s', 'FAILED PATH (anchor): ' + string(oldFile)); end
tok = regexp(fileread(oldFile), '^(?:BWFD|SG)-\w+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', 'tokens', 'lineanchors');
oldNum = cell2mat(cellfun(@(c) str2double(c), tok(:), 'UniformOutput', false));
newNum = zeros(0, 3);
for ci = 1:nCond
    for p = 1:nPipe
        vr = interp1(vgf_gen_vals, mean_vgf_rec(:,p,ci), 1.0 * K_RETIRED, 'linear', 'extrap');
        newNum(end+1, :) = [1.0 * K_RETIRED, vr, vr / K_RETIRED]; %#ok<AGROW>
    end
end
if ~isequal(size(oldNum), size(newNum)), error('plotVGF:AnchorSize', '%s', sprintf('anchor %s vs recomputed %s', mat2str(size(oldNum)), mat2str(size(newNum)))); end
if any(abs(oldNum(:,1) - newNum(:,1)) > 0.051) || any(abs(oldNum(:,2) - newNum(:,2)) > 0.051) || any(abs(oldNum(:,3) - newNum(:,3)) > 0.0006)
    error('plotVGF:Anchor', '%s', 'Retired-K recomputation does not reproduce VGFRecovery_v002_numbers.txt');
end
fprintf('REGRESSION ANCHOR passed: retired K reproduces v002 numbers (%d rows).\n', size(oldNum, 1));

%% 7. FIGURE  (3 columns = conditions, 2 rows = VGF / f0)
fig = figure('Units','centimeters','Position',[2 2 24 12],'Color','w', ...
    'Name','VGF Recovery by Noise Condition');

colTitles = conditions(:,1);

for ci = 1:nCond
    for row = 1:2   % row 1=VGF, row 2=f0
        axIdx = (row-1)*nCond + ci;
        ax = subplot(2, nCond, axIdx);
        hold(ax, 'on');

        if row == 1
            xVals = vgf_gen_vals;
            yMat  = mean_vgf_rec(:,:,ci);
            xLbl  = 'VGF_{gen} (mm/s^{2/3})';
            yLbl  = 'VGF_{rec} (mm/s^{2/3})';
        else
            xVals = f0_gen_vals;
            yMat  = f0_rec_mat(:,:,ci);
            xLbl  = 'f_0^{gen} (Hz)';
            yLbl  = 'f_0^{rec} (Hz)';
        end

        % Identity line
        lim = [min(xVals) max(xVals)];
        plot(ax, lim, lim, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        % Pipeline curves
        for p = 1:nPipe
            yp = yMat(:, p);
            if all(isnan(yp)), continue; end
            plot(ax, xVals, yp, ...
                pipelines{p,5}, 'Color', pipelines{p,4}, ...
                'LineWidth', pipelines{p,6}, 'DisplayName', pipelines{p,1});
        end

        xlabel(ax, xLbl, 'FontSize', 7);
        ylabel(ax, yLbl, 'FontSize', 7);
        axis(ax, 'square'); grid(ax, 'on'); box(ax, 'on');
        ax.FontSize = 7;

        if row == 1
            title(ax, colTitles{ci}, 'FontSize', 8, 'FontWeight', 'bold');
        end
        if ci == 1 && row == 1
            legend(ax, 'Location', 'northwest', 'FontSize', 6);
        end
    end
end

% Row labels
annotation('textbox',[0.01 0.74 0.05 0.1],'String','VGF', ...
    'EdgeColor','none','FontSize',8,'FontWeight','bold','Rotation',90);
annotation('textbox',[0.01 0.24 0.05 0.1],'String','f_0 (Hz)', ...
    'EdgeColor','none','FontSize',8,'FontWeight','bold','Rotation',90);

%% 8. SAVE
pngOut = fullfile(figDir, 'plotVGFRecovery_v003.png');
figOut = fullfile(figDir, 'plotVGFRecovery_v003.fig');
exportgraphics(fig, pngOut, 'Resolution', 300);
savefig(fig, figOut);
fprintf('Figure -> %s\n', pngOut);
fprintf('\nplotVGFRecovery_v003 COMPLETE.\n');

end
