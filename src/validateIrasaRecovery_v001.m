function validateIrasaRecovery_v001(opts)
% validateIrasaRecovery_v001  IRASA vs pmtm: synthetic validation + empirical agreement.
%
% TWO PARTS:
%
% PART 1: Synthetic validation (known-alpha signals)
%   Generates 1/f^alpha signals at alpha = 0:0.5:6 and at every project
%   sampling rate (60, 75, 100, 133, 240 Hz). For each combination, runs
%   both iraAlphaSigma_v001 and a direct pmtm log-log fit. Reports:
%   - Recovered vs true alpha for both methods
%   - |recovered - true| as function of alpha and fs
%   - Where IRASA and pmtm diverge from each other (bias ceiling)
%
% PART 2: Empirical agreement (real trial residuals)
%   Loads template-subtracted residuals from the five clinical datasets.
%   Runs IRASA and pmtm on the same residuals. Reports per-trial and
%   per-dataset:
%   - alpha_IRASA vs alpha_pmtm scatter (Brookshire-style figure)
%   - |alpha_IRASA - alpha_pmtm| distribution per dataset
%   - Agreement gate: what fraction of trials agree within 0.3?
%   - Whether agreement degrades at high alpha (Hickman datasets)
%
% RATIONALE:
%   Template subtraction removes dominant harmonics. If residuals are
%   genuinely fractal background, IRASA and plain pmtm should give the
%   same slope: the IRASA step (which was designed to remove oscillatory
%   peaks without knowing f0) is redundant when peaks are pre-removed.
%   Agreement between methods is therefore evidence that (a) template
%   subtraction was effective, and (b) the alpha estimate is trustworthy.
%   Disagreement flags residual oscillatory contamination or estimator
%   failure at extreme alpha.
%
%   Follows the validation logic of Brookshire (2022) who used multi-
%   method agreement as evidence for the aperiodic nature of broadband
%   spectral slopes.
%
% OUTPUT:
%   figures/iraValidation_synthetic_v001.png  : Part 1 recovery curves
%   figures/iraValidation_agreement_v001.png  : Part 2 Brookshire scatter
%   figures/iraValidation_agreementDist_v001.png: Part 2 gap distributions
%   src/iraValidationResults_v001.mat         : all numeric results
%
% USAGE:
%   validateIrasaRecovery_v001()               % full run (fast; ~5 min)
%   validateIrasaRecovery_v001(Part=1)         % synthetic only
%   validateIrasaRecovery_v001(Part=2)         % empirical only
%
% Fraser, D.S. (2026)  v001

    arguments
        opts.Part         (1,1) double  = 0      % 0=both, 1=synthetic, 2=empirical
        opts.FHighVals    (1,:) double  = [10 20] % upper bounds to compare (Hz)
        opts.NReps        (1,1) double  = 30     % synthetic replicates per condition
        opts.NTrialsCap   (1,1) double  = 200    % empirical: max trials per dataset
        opts.AgreeThr     (1,1) double  = 0.30   % |IRASA-pmtm| < threshold = agree
        opts.SaveFig      (1,1) logical = true
        opts.SaveMat      (1,1) logical = true
        opts.RngSeed      (1,1) double  = 1729
    end

    rng(opts.RngSeed, 'twister');

    srcDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(srcDir, 'functions')));
    figDir = fullfile(srcDir, '..', 'figures');
    if ~exist(figDir,'dir'), mkdir(figDir); end

    % IRASA parameters: identical to runLoopClosureFftnoise_v007
    IRA_FLOW = 1.0;
    IRA_HSET = 1.1:0.05:1.9;
    % fHigh = min(20, fs/2-1): computed per call

    fprintf('=== IRASA VALIDATION + BROOKSHIRE AGREEMENT ===\n\n');

    synResults = [];
    empResults = [];

    %% ================================================================
    %% PART 1: Synthetic validation
    %% ================================================================
    if opts.Part == 0 || opts.Part == 1
        fprintf('--- PART 1: Synthetic validation ---\n');

        alphaTrue = 0 : 0.5 : 6;
        fsVals    = [60, 75, 100, 133, 240];
        fHighVals = opts.FHighVals;
        nAlpha    = numel(alphaTrue);
        nFs       = numel(fsVals);
        nFH       = numel(fHighVals);
        nReps     = opts.NReps;
        TRIAL_DUR = 3.6;

        % Results: [nAlpha x nFs x nFH x nReps]
        iraRec  = NaN(nAlpha, nFs, nFH, nReps);
        pmtRec  = NaN(nAlpha, nFs, nFH, nReps);

        fprintf('  Conditions: %d alpha x %d fs x %d fHigh x %d reps\n', ...
            nAlpha, nFs, nFH, nReps);

        for fi = 1:nFs
            fs = fsVals(fi);
            N  = round(TRIAL_DUR * fs);
            for ai = 1:nAlpha
                a = alphaTrue(ai);
                for ri = 1:nReps
                    try
                        x = XuNoise_v002(N, a, 1.0, fs);
                    catch
                        x = spectralShapeWhite_local(N, a, fs);
                    end
                    x = x - mean(x);
                    for fhi = 1:nFH
                        fHi = min(fHighVals(fhi), fs/2 - 1);
                        if fHi <= IRA_FLOW, continue; end
                        [aI, ~] = iraAlphaSigma_v001(x(:), fs, IRA_FLOW, fHi, IRA_HSET);
                        aP      = pmtmAlpha_local(x(:), fs, IRA_FLOW, fHi);
                        iraRec(ai, fi, fhi, ri) = aI;
                        pmtRec(ai, fi, fhi, ri) = aP;
                    end
                end
            end
            fprintf('  fs=%d Hz done\n', fs);
        end

        % Summary over reps (collapse dim 4)
        iraMed  = median(iraRec,  4, 'omitnan');  % [nAlpha x nFs x nFH]
        iraSD   = std(iraRec,  0, 4, 'omitnan');
        pmtMed  = median(pmtRec, 4, 'omitnan');
        pmtSD   = std(pmtRec, 0, 4, 'omitnan');
        iraBias = iraMed - alphaTrue(:);           % broadcast over fs and fHigh
        pmtBias = pmtMed - alphaTrue(:);
        ipmGap  = iraMed - pmtMed;                % IRASA vs pmtm per condition

        % Sensitivity: difference between fHigh values (if 2 supplied)
        sensBias = [];
        if nFH >= 2
            sensBias = iraMed(:,:,2) - iraMed(:,:,1);  % fHigh(2) - fHigh(1)
            fprintf('\n  Band sensitivity (IRASA alpha change: fHigh=%d vs fHigh=%d Hz):\n', ...
                fHighVals(2), fHighVals(1));
            fprintf('  %-8s', 'alpha');
            for fi = 1:nFs, fprintf('  %6s', sprintf('%dHz',fsVals(fi))); end
            fprintf('\n%s\n', repmat('-',1,8+nFs*8));
            for ai = 1:nAlpha
                fprintf('  %-8.1f', alphaTrue(ai));
                for fi = 1:nFs
                    fprintf('  %6.3f', sensBias(ai,fi));
                end
                fprintf('\n');
            end
            maxSens = max(abs(sensBias(:)),[],'omitnan');
            fprintf('  Max |sensitivity| across all conditions: %.3f\n', maxSens);
            if maxSens < 0.3
                fprintf('  --> Band choice is NOT driving results (max effect < 0.3)\n');
            else
                fprintf('  --> Band choice IS materially affecting some conditions\n');
            end
        end

        % Bias ceiling per fs and fHigh
        for fhi = 1:nFH
            fprintf('\n  fHigh=%.0f Hz | alpha where |IRASA bias|>0.3:\n', fHighVals(fhi));
            for fi = 1:nFs
                bv = abs(iraBias(:,fi,fhi));
                th = find(bv > 0.3, 1);
                if isempty(th)
                    fprintf('    fs=%d: unbiased across all alpha\n', fsVals(fi));
                else
                    fprintf('    fs=%d: bias ceiling at alpha=%.1f\n', ...
                        fsVals(fi), alphaTrue(th));
                end
            end
        end

        synResults = struct('alphaTrue',alphaTrue,'fsVals',fsVals, ...
            'fHighVals',fHighVals, ...
            'iraRec',iraRec,'pmtRec',pmtRec, ...
            'iraMed',iraMed,'iraSD',iraSD,'pmtMed',pmtMed,'pmtSD',pmtSD, ...
            'iraBias',iraBias,'pmtBias',pmtBias,'ipmGap',ipmGap, ...
            'sensBias',sensBias);

        if opts.SaveFig
            plotSynthetic_local(synResults, figDir, alphaTrue, fsVals, fHighVals);
        end
    end

    %% ================================================================
    %% PART 2: Empirical agreement (Brookshire-style)
    %% ================================================================
    if opts.Part == 0 || opts.Part == 2
        fprintf('\n--- PART 2: Empirical Brookshire-style agreement ---\n');

        datasets = {
            'Pilot',        'noiseCharacterisation_pilot.mat',       240, 1/9.73, ...
                @() importDB_pilot_v001(fullfile(srcDir,'..','data','pilot'), ...
                    'Shapes',3,'Visits',[1 2],'Verbose',false);
            'Cook CTRL',    'noiseCharacterisation_cook.mat',         133, 0.248, ...
                @() importDB_cook_v002('Group','CTRL','Tasks',7,'Verbose',false);
            'Cook ASD',     'noiseCharacterisation_cookASD.mat',      133, 0.248, ...
                @() importDB_cook_v002('Group','ASD','Tasks',7,'Verbose',false);
            'Hickman PLAC', 'noiseCharacterisation_hickmanPLAC.mat',  133, 0.248, ...
                @() importDB_hickman_v003('Study',2,'Group','PLAC','Verbose',false);
            'Hickman HALO', 'noiseCharacterisation_hickmanHALO.mat',  133, 0.248, ...
                @() importDB_hickman_v003('Study',2,'Group','HALO','Verbose',false);
        };
        fHighVals = opts.FHighVals;
        nFH       = numel(fHighVals);
        nDS       = size(datasets,1);

        % Dataset colours (consistent with project palette)
        dsCols = [
            0.10 0.65 0.60;   % Pilot
            0.12 0.42 0.78;   % Cook CTRL
            0.85 0.18 0.18;   % Cook ASD
            0.08 0.58 0.35;   % Hick PLAC
            0.75 0.35 0.10;   % Hick HALO
        ];

        % Storage: empTbl(di, fhi) : dataset x fHigh
        empTbl = repmat(struct('name','','fs',0,'fHigh',0,'iraVec',[], ...
            'pmtVec',[],'gap',[],'agreeRate',NaN), nDS, nFH);
        allIRA = cell(nFH,1);
        allPMT = cell(nFH,1);
        allDS  = cell(nFH,1);
        for fhi = 1:nFH
            allIRA{fhi} = []; allPMT{fhi} = []; allDS{fhi} = [];
        end

        for di = 1:nDS
            dsName   = datasets{di,1};
            ncFile   = fullfile(srcDir, datasets{di,2});
            fs       = datasets{di,3};
            sigToMM  = datasets{di,4};
            importFn = datasets{di,5};

            fprintf('  Loading %s (fs=%d Hz)...\n', dsName, fs);
            bio      = load(ncFile,'bioResults').bioResults;
            trials   = importFn();
            trialIDs = string({trials.trialID}');
            selIdx   = randperm(height(bio), min(height(bio), opts.NTrialsCap));
            nSel     = numel(selIdx);

            % iraAll, pmtAll: [nSel x nFH]
            iraAll  = NaN(nSel, nFH);
            pmtAll  = NaN(nSel, nFH);
            % Major/minor axis alpha check (using project fHigh convention)
            % to verify: mean(alphaMaj, alphaMin) ≈ stored ira_alphaMean
            alphaMajVec  = NaN(nSel, 1);
            alphaMinVec  = NaN(nSel, 1);
            alphaStorVec = NaN(nSel, 1);   % from noiseCharacterisation

            for s = 1:nSel
                row = selIdx(s);
                tid = bio.trialID(row);
                m   = find(trialIDs == string(tid), 1);
                if isempty(m), continue; end

                tr   = trials(m);
                x240 = tr.x(:);  y240 = tr.y(:);

                f0 = estimateF0_local(double(x240), double(y240), fs);
                if ~isfinite(f0) || f0 <= 0, continue; end
                [rX, rY] = templateSubtract_local( ...
                    double(x240), double(y240), fs, f0, 4, 4);
                if numel(rX) < 20, continue; end

                xMM = double(x240)*sigToMM; xMM = xMM - mean(xMM);
                yMM = double(y240)*sigToMM; yMM = yMM - mean(yMM);
                [V, D] = eig(cov(xMM, yMM));
                [~, ord] = sort(diag(D),'descend'); V = V(:,ord);
                th   = atan2(V(2,1), V(1,1));
                rMaj = rX*cos(th) + rY*sin(th);
                rMin = -rX*sin(th) + rY*cos(th);   % minor axis residual

                % Both estimators at every fHigh : residuals computed once
                for fhi = 1:nFH
                    fHi = min(fHighVals(fhi), fs/2 - 1);
                    if fHi <= IRA_FLOW, continue; end
                    [aI, ~]       = iraAlphaSigma_v001(rMaj(:), fs, IRA_FLOW, fHi, IRA_HSET);
                    aP            = pmtmAlpha_local(rMaj(:), fs, IRA_FLOW, fHi);
                    iraAll(s,fhi) = aI;
                    pmtAll(s,fhi) = aP;
                end

                % Major/minor axis check at project fHigh (=min(20,fs/2-1))
                fHiProj = min(20.0, fs/2 - 1);
                [aMaj, ~] = iraAlphaSigma_v001(rMaj(:), fs, IRA_FLOW, fHiProj, IRA_HSET);
                [aMin, ~] = iraAlphaSigma_v001(rMin(:), fs, IRA_FLOW, fHiProj, IRA_HSET);
                alphaMajVec(s)  = aMaj;
                alphaMinVec(s)  = aMin;
                alphaStorVec(s) = bio.ira_alphaMean(row);   % stored project value
            end

            % Pack results per fHigh
            for fhi = 1:nFH
                fin = isfinite(iraAll(:,fhi)) & isfinite(pmtAll(:,fhi));
                gap = abs(iraAll(fin,fhi) - pmtAll(fin,fhi));
                empTbl(di,fhi).name      = dsName;
                empTbl(di,fhi).fs        = fs;
                empTbl(di,fhi).fHigh     = fHighVals(fhi);
                empTbl(di,fhi).iraVec    = iraAll(fin,fhi);
                empTbl(di,fhi).pmtVec    = pmtAll(fin,fhi);
                empTbl(di,fhi).gap       = gap;
                empTbl(di,fhi).agreeRate = mean(gap < opts.AgreeThr);
                allIRA{fhi} = [allIRA{fhi}; iraAll(fin,fhi)];
                allPMT{fhi} = [allPMT{fhi}; pmtAll(fin,fhi)];
                allDS{fhi}  = [allDS{fhi};  repmat(di,sum(fin),1)];
            end

            % Per-dataset sensitivity console report
            if nFH >= 2
                fin2   = isfinite(iraAll(:,1)) & isfinite(iraAll(:,2));
                diff12 = iraAll(fin2,2) - iraAll(fin2,1);
                fprintf('  %-14s fs=%d: med[%dHz]=%.3f  med[%dHz]=%.3f  |diff| med=%.3f max=%.3f\n', ...
                    dsName, fs, fHighVals(1), median(iraAll(fin2,1)), ...
                    fHighVals(2), median(iraAll(fin2,2)), ...
                    median(abs(diff12)), max(abs(diff12)));
            else
                fin = isfinite(iraAll(:,1));
                fprintf('  %-14s: N=%d  IRASA med=%.3f  pmtm med=%.3f\n', ...
                    dsName, sum(fin), median(iraAll(fin,1)), median(pmtAll(fin,1)));
            end

            % Major/minor axis decomposition report
            finAx = isfinite(alphaMajVec) & isfinite(alphaMinVec) & isfinite(alphaStorVec);
            if any(finAx)
                meanBothAxes = mean([alphaMajVec(finAx), alphaMinVec(finAx)], 2);
                fprintf('    Axis decomp (project fHigh=min(20,fs/2-1)):\n');
                fprintf('      alphaMaj med=%.3f  alphaMin med=%.3f  mean(maj+min) med=%.3f\n', ...
                    median(alphaMajVec(finAx)), median(alphaMinVec(finAx)), ...
                    median(meanBothAxes));
                fprintf('      stored ira_alphaMean med=%.3f  |diff vs mean(axes)| med=%.4f\n', ...
                    median(alphaStorVec(finAx)), ...
                    median(abs(meanBothAxes - alphaStorVec(finAx))));
            end
        end

        % Global sensitivity summary
        fprintf('\n  === Band sensitivity summary ===\n');
        for fhi = 1:nFH
            finAll = isfinite(allIRA{fhi});
            ipmGapAll = abs(allIRA{fhi}(finAll) - allPMT{fhi}(finAll));
            fprintf('  fHigh=%2d Hz: IRASA med=%.3f  pmtm med=%.3f  |IRASA-pmtm| med=%.3f  agree=%.0f%%\n', ...
                fHighVals(fhi), median(allIRA{fhi}(finAll)), ...
                median(allPMT{fhi}(finAll)), median(ipmGapAll), ...
                100*mean(ipmGapAll < opts.AgreeThr));
        end
        if nFH >= 2
            fin12 = isfinite(allIRA{1}) & isfinite(allIRA{2});
            r12   = corr(allIRA{1}(fin12), allIRA{2}(fin12));
            md12  = median(abs(allIRA{1}(fin12) - allIRA{2}(fin12)));
            mx12  = prctile(abs(allIRA{1}(fin12) - allIRA{2}(fin12)), 95);
            fprintf('\n  Cross-band IRASA (fHigh %d vs %d Hz): r=%.4f  |diff| med=%.3f  95pct=%.3f\n', ...
                fHighVals(1), fHighVals(2), r12, md12, mx12);
            if md12 < 0.3
                fprintf('  --> Band choice NOT driving empirical results (|diff| med < 0.3)\n');
            else
                fprintf('  --> Band choice IS materially affecting results (|diff| med >= 0.3)\n');
            end
        end

        empResults = struct('datasets',{empTbl},'allIRA',{allIRA},'allPMT',{allPMT}, ...
            'allDS',{allDS},'fHighVals',fHighVals,'agreeThresh',opts.AgreeThr);

        if opts.SaveFig
            for fhi = 1:nFH
                plotEmpirical_local(empResults, fhi, figDir, dsCols, opts.AgreeThr);
            end
            if nFH >= 2
                plotSensitivity_local(empResults, figDir, dsCols);
            end
        end
    end

    %% ---- Save ----------------------------------------------------------
    if opts.SaveMat
        matPath = fullfile(srcDir, 'iraValidationResults_v001.mat');
        save(matPath, 'synResults', 'empResults', '-v7.3');
        fprintf('\nMAT saved: %s\n', matPath);
    end

    fprintf('\n=== DONE ===\n');
    fprintf('Defensibility statement:\n');
    if ~isempty(empResults)
        fprintf('  IRASA-pmtm agreement per fitting band:\n');
        for fhi = 1:numel(empResults.fHighVals)
            allGap = vertcat(empResults.datasets(:,fhi).gap);
            fprintf('  fHigh=%2d Hz: |IRASA-pmtm| median=%.4f  95pct=%.4f  agree(%.1f)=%.0f%%\n', ...
                empResults.fHighVals(fhi), median(allGap), prctile(allGap,95), ...
                empResults.agreeThresh, 100*mean(allGap < empResults.agreeThresh));
        end
    end
end

%% ==========================================================================
function alpha = pmtmAlpha_local(x, fs, fLow, fHigh)
% Direct pmtm log-log alpha estimate (no resampling; no harmonic removal).
% TBP=4 matches iraAlphaSigma_v001. Used as the Brookshire comparison method.
    alpha = NaN;
    try
        [pxx, f] = pmtm(x(:), 4, [], fs);
        mask = f >= fLow & f <= fHigh & pxx > 0;
        if sum(mask) < 3, return; end
        p = polyfit(log10(f(mask)), log10(pxx(mask)), 1);
        alpha = -p(1);
    catch
    end
end

%% ==========================================================================
function x = spectralShapeWhite_local(N, alpha, fs)
% Fallback 1/f^alpha generator via spectral shaping (if XuNoise unavailable).
    w  = randn(N,1);
    W  = fft(w);
    f  = (0:N-1)' * fs / N;
    f(1) = f(2);   % avoid 0 Hz
    S  = f .^ (-alpha/2);
    S(1) = 0;
    x  = real(ifft(W .* S));
    x  = x / std(x);
end

%% ==========================================================================
function plotSynthetic_local(syn, figDir, alphaTrue, fsVals, fHighVals)
    nAlpha = numel(alphaTrue);
    nFs    = numel(fsVals);
    nFH    = numel(fHighVals);
    cols   = lines(nFs);

    % Use first fHigh value for the main 6-panel figure (default fHigh=20)
    fhi1 = nFH;  % last entry = largest fHigh (= 20 if [10 20])
    iraMed1 = syn.iraMed(:,:,fhi1);  iraSD1 = syn.iraSD(:,:,fhi1);
    pmtMed1 = syn.pmtMed(:,:,fhi1);  pmtSD1 = syn.pmtSD(:,:,fhi1);
    iraBias1 = syn.iraBias(:,:,fhi1); pmtBias1 = syn.pmtBias(:,:,fhi1);
    ipmGap1  = syn.ipmGap(:,:,fhi1);

    fig = figure('Units','centimeters','Position',[2 2 24 16],'Color','w');
    tl  = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    % Panel 1: IRASA recovery
    ax1 = nexttile(tl);
    hold(ax1,'on');
    plot(ax1,[0 6],[0 6],'k--','LineWidth',0.8,'DisplayName','ideal');
    for fi = 1:nFs
        errorbar(ax1, alphaTrue, iraMed1(:,fi), iraSD1(:,fi), ...
            'o-','Color',cols(fi,:),'MarkerFaceColor',cols(fi,:), ...
            'MarkerSize',5,'LineWidth',1.1, ...
            'DisplayName',sprintf('%d Hz',fsVals(fi)));
    end
    set(ax1,'Box','off','TickDir','out','FontSize',8);
    xlabel(ax1,'True \alpha'); ylabel(ax1,'IRASA recovered \alpha');
    title(ax1,sprintf('IRASA recovery (fHigh=%d Hz)',fHighVals(fhi1)));
    legend(ax1,'FontSize',7,'Box','off');

    % Panel 2: pmtm recovery
    ax2 = nexttile(tl);
    hold(ax2,'on');
    plot(ax2,[0 6],[0 6],'k--','LineWidth',0.8);
    for fi = 1:nFs
        errorbar(ax2, alphaTrue, pmtMed1(:,fi), pmtSD1(:,fi), ...
            's-','Color',cols(fi,:),'MarkerFaceColor',cols(fi,:), ...
            'MarkerSize',5,'LineWidth',1.1);
    end
    set(ax2,'Box','off','TickDir','out','FontSize',8);
    xlabel(ax2,'True \alpha'); ylabel(ax2,'pmtm log-log \alpha');
    title(ax2,'pmtm log-log recovery');

    % Panel 3: IRASA bias
    ax3 = nexttile(tl);
    hold(ax3,'on');
    plot(ax3,[0 6],[0 0],'k--','LineWidth',0.8);
    plot(ax3,[0 6],[ 0.3  0.3],'k:','LineWidth',0.6);
    plot(ax3,[0 6],[-0.3 -0.3],'k:','LineWidth',0.6);
    for fi = 1:nFs
        plot(ax3, alphaTrue, iraBias1(:,fi), 'o-', ...
            'Color',cols(fi,:),'MarkerFaceColor',cols(fi,:),'MarkerSize',5,'LineWidth',1.1);
    end
    set(ax3,'Box','off','TickDir','out','FontSize',8,'YLim',[-1.5 1.5]);
    xlabel(ax3,'True \alpha'); ylabel(ax3,'IRASA bias');
    title(ax3,'IRASA bias (dotted: \pm0.3)');

    % Panel 4: pmtm bias
    ax4 = nexttile(tl);
    hold(ax4,'on');
    plot(ax4,[0 6],[0 0],'k--','LineWidth',0.8);
    plot(ax4,[0 6],[ 0.3  0.3],'k:','LineWidth',0.6);
    plot(ax4,[0 6],[-0.3 -0.3],'k:','LineWidth',0.6);
    for fi = 1:nFs
        plot(ax4, alphaTrue, pmtBias1(:,fi), 's-', ...
            'Color',cols(fi,:),'MarkerFaceColor',cols(fi,:),'MarkerSize',5,'LineWidth',1.1);
    end
    set(ax4,'Box','off','TickDir','out','FontSize',8,'YLim',[-1.5 1.5]);
    xlabel(ax4,'True \alpha'); ylabel(ax4,'pmtm bias');
    title(ax4,'pmtm bias');

    % Panel 5: IRASA vs pmtm gap
    ax5 = nexttile(tl);
    hold(ax5,'on');
    plot(ax5,[0 6],[0 0],'k--','LineWidth',0.8);
    for fi = 1:nFs
        plot(ax5, alphaTrue, ipmGap1(:,fi), '^-', ...
            'Color',cols(fi,:),'MarkerFaceColor',cols(fi,:),'MarkerSize',5,'LineWidth',1.1, ...
            'DisplayName',sprintf('%d Hz',fsVals(fi)));
    end
    set(ax5,'Box','off','TickDir','out','FontSize',8);
    xlabel(ax5,'True \alpha'); ylabel(ax5,'IRASA - pmtm');
    title(ax5,'Estimator agreement (synthetic)');
    legend(ax5,'FontSize',7,'Box','off','Location','northwest');

    % Panel 6: band sensitivity (IRASA alpha change: fHigh2 - fHigh1)
    ax6 = nexttile(tl);
    hold(ax6,'on');
    if nFH >= 2 && ~isempty(syn.sensBias)
        plot(ax6,[0 6],[0 0],'k--','LineWidth',0.8);
        plot(ax6,[0 6],[ 0.3  0.3],'k:','LineWidth',0.6);
        plot(ax6,[0 6],[-0.3 -0.3],'k:','LineWidth',0.6);
        for fi = 1:nFs
            plot(ax6, alphaTrue, syn.sensBias(:,fi), 'o-', ...
                'Color',cols(fi,:),'MarkerFaceColor',cols(fi,:),'MarkerSize',5,'LineWidth',1.1, ...
                'DisplayName',sprintf('%d Hz',fsVals(fi)));
        end
        title(ax6,sprintf('Band sensitivity: IRASA fHigh %d vs %d Hz',fHighVals(1),fHighVals(2)));
    else
        text(ax6,0.5,0.5,'N/A (only one fHigh value)','Units','normalized','HorizontalAlignment','center');
        title(ax6,'Band sensitivity');
    end
    set(ax6,'Box','off','TickDir','out','FontSize',8,'YLim',[-1.0 1.0]);
    xlabel(ax6,'True \alpha'); ylabel(ax6,'\Delta\alpha (fHigh change)');
    legend(ax6,'FontSize',7,'Box','off');

    title(tl,'IRASA and pmtm synthetic validation: iraAlphaSigma\_v001', ...
        'FontSize',11,'FontWeight','bold');

    outP = fullfile(figDir,'iraValidation_synthetic_v001.png');
    exportgraphics(fig, outP,'Resolution',300);
    fprintf('  Fig saved: %s\n', outP);
end

%% ==========================================================================
function plotEmpirical_local(emp, fhi, figDir, dsCols, agrThr)
    nDS     = size(emp.datasets, 1);
    fHighV  = emp.fHighVals(fhi);
    suffix  = sprintf('fHigh%d', fHighV);

    %% Figure A: Brookshire scatter (IRASA vs pmtm per trial)
    figA = figure('Units','centimeters','Position',[2 2 18 14],'Color','w');
    ax   = axes(figA,'Units','normalized','Position',[0.12 0.12 0.74 0.76]);
    hold(ax,'on');

    maxV = 0;
    lgH  = gobjects(nDS,1);
    for di = 1:nDS
        d   = emp.datasets(di, fhi);
        col = dsCols(di,:);
        lgH(di) = scatter(ax, d.iraVec, d.pmtVec, 18, col, 'filled', ...
            'MarkerFaceAlpha',0.4,'DisplayName', ...
            sprintf('%s (N=%d, agree=%.0f%%)', d.name, numel(d.iraVec), ...
            100*d.agreeRate));
        maxV = max(maxV, max([d.iraVec; d.pmtVec],[],'omitnan'));
    end

    % Identity line and agreement bounds
    xL = [0 ceil(maxV)+0.5];
    plot(ax,xL,xL,'k-','LineWidth',1.2,'HandleVisibility','off');
    plot(ax,xL,xL+agrThr,'k:','LineWidth',0.7,'HandleVisibility','off');
    plot(ax,xL,xL-agrThr,'k:','LineWidth',0.7,'HandleVisibility','off');
    text(ax, xL(2), xL(2)+agrThr+0.1, sprintf('+%.1f',agrThr), ...
        'FontSize',7,'Color',[0.4 0.4 0.4]);

    % Correlation across all datasets
    finAll = isfinite(emp.allIRA{fhi}) & isfinite(emp.allPMT{fhi});
    r = corr(emp.allIRA{fhi}(finAll), emp.allPMT{fhi}(finAll));
    allGap = abs(emp.allIRA{fhi}(finAll) - emp.allPMT{fhi}(finAll));
    text(ax,0.05,0.92,sprintf('r=%.3f   |gap| med=%.3f', r, median(allGap)), ...
        'Units','normalized','FontSize',9,'FontWeight','bold');

    set(ax,'Box','off','TickDir','out','FontSize',9,'XLim',xL,'YLim',xL);
    xlabel(ax,'\alpha_{IRASA}','FontSize',11);
    ylabel(ax,'\alpha_{pmtm log-log}','FontSize',11);
    title(ax,sprintf('IRASA vs pmtm: empirical residuals, fHigh=%d Hz (Brookshire-style)', fHighV), ...
        'FontSize',10,'FontWeight','bold');
    subtitle(ax,sprintf('Dotted: \\pm%.1f band; both methods on same template-subtracted residuals',agrThr),'FontSize',8);
    legend(ax,lgH,'FontSize',7.5,'Box','off','Location','southeast');

    outA = fullfile(figDir,sprintf('iraValidation_agreement_%s_v001.png', suffix));
    exportgraphics(figA, outA,'Resolution',300);
    fprintf('  Fig A saved: %s\n', outA);

    %% Figure B: |gap| distributions per dataset
    figB = figure('Units','centimeters','Position',[2 16 18 8],'Color','w');
    axB  = axes(figB,'Units','normalized','Position',[0.11 0.15 0.85 0.72]);
    hold(axB,'on');

    for di = 1:nDS
        d = emp.datasets(di, fhi);
        if isempty(d.gap), continue; end
        % Empirical CDF of |gap|
        sorted = sort(d.gap);
        cdf    = (1:numel(sorted))' / numel(sorted);
        plot(axB, sorted, cdf*100, '-', 'Color',dsCols(di,:), ...
            'LineWidth',1.8,'DisplayName', ...
            sprintf('%s (med=%.2f)', d.name, median(d.gap)));
    end
    plot(axB,[agrThr agrThr],[0 100],'k:','LineWidth',0.9,'HandleVisibility','off');
    text(axB, agrThr+0.02, 15, sprintf('%.1f threshold',agrThr), ...
        'FontSize',8,'Color',[0.3 0.3 0.3]);

    set(axB,'Box','off','TickDir','out','FontSize',9,'YLim',[0 100],'XLim',[0 2]);
    xlabel(axB,'|alpha_{IRASA} - alpha_{pmtm}|','FontSize',10);
    ylabel(axB,'Cumulative %','FontSize',10);
    title(axB,sprintf('|IRASA - pmtm| distribution, fHigh=%d Hz',fHighV),'FontSize',10,'FontWeight','bold');
    legend(axB,'FontSize',7.5,'Box','off','Location','southeast');

    outB = fullfile(figDir,sprintf('iraValidation_agreementDist_%s_v001.png', suffix));
    exportgraphics(figB, outB,'Resolution',300);
    fprintf('  Fig B saved: %s\n', outB);
end

%% ==========================================================================
function f0 = estimateF0_local(x, y, fs)
    N = numel(x); nfft = 2^nextpow2(4*N);
    Xf = abs(fft(detrend(x(:),'linear'),nfft));
    Yf = abs(fft(detrend(y(:),'linear'),nfft));
    fAx = (0:nfft-1)'*fs/nfft;
    band = fAx > 0.1 & fAx < 5; idx = find(band);
    [~,pkX] = max(Xf(band)); [~,pkY] = max(Yf(band));
    f0x = fAx(idx(pkX)); f0y = fAx(idx(pkY));
    if abs(f0x-f0y)/max(f0x,0.01) < 0.2
        f0 = mean([f0x f0y]);
    elseif max(Xf(band)) > max(Yf(band))
        f0 = f0x;
    else
        f0 = f0y;
    end
end

%% ==========================================================================
%% ==========================================================================
function plotSensitivity_local(emp, figDir, dsCols)
% Figure C: IRASA alpha at fHigh(1) vs fHigh(2) per trial (sensitivity scatter).
% The key test: do the two band choices give the same answer?

    fHighVals = emp.fHighVals;
    nDS       = size(emp.datasets, 1);

    fig = figure('Units','centimeters','Position',[2 2 18 8],'Color','w');
    tl  = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    %% Panel 1: IRASA fHigh(1) vs fHigh(2) scatter
    ax1 = nexttile(tl);
    hold(ax1,'on');

    allA = emp.allIRA{1};
    allB = emp.allIRA{2};
    fin  = isfinite(allA) & isfinite(allB);
    allA = allA(fin); allB = allB(fin);
    dsIdx = emp.allDS{1}(fin);

    lgH = gobjects(nDS,1);
    for di = 1:nDS
        mask = dsIdx == di;
        lgH(di) = scatter(ax1, allA(mask), allB(mask), 12, dsCols(di,:), 'filled', ...
            'MarkerFaceAlpha',0.4,'DisplayName',emp.datasets(di,1).name);
    end

    xL = [0 max([allA;allB])*1.05+0.5];
    plot(ax1,xL,xL,'k-','LineWidth',1.0,'HandleVisibility','off');
    plot(ax1,xL,xL+0.3,'k:','LineWidth',0.7,'HandleVisibility','off');
    plot(ax1,xL,xL-0.3,'k:','LineWidth',0.7,'HandleVisibility','off');

    r    = corr(allA, allB);
    mdif = median(abs(allA-allB));
    text(ax1,0.05,0.90,sprintf('r=%.4f   |diff| med=%.3f',r,mdif), ...
        'Units','normalized','FontSize',8,'FontWeight','bold');

    set(ax1,'Box','off','TickDir','out','FontSize',9,'XLim',xL,'YLim',xL);
    xlabel(ax1,sprintf('\\alpha_{IRASA} (fHigh=%d Hz)',fHighVals(1)),'FontSize',10);
    ylabel(ax1,sprintf('\\alpha_{IRASA} (fHigh=%d Hz)',fHighVals(2)),'FontSize',10);
    title(ax1,'Band sensitivity: IRASA only','FontSize',9,'FontWeight','bold');
    legend(ax1,lgH,'FontSize',7,'Box','off','Location','southeast');

    %% Panel 2: |diff| CDF per dataset
    ax2 = nexttile(tl);
    hold(ax2,'on');
    for di = 1:nDS
        mask = emp.allDS{1}(fin) == di;
        if ~any(mask), continue; end
        d    = sort(abs(allA(mask) - allB(mask)));
        cdf  = (1:numel(d))'/numel(d)*100;
        plot(ax2, d, cdf, '-','Color',dsCols(di,:),'LineWidth',1.8, ...
            'DisplayName',sprintf('%s (med=%.2f)', emp.datasets(di,1).name, median(d)));
    end
    plot(ax2,[0.3 0.3],[0 100],'k:','LineWidth',0.9,'HandleVisibility','off');
    text(ax2,0.32,10,'0.3 threshold','FontSize',8,'Color',[0.3 0.3 0.3]);
    set(ax2,'Box','off','TickDir','out','FontSize',9,'YLim',[0 100],'XLim',[0 1.5]);
    xlabel(ax2,sprintf('|\\alpha_{%dHz} - \\alpha_{%dHz}|',fHighVals(1),fHighVals(2)),'FontSize',10);
    ylabel(ax2,'Cumulative %','FontSize',10);
    title(ax2,'Distribution of band difference','FontSize',9,'FontWeight','bold');
    legend(ax2,'FontSize',7,'Box','off','Location','southeast');

    title(tl,sprintf('Band sensitivity: fHigh=%d Hz vs fHigh=%d Hz (IRASA on real residuals)', ...
        fHighVals(1),fHighVals(2)), 'FontSize',10,'FontWeight','bold');

    outC = fullfile(figDir,'iraValidation_sensitivity_v001.png');
    exportgraphics(fig, outC,'Resolution',300);
    fprintf('  Fig C saved: %s\n', outC);
end
%% ==========================================================================
function [resX, resY, fitX, fitY] = templateSubtract_local(x, y, fs, f0, nH, nCW)
    N = numel(x);
    win = round(nCW/f0*fs); hop = round(win/2);
    win = min(win,N); win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    s = 1;
    while s+win-1 <= N
        e = s+win-1; tw = (0:win-1)'/fs;
        han = 0.5*(1-cos(2*pi*(0:win-1)'/(win-1)));
        D = [ones(win,1),tw];
        for h = 1:nH
            D(:,end+1) = cos(2*pi*h*f0*tw);  %#ok<AGROW>
            D(:,end+1) = sin(2*pi*h*f0*tw);  %#ok<AGROW>
        end
        Dw = D.*han;
        bX = Dw\(x(s:e).*han); bY = Dw\(y(s:e).*han);
        resX(s:e) = resX(s:e) + (x(s:e)-D*bX).*han;
        resY(s:e) = resY(s:e) + (y(s:e)-D*bY).*han;
        ww(s:e) = ww(s:e)+han; s = s+hop;
    end
    ok = ww>0;
    resX(ok) = resX(ok)./ww(ok);
    resY(ok) = resY(ok)./ww(ok);
    cl = max(round(win/4),1); cs = cl; ce = min(max(N-cl,cs+100),N);
    fitX = x(cs:ce)-resX(cs:ce); fitY = y(cs:ce)-resY(cs:ce);
    resX = resX(cs:ce); resY = resY(cs:ce);
end
