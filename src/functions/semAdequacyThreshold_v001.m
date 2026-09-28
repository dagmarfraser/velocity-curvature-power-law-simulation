function [semAdequate, semMarginal, mdc] = semAdequacyThreshold_v001()
% semAdequacyThreshold_v001  Canonical SEM adequacy thresholds for this project.
%
% Single source of truth. Every script that classifies a coordinate, draws an
% adequacy contour, or labels an adequacy line must call this rather than
% hardcoding a decimal. Nineteen hardcoded sites across src/ were found on
% 2026-09-09 carrying three mutually inconsistent values (an exact MDC/2.77, a
% truncated 0.0108, and a rounded 0.011); this function exists so that never
% recurs.
%
% OUTPUTS
%   semAdequate  MDC/2.77 = 0.0108303249...  SEM below this gives 95% confidence
%                in detecting a change of MDC (Haley & Fragala-Pinkham 2006;
%                Beckerman et al. 2001). Pre-registration Sec 3.2.
%   semMarginal  MDC/1.5 = 0.02             Upper bound of the pre-registered
%                marginal band. semAdequate <= SEM < semMarginal gives roughly
%                80% power at MDC, or 95% confidence for effects >= 0.055.
%                SEM >= semMarginal is the inadequate band.
%                Pre-registration Sec 4.2 and Sec 8.2 deliverable 4.
%   mdc          0.03                       Minimal Detectable Change, the
%                |delta-beta| divergence reported between autistic and
%                non-autistic cohorts (Cook et al. 2026; Fourie et al. 2024).
%
% REPORTING NOTE. Report the criterion as SEM < MDC/2.77 and gloss the decimal
% to five places (0.01083). Do not report a rounded 0.011 while computing
% against MDC/2.77: a threshold is a stipulation, not an estimate, so it
% carries no sampling error and must be stated at the precision actually used.
%
% MATERIALITY, verified 2026-09-09. Across perCoordinateSEM_v2_001.mat
% (3,493,466 coordinate-pipeline rows) only 10,690 rows (0.31%) fall in the
% [MDC/2.77, 0.011) sliver where the exact and rounded constants disagree. All
% six v004 empirical centroids sit at SEM 0.003 to 0.0095 across all six
% pipelines, and Dhieb sits near 0.019, so no dataset-pipeline cell changes
% verdict under any of the three constants. The inconsistency is a hygiene
% defect, not a results defect.
%
% USAGE
%   semAdequate = semAdequacyThreshold_v001;
%   [semAdequate, semMarginal] = semAdequacyThreshold_v001;
%
% Fraser, D.S. (2026)  v001

MDC = 0.03;          % pre-registration Sec 3.2
Z_FACTOR = 2.77;     % 1.96 * sqrt(2), the 95% MDC multiplier
MARGINAL_FACTOR = 1.5;

mdc         = MDC;
semAdequate = MDC / Z_FACTOR;
semMarginal = MDC / MARGINAL_FACTOR;

end
