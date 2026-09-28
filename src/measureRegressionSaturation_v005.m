function measureRegressionSaturation_v005()
% measureRegressionSaturation_v005  Real dataset geometries, LC alpha.
%
% Each geometry uses its dataset-specific LC (template-subtracted) alpha.
% Maoz is included as a reference. This isolates the effect of geometry
% (semi-axes, f0, FS, M) on the saturation surface while holding the noise
% model (LC alpha, Xu coloured noise) constant to each dataset.
%
% Compare to v004 (Maoz geometry, LC alpha sweep) to separate geometry
% effects from noise-colour effects.
%
% Dataset geometry values from loop closure results (median across trials):
%   Cook CTRL:  a=57.5mm, b=27.1mm, f0=0.369Hz, FS=133Hz, LC alpha=4.86
%   Hickman:    a=57.0mm, b=26.5mm, f0=0.828Hz, FS=133Hz, LC alpha=5.5
%   Pilot:      a=30.8mm, b=13.8mm, f0=0.732Hz, FS=240Hz, LC alpha=4.4
%   Zarandi:    a=64.0mm, b=15.6mm, f0=0.647Hz, FS=100Hz, LC alpha=3.4
%                         (ecc=0.970; mean of maj=4.41 and min=2.43)
%
% The sigma=0 (noiseless) column in the output is the key Zarandi diagnostic:
%   - if beta_analytic tracks beta_gen cleanly -> CCC collapse in real Zarandi
%     data is the Wann 1988 metronome artefact, not a pipeline geometry failure
%   - if beta_analytic is flat or distorted -> SG pipeline breaks at this
%     sampling geometry (Zarandi ecc=0.970, 6.5 cycles in 10s at FS=100)
%
% Fraser, D.S. (2026)  v005

    srcDir = fileparts(mfilename('fullpath'));

    %% --- Geometries (10 seconds at each dataset FS) ---------------------
    G = {
        'Maoz',      50.0, 25.0, 1.000, 100, 10,  0.0,  NaN;
        'Cook CTRL', 57.5, 27.1, 0.369, 133, 10,  4.86, 5.76;
        'Hickman',   57.0, 26.5, 0.828, 133, 10,  5.5,  5.68;
        'Pilot',     30.8, 13.8, 0.732, 240, 10,  4.4,  2.21;
        'Zarandi',   64.0, 15.6, 0.647, 100, 10,  3.4,  4.44;
    };
    % Columns: name, a, b, f0, FS, nCycles, alpha, sigmaEmp

    geom = struct();
    for k = 1:size(G,1)
        geom(k).name     = G{k,1};
        geom(k).a        = G{k,2};
        geom(k).b        = G{k,3};
        geom(k).f0       = G{k,4};
        geom(k).FS       = G{k,5};
        geom(k).nCycles  = G{k,6};
        geom(k).alpha    = G{k,7};
        geom(k).sigmaEmp = G{k,8};
    end

    betaGenSweep = linspace(0, 0.75, 10);
    sigmaSweep   = [0, 1, 2, 5, 10, 20, 35, 50];
    N_REPS       = 20;
    outFile      = fullfile(srcDir, 'measureRegressionSaturation_v005.mat');

    saturationSweepEngine_v001(geom, betaGenSweep, sigmaSweep, N_REPS, outFile);
end
