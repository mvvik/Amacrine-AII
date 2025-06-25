clear
addpath('../COMMON/');
ComputeCalciumPDE;

% --- Load generated calcium profiles -------------------------------------

tArray                    = [1 20];      % Simulation durations of interest
[CaArrayEGTA2mM,  rArray] = ReadCaSliceAtTime(fileEGTA2mM_ICa,  tArray);
[CaArrayEGTA10mM, ~     ] = ReadCaSliceAtTime(fileEGTA10mM_ICa, tArray);
[CaArrayBAPTA1mM, ~     ] = ReadCaSliceAtTime(fileBAPTA1mM_ICa, tArray);

% =========================================================================
%                    Excess Buffer Approximation (EBA) for 1mM BAPTA
% =========================================================================

pA      = 5.182134;                    % Conversion: pA to µM·µm^3/ms
SCa     = ICa * pA;                    % Source strength (Ca2+ influx)
Cinf    = 0.05;                        % Baseline Ca2+ concentration (µM)
DC      = 0.22;                        % Ca2+ diffusion coefficient (µm^2/ms)
maxR    = 1000 * Rpatch;               % Max distance from Ca2+ channel (nm)

% --- BAPTA buffer parameters ---------------------------------------------
KD      = 0.176;                       % Dissociation constant (µM)
kminus  = 0.079;                       % Unbinding rate (ms^-1)
kplus   = kminus / KD;                 % Binding rate (µM^-1 ms^-1)
BT      = 1000;                        % Total BAPTA concentration (µM)
DB      = 0.2;                         % BAPTA diffusion coefficient (µm^2/ms)
Binf    = BT * KD / (KD + Cinf);       % Effective mobile buffer concentration

% --- EBA-related derived quantities --------------------------------------
lambda  = sqrt(DC / (Binf * kplus));   % Space constant (nm)
kappa   = KD * BT / (KD + Cinf)^2;     % Buffering capacity

% --- Free Ca2+ and EBA prediction profiles -------------------------------
CaFree      = SCa / (4 * pi * DC) ./ rArray + Cinf;
CaEBA_BAPTA = SCa / (4 * pi * DC) * exp(-rArray / lambda) ./ rArray + Cinf;

% =========================================================================
%                      Plot Settings and Style
% =========================================================================

rArray = 1000 * rArray;                % Convert µm → nm for plotting

CaLabel     = '[Ca^{2+}]';
CaLabelStr  = '[Ca^{2+}] (\muM)';

% --- Line styles ---------------------------------------------------------
lt1    = '-';  lt10   = '-';  ltFree = ':';  ltEBA  = ':';  ltFit1 = ':';  ltFit2 = ':';

% --- Colors --------------------------------------------------------------
clrEBA  = 0.7 * [1.0 1.0 0.2];
clrFit1 = 0.9 * [1.0 0.7 0.7];
clrFit2 = 0.9 * [0.6 1.0 1.0];
clr1    = 'm';
clr10   = 'b';
clrFree = 0.8 * [1 1 1];

% --- Line widths and font sizes ------------------------------------------
lw1 = 1.5; lw10 = 1.5; lwFit1 = 4; lwFit2 = 4; lwEBA = 4; lwFree = 3.5;
lfs = 11; tfs = 14;

% =========================================================================
%                      Plot 1: 2mM EGTA Case
% =========================================================================

figure(12); set(gcf, 'Position', [388 994 413 337]); set(gca, 'fontsize', lfs); hold off;

[x1, r, y1, ~] = TraceFit(rArray, CaArrayEGTA2mM(:,1));
[x2, r, y2, ~] = TraceFit(rArray, CaArrayEGTA2mM(:,2));

plot(rArray, CaArrayEGTA2mM(:,1), lt1,    'color', clr1,    'linewidth', lw1);     hold on;
plot(r,      y1,                  ltFit1, 'color', clrFit1, 'linewidth', lwFit1);
plot(rArray, CaArrayEGTA2mM(:,2), lt10,   'color', clr10,   'linewidth', lw10);
plot(r,      y2,                  ltFit2, 'color', clrFit2, 'linewidth', lwFit2);
plot(rArray, CaFree,              ltFree, 'color', clrFree, 'linewidth', lwFree);
plot(42,     102,                 '*',    'color', [1 1 1]);   % Dummy entry for spacing

legend( ...
    ['Simulated [Ca$^{2+}$] at t =$\bf ', num2str(tArray(1)), '$ms' ], ...
    sprintf('Best Fit: $%d \\ \\exp(-r / {\\bf%d}) / r + %.2f$', round(x1(1)), round(1/x1(2)), x1(3)), ...
    ['Simulated [Ca$^{2+}$] at t =$\bf ', num2str(tArray(2)), '$ms' ], ...
    sprintf('Best Fit: $%d \\ \\exp(-r / {\\bf%d}) / r + %.2f$', round(x2(1)), round(1/x2(2)), x2(3)), ...
    'Free $[Ca^{2+}] = I_{Ca} / (4\pi F D_{Ca} r) + Ca_{\infty}$', ...
    '               $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ = \ 375 \ / \ r + 0.05$', ...
    'location', 'best', 'fontsize', lfs, 'interpreter', 'latex');

set(gca, 'yscale', 'log');
xlabel('Distance from Ca^{2+} channel (nm)');
ylabel(CaLabelStr);
set(gca, 'TickLength', [0.04 0.02], 'box', 'off');
axis([1 maxR 0.05 400]);

% =========================================================================
%                      Plot 2: 10mM EGTA Case
% =========================================================================

figure(13); set(gcf, 'Position', [384 571 413 343]); set(gca, 'fontsize', lfs); hold off;

[x1, r, y1, ~] = TraceFit(rArray, CaArrayEGTA10mM(:,1));
[x2, r, y2, ~] = TraceFit(rArray, CaArrayEGTA10mM(:,2));

plot(rArray, CaArrayEGTA10mM(:,1), lt1,    'color', clr1,    'linewidth', lw1);     hold on;
plot(r,      y1,                   ltFit1, 'color', clrFit1, 'linewidth', lwFit1);
plot(rArray, CaArrayEGTA10mM(:,2), lt10,   'color', clr10,   'linewidth', lw10);
plot(r,      y2,                   ltFit2, 'color', clrFit2, 'linewidth', lwFit2);
plot(rArray, CaFree,               ltFree, 'color', clrFree, 'linewidth', lwFree);

legend( ...
    ['Simulated [Ca$^{2+}$] at t =$\bf ', num2str(tArray(1)), '$ms' ], ...
    sprintf('Best Fit: $%d \\ \\exp(-r / {\\bf%d}) / r + %.2f$', round(x1(1)), round(1/x1(2)), x1(3)), ...
    ['Simulated [Ca$^{2+}$] at t =$\bf ', num2str(tArray(2)), '$ms' ], ...
    sprintf('Best Fit: $%d \\ \\exp(-r / {\\bf%d}) / r + %.2f$', round(x2(1)), round(1/x2(2)), x2(3)), ...
    'Free $[Ca^{2+}] = 375 / r + 0.05$', ...
    'location', 'best', 'fontsize', lfs, 'interpreter', 'latex');

set(gca, 'yscale', 'log');
xlabel('Distance from Ca^{2+} channel (nm)');
ylabel(CaLabelStr);
set(gca, 'TickLength', [0.04 0.02], 'box', 'off');
axis([1 maxR 0.05 400]);

% =========================================================================
%                      Plot 3: 1mM BAPTA and EBA Comparison
% =========================================================================

figure(14); set(gcf, 'Position', [381 144 419 345]); set(gca, 'fontsize', lfs); hold off;

[x2, r, y2, ~] = TraceFit(rArray, CaArrayBAPTA1mM(:,2)); % Only fit 20ms data

plot(rArray, CaArrayBAPTA1mM(:,1), lt1,    'color', clr1,    'linewidth', lw1);      hold on;
plot(rArray, CaEBA_BAPTA,          ltEBA,  'color', clrEBA,  'linewidth', lwEBA); 
plot(42,             102,           '*',   'color', [1 1 1]);  % Dummy entries to space legend
plot(40,             100,           '*',   'color', [1 1 1]);
plot(rArray, CaArrayBAPTA1mM(:,2), lt10,   'color', clr10,   'linewidth', lw10);
plot(r,      y2,                   ltFit2, 'color', clrFit2, 'linewidth', lwFit2);
plot(rArray, CaFree,              ltFree, 'color', clrFree,  'linewidth', lwFree);

legend( ...
    ['Simulated [Ca$^{2+}$] at t =$\bf ', num2str(tArray(1)), '$ms' ], ...
    'EBA: $[Ca^{2+}] = A \exp(-r / \lambda) / r + Ca_{\infty}$', ...
    '$A = I_{Ca} / (4\pi F D_{Ca}) = 375\ \mu M\cdot nm$', ...
    '$\lambda = \sqrt{D_{Ca} / (B_{\infty}k_{on})} = {\bf 25}\ nm$', ...
    ['Simulated [Ca$^{2+}$] at t =$\bf ', num2str(tArray(2)), '$ms' ], ...
    sprintf('Best Fit: $%d \\ \\exp(-r / {\\bf%d}) / r + %.2f$', round(x2(1)), round(1/x2(2)), x2(3)), ...
    'Free $[Ca^{2+}] = 375 / r + 0.05$', ...
    'location', 'best', 'fontsize', lfs, 'interpreter', 'latex');

set(gca, 'yscale', 'log');
xlabel('Distance from Ca^{2+} channel (nm)');
ylabel(CaLabelStr);
set(gca, 'TickLength', [0.04 0.02], 'box', 'off');
axis([1 maxR 0.05 400]);

% ================================ END ====================================
