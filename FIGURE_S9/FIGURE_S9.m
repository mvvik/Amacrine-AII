clear
addpath('../COMMON/');
ComputeCalciumPDE;

% --- Load time-resolved calcium profiles from files ---------------------
[~,      CaGrid1, tArray1] = ReadDataVsTime(fileControl);  % Control (no buffer)
[~,      CaGrid2, tArray2] = ReadDataVsTime(fileEGTA10mM); % EGTA 10 mM
[rArray, CaGrid3, tArray3] = ReadDataVsTime(fileBAPTA1mM); % BAPTA 1 mM

% --- Define radial range and simulation time ----------------------------
Rmin = 0.002;                      % in microns (2 nm)
Rmax = 0.080;                      % in microns (80 nm)
RR   = linspace(Rmin, Rmax, 100);  % evaluation points
nR   = numel(RR);
DT   = totalTime;                  % total integration time (from ComputeCalciumPDE)

% --- Allocate space for CaM binding results -----------------------------
Result  = zeros(3, nR);   % Fully bound CaM
ResultN = zeros(3, nR);   % N-lobe bound CaM
ResultC = zeros(3, nR);   % C-lobe bound CaM

% --- Compute CaM binding for all three conditions -----------------------
for nn = 1:nR
    R = RR(nn);
    [Result(1,nn), ResultN(1,nn), ResultC(1,nn)] = Get_CaM_binding(CaGrid1, tArray1, rArray, R, DT, @ode15s, 1e-5); % Control
    [Result(2,nn), ResultN(2,nn), ResultC(2,nn)] = Get_CaM_binding(CaGrid2, tArray2, rArray, R, DT, @ode15s, 1e-5); % EGTA
    [Result(3,nn), ResultN(3,nn), ResultC(3,nn)] = Get_CaM_binding(CaGrid3, tArray3, rArray, R, DT, @ode15s, 1e-5); % BAPTA
end

% --- Compute CaM binding ratios -----------------------------------------
EGTA10_to_Control    = Result(2,:)  ./ Result(1,:);
BAPTA1_to_Control    = Result(3,:)  ./ Result(1,:);
BAPTA1_to_EGTA10     = Result(3,:)  ./ Result(2,:);

EGTA10_to_Control_N  = ResultN(2,:) ./ ResultN(1,:);
BAPTA1_to_Control_N  = ResultN(3,:) ./ ResultN(1,:);
BAPTA1_to_EGTA10_N   = ResultN(3,:) ./ ResultN(2,:);

EGTA10_to_Control_C  = ResultC(2,:) ./ ResultC(1,:);
BAPTA1_to_Control_C  = ResultC(3,:) ./ ResultC(1,:);
BAPTA1_to_EGTA10_C   = ResultC(3,:) ./ ResultC(2,:);

% --- Plot CaM binding ratio: EGTA 10mM vs Control -----------------------
nRows = 1; nCols = 2;
figure;

subplot(nRows, nCols, 1); hold on;
plot(1000*RR, EGTA10_to_Control,   'r-', 'LineWidth', 2 );  % Fully bound
plot(1000*RR, EGTA10_to_Control_N, 'k-', 'LineWidth', 1 );  % N-lobe
plot(1000*RR, EGTA10_to_Control_C, 'b-', 'LineWidth', 1 );  % C-lobe

title(['\bfA  Bound CaM ratio: EGTA 10mM / Control', newline]);
legend('Fully bound CaM', 'Bound N-lobe', 'Bound C-lobe');
xlabel('Distance from channel (nm)');
ylabel('CaM binding ratio');
axis([0 1000*RR(end) 0 1]);
grid on;

% --- Plot CaM binding ratio: BAPTA 1mM vs Control -----------------------
subplot(nRows, nCols, 2); hold on;
plot(1000*RR, BAPTA1_to_Control,   'r-', 'LineWidth', 2 );  % Fully bound
plot(1000*RR, BAPTA1_to_Control_N, 'k-', 'LineWidth', 1 );  % N-lobe
plot(1000*RR, BAPTA1_to_Control_C, 'b-', 'LineWidth', 1 );  % C-lobe

title(['\bfB  Bound CaM ratio: BAPTA 1mM / Control', newline]);
legend('Fully bound CaM', 'Bound N-lobe', 'Bound C-lobe');
xlabel('Distance from channel (nm)');
ylabel('CaM binding ratio');
axis([0 1000*RR(end) 0 1]);
grid on;
