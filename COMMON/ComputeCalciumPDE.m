
% --- Load global model parameters ----------------------------------------

addpath('../COMMON/');
mkdir('../DATA');
CommonParameters;

% --- Simulation call setup -----------------------------------------------

Script = '  ../COMMON/GenerateCalciumData.par';

if contains(computer, 'WIN')
    prog   = ['..\CALC\cwin6107x64.exe ', Script];
else
    prog   = ['../CALC/cmac6107x11 '    , Script];
end

% --- Output file names for each buffer condition -------------------------

fileControl  = '../DATA/CaControl';
fileEGTA2mM  = '../DATA/CaEGTA2mM';
fileEGTA10mM = '../DATA/CaEGTA10mM';
fileBAPTA1mM = '../DATA/CaBAPTA1mM';

% --- Output file names for Fig S5: non-inactivating Ca channel  ----------

fileEGTA2mM_ICa  = '../DATA/CaEGTA2mM_constICa';
fileEGTA10mM_ICa = '../DATA/CaEGTA10mM_constICa';
fileBAPTA1mM_ICa = '../DATA/CaBAPTA1mM_constICa';

%  Call CalC simulator with different buffer parameters and const ICa flag 

Y1 = ''; Y2 = ''; Y3 = ''; Y4 = ''; Y5 = ''; Y6 = ''; Y7 = '';

if ~isfile(fileControl)      [~, Y1] = system( [ prog, ' ', num2str([0,    0, 1]), ' ', fileControl  ]); end % Generate control data
if ~isfile(fileEGTA2mM)      [~, Y2] = system( [ prog, ' ', num2str([2e3,  0, 1]), ' ', fileEGTA2mM  ]); end % Generate Ca2+ with EGTA  = 2mM
if ~isfile(fileEGTA10mM)     [~, Y3] = system( [ prog, ' ', num2str([1e4,  0, 1]), ' ', fileEGTA10mM ]); end % Generate Ca2+ with EGTA  = 10mM
if ~isfile(fileBAPTA1mM)     [~, Y4] = system( [ prog, ' ', num2str([0,  1e3, 1]), ' ', fileBAPTA1mM ]); end % Generate Ca2+ with BAPTA = 1mM

if ~isfile(fileEGTA2mM_ICa)  [~, Y5] = system( [ prog, ' ', num2str([2e3,  0, 0]), ' ', fileEGTA2mM_ICa ]); end % EGTA 2mM, const ICa
if ~isfile(fileEGTA10mM_ICa) [~, Y6] = system( [ prog, ' ', num2str([1e4,  0, 0]), ' ', fileEGTA10mM_ICa]); end % EGTA 10mM, const ICa
if ~isfile(fileBAPTA1mM_ICa) [~, Y7] = system( [ prog, ' ', num2str([0,  1e3, 0]), ' ', fileBAPTA1mM_ICa]); end % BAPTA 1mM, const ICa

% --- Error message if CalC simulation failed -----------------------------

if numel(Y1) + numel(Y2) + numel(Y3) + numel(Y4) + numel(Y5) + numel(Y6) + numel(Y7) > 0
    fprintf(" %s\n %s\n %s\n %s\n %s\n %s\n %s\n\n", Y1, Y2, Y3, Y4, Y5, Y6, Y7);
    fprintf("Error computing [Ca2+] using CalC: check README instructions\n");
    return;
else
    fprintf("*** CalC data generated successfully ***\n");
end

% ================================ END ====================================
