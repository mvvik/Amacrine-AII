
% --- Load global model parameters ----------------------------------------

addpath('../COMMON/');
if ~isfolder('../DATA/'), mkdir('../DATA/'); end
CommonParameters;

% --- Simulation call setup -----------------------------------------------

if contains(computer, 'WIN')
     prog = '..\CALC\cwin6107x64.exe ';
elseif contains(computer('arch'), 'maca64')
     system('chmod +x ../CALC/cmac6107xM1');
     prog = '../CALC/cmac6107xM1 ';
elseif contains(computer('arch'), 'maci64')
     system('chmod +x ../CALC/cmac6107x86');
     prog = '../CALC/cmac6107x86 ';
else
     fprintf('Do not recognize this architecture: %s\n', computer('arch'));
     fprintf('Compile CalC first: github.com/mvvik/CalC-simple-buffer');
     return;
end

prog = [ prog, ' ../COMMON/GenerateCalciumData.par'];

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

try
    if ~isfile(fileControl)      [~, Y1] = system( [ prog, ' ', num2str([0,    0, 1]), ' ', fileControl  ]); end % Generate control data
    if ~isfile(fileEGTA2mM)      [~, Y2] = system( [ prog, ' ', num2str([2e3,  0, 1]), ' ', fileEGTA2mM  ]); end % Generate Ca2+ with EGTA  = 2mM
    if ~isfile(fileEGTA10mM)     [~, Y3] = system( [ prog, ' ', num2str([1e4,  0, 1]), ' ', fileEGTA10mM ]); end % Generate Ca2+ with EGTA  = 10mM
    if ~isfile(fileBAPTA1mM)     [~, Y4] = system( [ prog, ' ', num2str([0,  1e3, 1]), ' ', fileBAPTA1mM ]); end % Generate Ca2+ with BAPTA = 1mM
    
    if ~isfile(fileEGTA2mM_ICa)  [~, Y5] = system( [ prog, ' ', num2str([2e3,  0, 0]), ' ', fileEGTA2mM_ICa ]); end % EGTA 2mM, const ICa
    if ~isfile(fileEGTA10mM_ICa) [~, Y6] = system( [ prog, ' ', num2str([1e4,  0, 0]), ' ', fileEGTA10mM_ICa]); end % EGTA 10mM, const ICa
    if ~isfile(fileBAPTA1mM_ICa) [~, Y7] = system( [ prog, ' ', num2str([0,  1e3, 0]), ' ', fileBAPTA1mM_ICa]); end % BAPTA 1mM, const ICa
catch ERR
    if contains(computer('arch'), 'maca')
        fprintf('Give permissions to the CALC exectuable in settings\n');
        system('open \"x-apple.systempreferences:com.apple.preference.securoty?Privacy\"');
    end
    fprintf('Error in executing CalC code: %s\n%s\n', prog, ERR);
    return;
end

% --- Error message if CalC simulation failed -----------------------------

if numel(Y1) + numel(Y2) + numel(Y3) + numel(Y4) + numel(Y5) + numel(Y6) + numel(Y7) > 0
    fprintf(" %s\n %s\n %s\n %s\n %s\n %s\n %s\n\n", Y1, Y2, Y3, Y4, Y5, Y6, Y7);
    fprintf("Error computing [Ca2+] using CalC: check README instructions\n");
    return;
else
    fprintf("*** CalC data generated successfully ***\n");
end

% --- Generate exocytosis data curve that the model will attempt fitting:

DT                = logspace(log10(0.5), log10(totalTime), nPoints); 
[DataIn, ErrorIn] = GenerateExocytosisData( DT, 0);

% --- Read computed [Ca2+] vs distance and time, for each buffering condition:

[  ~,    CaGrid1, tArray1] = ReadDataVsTime(fileEGTA2mM );
[  ~,    CaGrid2, tArray2] = ReadDataVsTime(fileEGTA10mM);
[rArray, CaGrid3, tArray3] = ReadDataVsTime(fileBAPTA1mM);

% ================================ END ====================================
