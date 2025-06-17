DefineLocalParameters;

DT                = logspace(log10(0.5), log10(totalTime), nPoints); 
[DataIn, ErrorIn] = GenerateExocytosisData( DT, 0);

Kappa      = 50;
BKD        = 1;
Bkplus     = 0.4;
ICa        = 0.2;               
CaTau0     = 10;
CaTau      = CaTau0 / (1 + Kappa);
nGrid      = 100;
LLR_CutOff = 8.0;

fStr = '';
for ind = 1 : nPars
   fStr = [fStr, '%g '];
end

prefix       = ['Rpatch_', num2str(1000*Rpatch), 'nm_ICa_', num2str(ICa), 'pA_CaTau', num2str(CaTau0), 'ms_N', num2str(nGrid), '_'];
prefixN      = [ prefix, 'CaSites', num2str(nCaSites), '_'];

nameEGTA2mM  = ['Data/', prefix, 'EGTA2mM'];
nameEGTA10mM = ['Data/', prefix, 'EGTA10mM'];
nameBAPTA1mM = ['Data/', prefix, 'BAPTA1mM'];

Script = '  GenerateCalciumData.par';

if contains(computer, 'WIN')
    prog   = ['..\cwin694x64.exe ', Script];
else
    prog   = ['../cmac697x11 '    , Script];
end

        
         