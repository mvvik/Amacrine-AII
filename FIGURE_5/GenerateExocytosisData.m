function [DataIn, ErrorIn] = GenerateExocytosisData( DT, plotFlag )

    arrayY0   = [-24.75   -21.49  -9.919];
    arrayPlat = [81.97     77.16     120];
    arrayTau1 = [0.8858   0.9133    1.97];
    arrayTau2 = [39.38     30.25   458.5];
    arrayFast = [66.75     74.11   33.07]/100;

    if plotFlag
        figure; hold on;
        legendStr = {'2mM EGTA', '10mM EGTA', '1mM BAPTA'};
        clr       = [0 0 0; 0 0 1; 1 0 0];
    end

    DataIn = zeros(3, numel(DT));

    for k = 1 : 3
        Y0 = arrayY0(k);
        Plat = arrayPlat(k);
        Tau1 = arrayTau1(k);
        Tau2 = arrayTau2(k);
        Fast = arrayFast(k);

        Y = Plat + (Y0-Plat) * (Fast*exp(-DT/Tau1) + (1-Fast)*exp(-DT/Tau2));
        if plotFlag
            plot(DT, Y, 'o-', 'color', clr(k,:), 'linewidth', 2);
        end
        DataIn(k,:) = Y;
    end

    In1 = (DataIn(1,:) + DataIn(2,:)) / 2;
    DataEGTA  = In1;
    DataBAPTA = DataIn(3,:);
    DataIn    = [DataEGTA, DataEGTA, DataBAPTA];

    ErrorEGTA  = 1.5 + DT * 5.5 / 20;
    ErrorBAPTA = 1 + tanh(0.2*(DT-0.5))*5.1;

    ErrorIn = [ErrorEGTA, ErrorEGTA, ErrorBAPTA].^(-2);

    if plotFlag
        set(gca, 'xscale', 'log');
        legend(legendStr{1}, legendStr{2}, legendStr{3});
    end
end