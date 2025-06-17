clear;
nCaSites    = 5;
GetModel    = @GetModelOutputShorter5;
GetEquil    = @getEquilibriumNonCoop5;
GetAffinity = @findHalfMaxAffinity5;

actualAffinityFlag =  0;
flagCorr           =  1;
nBins              = 40; 
clr                = [0.6, 0.0, 0.6];

CollectParForData;

    %     1     2       3       4     5     6      7   
    %    R1     R2     SCaKD  SCakp  coop  gamma  frac

%            1       2       3        4                5        6         7        8            
labels = { 'R_X',  'R_Y' , 'K_D',   'k_{on}',       '\beta', '\gamma', 'Pool_X', 'Pool_Y'};
Units  = { 'nm',   'nm',   '\muM', '(\muM ms)^{-1}',  '',    'ms^{-1}',   '%',     '%' }; 

Order  = [ 1 3 7 5  ...
           2 4 8 6 ];
nRows = 2;
nCols = 4;

plfs = 12;
tfs = 12;
lfs = 12;

ind1 = 1; % Special parameter index 1
ind2 = 2; % Special parameter index 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
indFrac    = nPars;  % index of pool fraction parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ResultsOut(:, 1:nPars) = SetParamBounds( ResultsOut(:, 1:nPars), Rpatch, 0 );

[Vmin, Imin] = min( ResultsOut(:,end) );
ParMin       = ResultsOut(Imin, 1:nPars);
fprintf(['min = %g \n Pars = %g ', fStr, '\n'], Vmin, ParMin);

MaxError           = min(ResultsOut(:,end)) + LLR_CutOff;
indTotal           = find( ResultsOut(:,end) <= MaxError);
ResultsOut         = ResultsOut(indTotal, :);
ResultsOut(:, end) = (ResultsOut(:,end) - Vmin);

nData      = size(ResultsOut, 1);
fprintf('Filtered set: N = %d\n', nData );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if actualAffinityFlag
    for kk = 1 : numel(indTotal)
        K1 = ResultsOut(kk, 3);
        cp = ResultsOut(kk, 5);
        ResultsOut(kk, 3) = GetAffinity(K1, cp);
    end
end

ResultsOut(:,1)     = ResultsOut(:,1)     * 1000;
ResultsOut(:,2)     = ResultsOut(:,2)     * 1000;
ResultsOut(:,nPars) = ResultsOut(:,nPars) * 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nPars = nPars + 1;
ResultsNew              = zeros( nData, nPars + 1);
ResultsNew(:, 1:nPars ) = ResultsOut;
ResultsNew(:, nPars)    = 100 - ResultsOut(:, nPars-1);
ResultsNew(:, end  )    = ResultsOut(:, end);
ResultsOut              = ResultsNew;
ParMin                  = ResultsOut(1, 1:nPars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN  = size(ResultsOut, 1);

% ResultsOut(:, 9) = ResultsOut(:, 9) * 1000;  % Convert pA to nA
for jj = 1 : NN
    
    R1  = ResultsOut(jj, 1);
    R2  = ResultsOut(jj, 2);
    f1  = ResultsOut(jj, nPars-1);
    f2  = ResultsOut(jj, nPars  );
    
    if R2 < R1
        ResultsOut(jj, 1) = R2;
        ResultsOut(jj, 2) = R1;
        ResultsOut(jj, nPars-1) = f2;
        ResultsOut(jj, nPars  ) = f1;
    end
end

for jj = 1 : NN
    % ResultsOut(jj, 1:nPars) = SetParamBounds( ResultsOut(jj, 1:nPars), 0); %%%%%%%%
    ResultsOut(jj, 1:nPars) = ResultsOut(jj, Order); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flagCorr
    NN  = size(ResultsOut, 1);
    M1  = zeros(nPars, 1);
    M2  = zeros(nPars, 1);
    STD = zeros(nPars, 1);
    
    for mm = 1 : nPars
            M1(mm) = sum(ResultsOut(:, mm)  ) / NN;
            M2(mm) = sum(ResultsOut(:, mm).^2) / NN;
            STD(mm) = sqrt(M2(mm) - M1(mm)^2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    CRR     = ones(nPars, nPars);
    
    for mm = 2 : nPars
         for nn = 1 : mm - 1
             C = sum( (ResultsOut(:, mm) - M1(mm)) .* (ResultsOut(:, nn) - M1(nn)) );
             C = C / STD(mm) / STD(nn) / NN;
             CRR(mm,nn) = C;
             CRR(nn,mm) = C;
             CorrStr = string([labels{Order(nn)}, ' vs ', labels{Order(mm)}]);
             if abs(C) > 0.1
                fprintf('%22s %g\n', CorrStr, C);
             end
        end
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f = figure;
    hold off;
    
    PCORR = [CRR,   zeros(nPars,1) ];
    PCORR = [PCORR; zeros(1,nPars+1) ];
    pcolor(1:(nPars+1), 1:(nPars+1), real(PCORR) ); hold on;
    colormap gray;
    colormap(f, flipud(colormap(f)))
    colorbar('NorthOutside');
    hold on;
    
    shading flat;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    
    for ii = 1 : nPars
        text(ii+0.22, 0.5, labels{Order(ii)}, 'fontsize', plfs); 
        text(0.3,  ii+0.4, labels{Order(ii)}, 'fontsize', plfs);
        plot([ii ii], [0 nPars+1],  '-', 'color', [1 1 1], 'linewidth', 1);
        plot([0 nPars+1], [ii ii],  '-', 'color', [1 1 1], 'linewidth', 1);
    end
    
    title('{\bf (B)}   Parameter Correlations', 'fontsize', tfs);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     [coeff,score,latent] = pca(ResultsOut(:,1:nPars));
%     figure; 
%     plot(log(latent), 'k-o');
%     xlabel('Principal Component ID');
%     ylabel('Log(Principal Component Variance)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; 

Vmin     = 48.539360294084780;
MaxError = Vmin + LLR_CutOff;

for kk = 1 : nPars

    Xmin = min( ResultsOut(:, kk) );
    Xmax = max( ResultsOut(:, kk) );
    if kk == 8
            Xmin = log( Xmin );
            Xmax = log( Xmax );
    end
    valPars   = zeros(1, nBins);
    valErrs   = zeros(1, nBins);

    for indBin = 1 : nBins

        valPar0 = Xmin + (indBin - 1) / nBins * (Xmax - Xmin);
        valPar1 = Xmin +  indBin      / nBins * (Xmax - Xmin);
        valPar  = 0.5 * (valPar0 + valPar1);

        if kk == 8
            valPar0 = exp(valPar0);
            valPar1 = exp(valPar1);
            valPar  = exp( valPar);
        end
        valPars(indBin) = valPar;

        %fprintf("Par = %d, val = %g .. %g \n", indPar, valPar0, valPar1);
        inds1 = find(ResultsOut(:,     kk) >= valPar0);
        if numel(inds1) > 0
            inds2 = find(ResultsOut(inds1, kk) <  valPar1);
        else
            inds2 = [];
        end

        if numel(inds2) > 0 
            RES             = ResultsOut(inds1(inds2), :);
            [xx, ind]       = min(RES(:,end));
            valErrs(indBin) = RES(ind, end);
        else
            [C, ind]        = min( (ResultsOut(:, kk) - valPar).^2 );
            valErrs(indBin) = ResultsOut(ind, end);
        end

        %fprintf(['finished %g --> %g \n  Pars: ', fStr, '\n\n'], cost0, cost, Y);
    end
   
    if kk == 8
        Xmin = exp( Xmin );
        Xmax = exp( Xmax );
    end
    
    subplot(nRows,nCols,kk); hold on;
    valPars = [ valPars(1),       valPars,   valPars(end) ];
    valErrs = [ LLR_CutOff,       valErrs,    LLR_CutOff];

    fill(valPars, valErrs,  clr, 'EdgeColor', clr, 'LineWidth', 1);

    %plot(valPars, valErrs, 'm-');
    axis([Xmin Xmax 0 LLR_CutOff]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for mm = 0 : 6

   ind1 = find(ResultsOut(:,    1) >= 20 + mm*2);
   ind2 = find(ResultsOut(ind1, 1) <  22 + mm*2);
   clrB  = (0.9 - mm/10);
   clrG  = 1 - clrB;
   clr   = [ clrG, clrB, 0 ];
   sz    = (9 - mm)/3;

    for nn = [1 5] %1 : nPars 
    
        XX = ResultsOut(ind1(ind2),end);
        YY = ResultsOut(ind1(ind2),nn);
    
        subplot(nRows,nCols,nn); hold on;
        plot(YY, XX, 'o' , 'markerfacecolor', clr, 'color', clr, 'MarkerSize', sz);
        %plot(ParTrue(kk)*[1 1], [min(XX), max(XX)], 'r-', 'linewidth', 2);
    end
end


for mm = 1 : 6

   ind1 = find(ResultsOut(:,    6) >= (mm + 3)/10);
   ind2 = find(ResultsOut(ind1, 6) <  (mm + 4)/10);
   clrB  = (0.9 - mm/10);
   clrG  = 1 - clrB;
   clr   = [ clrG, 0, clrB ];
   sz    = (9 - mm)/3;

    for nn = [6 8] %1 : nPars 
    
        XX = ResultsOut(ind1(ind2),end);
        YY = ResultsOut(ind1(ind2),nn);
    
        subplot(nRows,nCols,nn); hold on;
        plot(YY, XX, 'o' , 'markerfacecolor', clr, 'color', clr, 'MarkerSize', sz);
        %plot(ParTrue(kk)*[1 1], [min(XX), max(XX)], 'r-', 'linewidth', 2);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nn = 1 : nPars
        
        kk = Order(nn);
        subplot(nRows,nCols,nn); 
        plot( ParMin(kk)*[1 1], [0, LLR_CutOff], 'w:', 'linewidth', 2);
        
        if ParMin(kk) < 10
            STR1 = sprintf('%s=%.2g%s', labels{kk}, ParMin(kk), Units{kk} );
            %STR2 = sprintf('(true: %s=%.3g%s)', labels{kk}, ParTrue(kk), Units{kk} );
        else
            STR1 = sprintf('%s=%.0f%s', labels{kk}, ParMin(kk), Units{kk} );
            %STR2 = sprintf('(true: %s=%.0f%s)', labels{kk}, ParTrue(kk), Units{kk} );
        end
        title(STR1, 'fontsize',    tfs, 'color', 'k');
        %subtitle(STR2, 'fontsize', tfs, 'color', 'r');
        %tt = title({STR1, STR2}, 'fontsize', tfs);
        
        %ylabel('Log Likelihood Ratio', 'fontsize', lfs);
        %xlabel(labels{kk}, 'fontsize', lfs);
        set(gca, 'Ydir',   'reverse');
        %set(gca, 'Xscale', 'log'    );
        if nn == nPars
            set(gca, 'xscale', 'log');
            set(gca,'xtick',[3 5 10 20 50])
        end
        axis tight;
end

