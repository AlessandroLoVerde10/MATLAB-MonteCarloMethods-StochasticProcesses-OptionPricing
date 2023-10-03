 %% - ole 

clear all
clc
close all
[num] = xlsread('market_data')

MarketStrikes = num(1,9:end)/10000
MarketVolatilities =  num(2:8,9:end)'
CurrentForwardValues = num(2:8,3)
NumMaturities = 7
Beta = 0.5

Settle = datenum('15-Sep-2013');

ExerciseDate = [Settle + (364/4); Settle + (364/2); Settle + (364/4)*3; Settle + 365; Settle + 365*2; Settle + 365*5; Settle + 365*10]

% Calibrate Alpha, Rho, and Nu
MarketStrike = [MarketStrikes; MarketStrikes;MarketStrikes;MarketStrikes;MarketStrikes;MarketStrikes;MarketStrikes]'
Betas = repmat(Beta, NumMaturities, 1)
Alphas = zeros(NumMaturities, 1);
Rhos = zeros(NumMaturities, 1);
Nus = zeros(NumMaturities, 1);

%options = optimoptions('lsqnonlin','Display','none');

for k = 1:NumMaturities
    % Fit Rho and Nu (while converting at-the-money volatility into Alpha)
    objFun = @(X) MarketVolatilities(:,k) -  blackvolbysabr(X(1), Betas(k),  X(2),  X(3), Settle, ExerciseDate(k), CurrentForwardValues(k), MarketStrike(:,k));
    
    X = lsqnonlin(objFun, [0.5 0 0.5], [0 -1 0], [Inf 1 Inf]);
    
    Alphas(k) = X(1);
    Rhos(k) = X(2);
    Nus(k) = X(3);
end

CalibratedPrameters = [Alphas Betas Rhos Nus]

%% plot
PlottingStrikes = (25:5:150)'/10000;

ComputedVols = zeros(length(PlottingStrikes), NumMaturities);

for k = 1:NumMaturities
ComputedVols(:,k) = blackvolbysabr(CalibratedPrameters(k,1), ...
CalibratedPrameters(k,2), CalibratedPrameters(k,3), ....
CalibratedPrameters(k,4), Settle, ExerciseDate(k), CurrentForwardValues(k), PlottingStrikes);
end

%Errors = sum((((ComputedVols(1,:))' - MarketVolatilities(:,1)).^2 + ((ComputedVols(6,:))' - MarketVolatilities(:,2)).^2+ ((ComputedVols(16,:))' - MarketVolatilities(:,3)).^2 + (ComputedVols(26,:)' - MarketVolatilities(:,4)).^2)./4)
YearsToExercise = yearfrac(Settle, ExerciseDate, 1);
figure
surf(YearsToExercise, PlottingStrikes, ComputedVols)
xlim([0 10]); ylim([0.0025 0.0150]); zlim([0.1 0.5]);
view(113,32);
set(gca, 'Position', [0.13 0.11 0.775 0.815], 'PlotBoxAspectRatioMode', 'manual');
xlabel('Years to exercise', 'Fontweight', 'bold');
ylabel('Strike', 'Fontweight', 'bold');
zlabel('Implied Black volatility', 'Fontweight', 'bold');
hold on 
PlottingStrikes2 = [25, 50, 100, 150]'/10000;
MarketVolatiliti = MarketVolatilities'
MarketVolatilitiescatter = [MarketVolatiliti(:,1); MarketVolatiliti(:,2); MarketVolatiliti(:,3); MarketVolatiliti(:,4)];
YearsToExerciseScatter = [YearsToExercise;YearsToExercise;YearsToExercise;YearsToExercise];
PlottingStrikesScatter = [25;25;25;25;25;25;25;50;50;50;50;50;50;50;100;100;100;100;100;100;100;150;150;150;150;150;150;150]/10000;
scatter3( YearsToExerciseScatter, PlottingStrikesScatter ,MarketVolatilitiescatter, 'filled')
