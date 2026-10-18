clear;
clc;
close all;

%% =========================================================
%  BASELINE PARAMETERS
% ==========================================================

c0 = 0.10;       % scaled conversion coefficient
d  = 0.10;       % bacterial loss/death rate (per model time unit)
K  = 14.0;       % scaled carrying capacity
B0 = 0.10;       % scaled initial bacterial abundance

% Baseline uncertain parameters
c_base  = 0.90;  % inulin consumption/interaction coefficient
Ks_base = 2.00;  % Monod half-saturation parameter (relative fibre units)

%% =========================================================
%  INITIAL INULIN AVAILABILITY SCENARIOS
% ==========================================================

Fref = 3.0;
multipliers = [0.5 1.0 2.0];
F0_values = Fref .* multipliers;    % [1.5, 3.0, 6.0]

scenario_names = { ...
    'Low (0.5x baseline)', ...
    'Baseline (1x)', ...
    'High (2x baseline)'};

%% =========================================================
%  SIMULATION SETTINGS
% ==========================================================

tspan = [0 8];    % model time units

opts = odeset( ...
    'NonNegative',[1 2], ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

%% =========================================================
%  MAIN SIMULATIONS USING BASELINE c AND Ks
% ==========================================================

nScen = numel(F0_values);

t_main = cell(nScen,1);
y_main = cell(nScen,1);

PeakB_main      = zeros(nScen,1);
TimeToPeak_main = zeros(nScen,1);
FinalB_main     = zeros(nScen,1);
FinalF_main     = zeros(nScen,1);

for s = 1:nScen
    y0 = [B0; F0_values(s)];

    [t_main{s}, y_main{s}] = ode45( ...
        @(t,y) modelODE(t,y,c0,c_base,d,K,Ks_base), ...
        tspan, y0, opts);

    B = y_main{s}(:,1);
    F = y_main{s}(:,2);

    [PeakB_main(s), idx] = max(B);
    TimeToPeak_main(s)   = t_main{s}(idx);
    FinalB_main(s)       = B(end);
    FinalF_main(s)       = F(end);
end

%% =========================================================
%  FIGURE: BACTERIAL DYNAMICS
% ==========================================================

figure;
hold on;

for s = 1:nScen
    plot(t_main{s}, y_main{s}(:,1), 'LineWidth', 2);
end

xlabel('Time (model units)');
ylabel('Scaled B. adolescentis abundance');
title('Effect of Initial Inulin Availability on B. adolescentis');
legend(scenario_names, 'Location','best');
grid on;
box on;

%% =========================================================
%  FIGURE: INULIN DEPLETION
% ==========================================================

figure;
hold on;

for s = 1:nScen
    plot(t_main{s}, y_main{s}(:,2), 'LineWidth', 2);
end

xlabel('Time (model units)');
ylabel('Relative inulin availability');
title('Inulin Depletion Over Time');
legend(scenario_names, 'Location','best');
grid on;
box on;

%% =========================================================
%  DISPLAY MAIN RESULTS
% ==========================================================

Scenario = scenario_names(:);
InitialInulin = F0_values(:);

main_results = table( ...
    Scenario, InitialInulin, PeakB_main, TimeToPeak_main, ...
    FinalB_main, FinalF_main);

disp(' ');
disp('MAIN RESULTS');
disp(main_results);

%% =========================================================
%  SENSITIVITY ANALYSIS
% ==========================================================
c_values  = [0.70 0.90 1.00];
Ks_values = [1.00 2.00 4.00];

nC  = numel(c_values);
nKs = numel(Ks_values);

PeakB_sens = zeros(nC,nKs,nScen);
ordering_holds = false(nC,nKs);

for ic = 1:nC
    for ik = 1:nKs

        for s = 1:nScen
            y0 = [B0; F0_values(s)];

            [~, y] = ode45( ...
                @(t,y) modelODE(t,y,c0,c_values(ic),d,K,Ks_values(ik)), ...
                tspan, y0, opts);

            PeakB_sens(ic,ik,s) = max(y(:,1));
        end

        lowPeak  = PeakB_sens(ic,ik,1);
        basePeak = PeakB_sens(ic,ik,2);
        highPeak = PeakB_sens(ic,ik,3);

        ordering_holds(ic,ik) = ...
            (highPeak > basePeak) && (basePeak > lowPeak);
    end
end

%% =========================================================
%  PRINT SENSITIVITY RESULTS
% ==========================================================

fprintf('\nSENSITIVITY ANALYSIS\n');
fprintf('-------------------------------------------------------------\n');
fprintf(' c      Ks      Peak_low    Peak_base   Peak_high   Ordering\n');
fprintf('-------------------------------------------------------------\n');

for ic = 1:nC
    for ik = 1:nKs

        if ordering_holds(ic,ik)
            orderText = 'YES';
        else
            orderText = 'NO';
        end

        fprintf('%4.2f   %4.2f     %8.4f    %8.4f    %8.4f      %s\n', ...
            c_values(ic), ...
            Ks_values(ik), ...
            PeakB_sens(ic,ik,1), ...
            PeakB_sens(ic,ik,2), ...
            PeakB_sens(ic,ik,3), ...
            orderText);
    end
end

robust_percent = 100 * sum(ordering_holds(:)) / numel(ordering_holds);

fprintf('\nThe expected ordering High > Baseline > Low held in %.1f%% of tested parameter combinations.\n', ...
    robust_percent);

%% =========================================================
%  FIGURE: SENSITIVITY OF HIGH-FIBRE PEAK TO c AND Ks
% ==========================================================

highPeakMatrix = PeakB_sens(:,:,3);

figure;
imagesc(Ks_values, c_values, highPeakMatrix);
set(gca,'YDir','normal');
colorbar;

xlabel('K_s (relative inulin units)');
ylabel('c (scaled interaction coefficient)');
title('Sensitivity of Peak B. adolescentis Abundance');

%% =========================================================
%  FIGURE: PEAK ABUNDANCE ACROSS MAIN SCENARIOS
% ==========================================================

figure;
bar(1:nScen, PeakB_main);

xticks(1:nScen);
xticklabels({'Low','Baseline','High'});

ylabel('Peak scaled B. adolescentis abundance');
xlabel('Initial inulin scenario');
title('Peak Population Response to Initial Inulin Availability');
grid on;
box on;

%% =========================================================
%  MODEL FUNCTION
% ==========================================================

function dydt = modelODE(~, y, c0, c, d, K, Ks)

    B = y(1);
    F = y(2);

    % Monod-type resource-dependent growth contribution
    g = (c0 * c * F) / (1 + F / Ks);

    % Bacterial dynamics
    dBdt = B * (g - d) * (K - B);

    % Closed-system inulin depletion
    dFdt = -c * F * B;

    dydt = [dBdt;
            dFdt];
end
