% Base case values
basePerm = 100 * milli * darcy;
basePoro = 0.3;
baseInjectionRate = 200 * meter^3/day;

% Initialize MRST and add necessary modules
mrstModule add ad-core ad-blackoil ad-props mrst-gui

% Define the grid dimensions and physical parameters
nx = 60; ny = 60; nz = 1;
Lx = 1000; Ly = 1000; Lz = 10;
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
G = computeGeometry(G);

% Rock and fluid properties
rock.perm = basePerm * ones(G.cells.num, 1);
rock.poro = basePoro * ones(G.cells.num, 1);
fluid = initSimpleADIFluid('phases', 'WOG', ...
    'mu', [1, 10, 0.05] * centi * poise, ...
    'rho', [1000, 800, 600], ...
    'n', [2, 2, 2]);

% Initial state
state = initResSol(G, 200 * barsa, [0.2, 0.8, 0.0]);

% Well configuration
W = [];
W = addWell(W, G, rock, 1, 'Type', 'rate', 'Val', baseInjectionRate, ...
    'Radius', 0.1, 'Comp_i', [0, 1, 0]); % CO2 injection well
W = addWell(W, G, rock, G.cells.num, 'Type', 'bhp', 'Val', 150 * barsa, ...
    'Radius', 0.1, 'Comp_i', [0, 0, 1]); % Oil production well

% Simulation schedule
nSteps = 50;
dt = repmat(10 * day, nSteps, 1);
schedule = simpleSchedule(dt, 'W', W);

% Model setup
model = ThreePhaseBlackOilModel(G, rock, fluid);

% Run the base simulation
[wellSols, states] = simulateScheduleAD(state, model, schedule);

% Calculate base oil recovery efficiency
baseTotalOilProduced = sum(diff(cat(1, wellSols.oil)));
baseInitialOilVolume = sum(state.s(:, 2) .* rock.poro .* G.cells.volumes);
baseRecoveryFactor = baseTotalOilProduced / baseInitialOilVolume;

% Define a range of permeability values to test
permValues = [50, 100, 200] * milli * darcy;
recoveryFactorsPerm = zeros(length(permValues), 1);

for i = 1:length(permValues)
    % Update permeability
    rock.perm = permValues(i) * ones(G.cells.num, 1);

    % Run the simulation
    [wellSols, states] = simulateScheduleAD(state, model, schedule);

    % Calculate oil recovery efficiency
    totalOilProduced = sum(diff(cat(1, wellSols.oil)));
    recoveryFactorsPerm(i) = totalOilProduced / baseInitialOilVolume;
end

% Plot results
figure;
plot(permValues / milli, recoveryFactorsPerm, '-o');
xlabel('Permeability (mD)');
ylabel('Oil Recovery Factor');
title('Sensitivity of Oil Recovery to Permeability');

poroValues = [0.2, 0.3, 0.4];
recoveryFactorsPoro = zeros(length(poroValues), 1);

for i = 1:length(poroValues)
    rock.poro = poroValues(i) * ones(G.cells.num, 1);
    [wellSols, states] = simulateScheduleAD(state, model, schedule);
    totalOilProduced = sum(diff(cat(1, wellSols.oil)));
    recoveryFactorsPoro(i) = totalOilProduced / baseInitialOilVolume;
end

figure;
plot(poroValues, recoveryFactorsPoro, '-o');
xlabel('Porosity');
ylabel('Oil Recovery Factor');
title('Sensitivity of Oil Recovery to Porosity');

injectionRates = [100, 200, 300] * meter^3/day;
recoveryFactorsInjectionRate = zeros(length(injectionRates), 1);

for i = 1:length(injectionRates)
    W(1).val = injectionRates(i);
    [wellSols, states] = simulateScheduleAD(state, model, schedule);
    totalOilProduced = sum(diff(cat(1, wellSols.oil)));
    recoveryFactorsInjectionRate(i) = totalOilProduced / baseInitialOilVolume;
end

figure;
plot(injectionRates, recoveryFactorsInjectionRate, '-o');
xlabel('Injection Rate (m^3/day)');
ylabel('Oil Recovery Factor');
title('Sensitivity of Oil Recovery to Injection Rate');