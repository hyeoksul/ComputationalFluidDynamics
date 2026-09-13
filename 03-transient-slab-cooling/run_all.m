%RUN_ALL Reproduce and verify the transient slab-cooling solution.

project_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(project_dir, 'results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end
addpath(project_dir);

intervals = [50; 100; 200];
time_steps = [0.04; 0.02; 0.01];
number_of_cases = numel(intervals);
solutions = cell(number_of_cases, 1);

target_time = zeros(number_of_cases, 1);
midplane_temperature_C = zeros(number_of_cases, 1);
rms_analytical_error_K = zeros(number_of_cases, 1);
maximum_analytical_error_K = zeros(number_of_cases, 1);
symmetry_error_K = zeros(number_of_cases, 1);
maximum_energy_residual = zeros(number_of_cases, 1);
tdma_direct_difference_K = zeros(number_of_cases, 1);

for case_index = 1:number_of_cases
    solutions{case_index} = solve_slab_cooling(...
        intervals(case_index), time_steps(case_index));
    current = solutions{case_index};
    target_time(case_index) = current.time;
    midplane_temperature_C(case_index) = current.midplane_temperature - 273.15;
    rms_analytical_error_K(case_index) = current.rms_analytical_error;
    maximum_analytical_error_K(case_index) = current.maximum_analytical_error;
    symmetry_error_K(case_index) = current.symmetry_error;
    maximum_energy_residual(case_index) = current.maximum_energy_residual;
    tdma_direct_difference_K(case_index) = current.tdma_direct_difference;
end

verification = table(intervals, time_steps, target_time, ...
    midplane_temperature_C, rms_analytical_error_K, ...
    maximum_analytical_error_K, symmetry_error_K, ...
    maximum_energy_residual, tdma_direct_difference_K, ...
    'VariableNames', {'Intervals', 'TimeStepS', 'TargetTimeS', ...
    'MidplaneTemperatureC', 'RMSAnalyticalErrorK', ...
    'MaximumAnalyticalErrorK', 'SymmetryErrorK', ...
    'MaximumEnergyResidual', 'TDMADirectDifferenceK'});
writetable(verification, fullfile(results_dir, 'refinement_verification.csv'));

finest = solutions{end};

fig = figure('Color', 'w', 'Visible', 'off');
plot(finest.time_history, finest.surface_history - 273.15, ...
    'LineWidth', 1.8); hold on;
plot(finest.time_history, finest.midplane_history - 273.15, ...
    '--', 'LineWidth', 1.8);
yline(300, ':', 'Target surface temperature', 'LineWidth', 1.4);
grid on; box on;
xlabel('Time (s)'); ylabel('Temperature (degC)');
title('Surface and midplane temperatures during slab cooling');
legend('Surface', 'Midplane', 'Location', 'northeast');
exportgraphics(fig, fullfile(results_dir, 'cooling_history.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
plot(finest.x, finest.analytical_temperature - 273.15, ...
    'k-', 'LineWidth', 2.2); hold on;
plot(finest.x, finest.temperature - 273.15, 'o', ...
    'MarkerSize', 3.5, 'LineWidth', 1.0);
grid on; box on;
xlabel('Position through slab (m)'); ylabel('Temperature (degC)');
title(sprintf('Temperature profile when the surface reaches 300 degC (t = %.3f s)', ...
    finest.time));
legend('Analytical series', 'Fully implicit FVM', 'Location', 'northwest');
exportgraphics(fig, fullfile(results_dir, 'temperature_profile_at_target.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
yyaxis left;
plot(intervals, rms_analytical_error_K, 'o-', 'LineWidth', 1.8);
ylabel('RMS analytical error (K)');
yyaxis right;
plot(intervals, target_time, 's--', 'LineWidth', 1.8);
ylabel('Predicted target time (s)');
grid on; box on;
xlabel('Number of spatial intervals');
title('Refinement of analytical error and target-time prediction');
legend('RMS error', 'Target time', 'Location', 'best');
exportgraphics(fig, fullfile(results_dir, 'refinement_verification.png'), ...
    'Resolution', 200);
close(fig);

assert(all(symmetry_error_K < 1e-9), 'Symmetry verification failed.');
assert(all(maximum_energy_residual < 1e-9), ...
    'Discrete energy-conservation verification failed.');
assert(all(tdma_direct_difference_K < 1e-9), ...
    'TDMA does not agree with MATLAB direct solution.');
assert(abs(finest.temperature(1) - finest.target_temperature) < 1e-9, ...
    'Event interpolation did not reach the target surface temperature.');
assert(all(diff(rms_analytical_error_K) < 0), ...
    'Refinement did not reduce the analytical-solution error.');
assert(rms_analytical_error_K(end) < 1.0, ...
    'Finest-case analytical error exceeds the verification tolerance.');

expected_files = {
    'refinement_verification.csv'
    'cooling_history.png'
    'temperature_profile_at_target.png'
    'refinement_verification.png'
};
for index = 1:numel(expected_files)
    assert(isfile(fullfile(results_dir, expected_files{index})), ...
        'Expected result file was not generated: %s', expected_files{index});
end

disp(verification);
fprintf('\nVerification passed. Results written to:\n%s\n', results_dir);
