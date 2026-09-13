%RUN_ALL Reproduce and verify transient concrete-steel heating.

project_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(project_dir, 'results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end
addpath(project_dir);

concrete_cells = [30; 60; 120; 240];
steel_cells = [10; 20; 40; 80];
time_step_s = [4; 2; 1; 0.5];
number_of_cases = numel(concrete_cells);
solutions = cell(number_of_cases, 1);

surface_temperature_C = zeros(number_of_cases, 1);
temperature_at_50mm_C = zeros(number_of_cases, 1);
interface_temperature_C = zeros(number_of_cases, 1);
core_temperature_C = zeros(number_of_cases, 1);
maximum_energy_residual = zeros(number_of_cases, 1);
maximum_nonlinear_iterations = zeros(number_of_cases, 1);
tdma_direct_difference_K = zeros(number_of_cases, 1);
minimum_temperature_increment_K = zeros(number_of_cases, 1);

for case_index = 1:number_of_cases
    solutions{case_index} = solve_composite_heating(...
        concrete_cells(case_index), steel_cells(case_index), ...
        time_step_s(case_index));
    current = solutions{case_index};
    surface_temperature_C(case_index) = current.surface_temperature - 273.15;
    temperature_at_50mm_C(case_index) = current.mid_concrete_temperature - 273.15;
    interface_temperature_C(case_index) = current.interface_temperature - 273.15;
    core_temperature_C(case_index) = current.core_temperature - 273.15;
    maximum_energy_residual(case_index) = current.maximum_energy_residual;
    maximum_nonlinear_iterations(case_index) = current.maximum_nonlinear_iterations;
    tdma_direct_difference_K(case_index) = current.tdma_direct_difference;
    minimum_temperature_increment_K(case_index) = current.minimum_temperature_increment;
end

reference_vector = [surface_temperature_C(end), temperature_at_50mm_C(end), ...
    interface_temperature_C(end), core_temperature_C(end)];
reference_difference_C = zeros(number_of_cases, 1);
for case_index = 1:number_of_cases
    current_vector = [surface_temperature_C(case_index), ...
        temperature_at_50mm_C(case_index), interface_temperature_C(case_index), ...
        core_temperature_C(case_index)];
    reference_difference_C(case_index) = norm(current_vector - reference_vector, 2);
end

total_cells = concrete_cells + steel_cells;
refinement = table(concrete_cells, steel_cells, total_cells, time_step_s, ...
    surface_temperature_C, temperature_at_50mm_C, interface_temperature_C, ...
    core_temperature_C, reference_difference_C, maximum_energy_residual, ...
    maximum_nonlinear_iterations, tdma_direct_difference_K, ...
    'VariableNames', {'ConcreteCells', 'SteelCells', 'TotalCells', ...
    'TimeStepS', 'SurfaceTemperatureC', 'TemperatureAt50mmC', ...
    'InterfaceTemperatureC', 'CoreTemperatureC', ...
    'DifferenceFromFinestC', 'MaximumEnergyResidual', ...
    'MaximumNonlinearIterations', 'TDMADirectDifferenceK'});
writetable(refinement, fullfile(results_dir, 'refinement_verification.csv'));

finest = solutions{end};
final_summary = table(finest.surface_temperature - 273.15, ...
    finest.mid_concrete_temperature - 273.15, ...
    finest.interface_temperature - 273.15, ...
    finest.core_temperature - 273.15, finest.interface_flux, ...
    'VariableNames', {'SurfaceTemperatureC', 'TemperatureAt50mmC', ...
    'InterfaceTemperatureC', 'CoreTemperatureC', 'InterfaceFluxWPerM2'});
writetable(final_summary, fullfile(results_dir, 'final_summary.csv'));

fig = figure('Color', 'w', 'Visible', 'off');
plot([0; finest.centers; 0.1]*1000, ...
    [finest.surface_temperature; finest.temperature; finest.core_temperature] ...
    - 273.15, 'LineWidth', 2.0);
xline(75, '--', 'Concrete-steel interface', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
grid on; box on;
xlabel('Depth (mm)'); ylabel('Temperature (degC)');
title('Temperature distribution after 10,000 s');
exportgraphics(fig, fullfile(results_dir, 'final_temperature_profile.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
plot(finest.time_history, finest.surface_history - 273.15, ...
    'LineWidth', 1.8); hold on;
plot(finest.time_history, finest.mid_concrete_history - 273.15, ...
    '--', 'LineWidth', 1.8);
plot(finest.time_history, finest.core_history - 273.15, ...
    '-.', 'LineWidth', 1.8);
grid on; box on;
xlabel('Time (s)'); ylabel('Temperature (degC)');
title('Temperature histories in the composite');
legend('Exposed surface', 'Concrete at 50 mm', 'Adiabatic steel boundary', ...
    'Location', 'southeast');
exportgraphics(fig, fullfile(results_dir, 'temperature_histories.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
plot(finest.time_history, finest.convective_flux_history/1000, ...
    'LineWidth', 1.8); hold on;
plot(finest.time_history, finest.radiative_flux_history/1000, ...
    '--', 'LineWidth', 1.8);
grid on; box on;
xlabel('Time (s)'); ylabel('Inward heat flux (kW/m^2)');
title('Convective and radiative surface heating');
legend('Convection', 'Radiation', 'Location', 'northeast');
exportgraphics(fig, fullfile(results_dir, 'boundary_heat_flux.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
semilogy(total_cells(1:end-1), reference_difference_C(1:end-1), ...
    'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on; box on;
xlabel('Total number of control volumes');
ylabel('Four-temperature difference from finest case (degC)');
title('Coupled grid and time-step refinement');
exportgraphics(fig, fullfile(results_dir, 'refinement_verification.png'), ...
    'Resolution', 200);
close(fig);

assert(all(maximum_energy_residual < 1e-8), ...
    'Global transient energy balance failed.');
assert(all(tdma_direct_difference_K < 1e-9), ...
    'TDMA does not agree with MATLAB direct solution.');
assert(all(minimum_temperature_increment_K > -1e-9), ...
    'A monitored temperature decreased during heating.');
assert(all(diff(reference_difference_C(1:end-1)) < 0), ...
    'Coupled refinement did not approach the finest result monotonically.');
assert(all(maximum_nonlinear_iterations < 50), ...
    'Nonlinear boundary iteration did not converge.');

expected_files = {
    'refinement_verification.csv'
    'final_summary.csv'
    'final_temperature_profile.png'
    'temperature_histories.png'
    'boundary_heat_flux.png'
    'refinement_verification.png'
};
for index = 1:numel(expected_files)
    assert(isfile(fullfile(results_dir, expected_files{index})), ...
        'Expected result file was not generated: %s', expected_files{index});
end

disp(refinement);
disp(final_summary);
fprintf('\nVerification passed. Results written to:\n%s\n', results_dir);
