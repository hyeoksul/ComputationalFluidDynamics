%RUN_ALL Reproduce and verify the 2D steady-conduction project.

project_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(project_dir, 'results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end
addpath(project_dir);

%% Relaxation-factor benchmark on a fixed grid
relaxation_factor = [0.8; 1.0; 1.5; 1.8; 1.9];
number_of_factors = numel(relaxation_factor);
iterations = zeros(number_of_factors, 1);
final_residual = zeros(number_of_factors, 1);
center_temperature_C = zeros(number_of_factors, 1);
heat_output_W_per_m = zeros(number_of_factors, 1);

benchmark_solutions = cell(number_of_factors, 1);
for index = 1:number_of_factors
    benchmark_solutions{index} = solve_steady_conduction_2d(...
        24, 48, relaxation_factor(index), 1e-11);
    current = benchmark_solutions{index};
    iterations(index) = current.iterations;
    final_residual(index) = current.final_residual;
    center_temperature_C(index) = current.center_temperature;
    heat_output_W_per_m(index) = current.heat_output;
end

relaxation_benchmark = table(relaxation_factor, iterations, ...
    final_residual, center_temperature_C, heat_output_W_per_m, ...
    'VariableNames', {'RelaxationFactor', 'Iterations', 'FinalResidual', ...
    'CenterTemperatureC', 'RightBoundaryHeatRateWPerM'});
writetable(relaxation_benchmark, ...
    fullfile(results_dir, 'relaxation_benchmark.csv'));

%% Grid refinement using a near-optimal relaxation factor
nx_values = [8; 16; 32; 64];
ny_values = 2*nx_values;
number_of_grids = numel(nx_values);
rms_analytical_error_C = zeros(number_of_grids, 1);
maximum_analytical_error_C = zeros(number_of_grids, 1);
energy_imbalance = zeros(number_of_grids, 1);
direct_solution_difference_C = zeros(number_of_grids, 1);
grid_iterations = zeros(number_of_grids, 1);
grid_solutions = cell(number_of_grids, 1);

for index = 1:number_of_grids
    grid_solutions{index} = solve_steady_conduction_2d(...
        nx_values(index), ny_values(index), 1.9, 1e-11);
    current = grid_solutions{index};
    rms_analytical_error_C(index) = current.rms_analytical_error;
    maximum_analytical_error_C(index) = current.maximum_analytical_error;
    energy_imbalance(index) = current.energy_imbalance;
    direct_solution_difference_C(index) = current.direct_solution_difference;
    grid_iterations(index) = current.iterations;
end

grid_refinement = table(nx_values, ny_values, grid_iterations, ...
    rms_analytical_error_C, maximum_analytical_error_C, ...
    energy_imbalance, direct_solution_difference_C, ...
    'VariableNames', {'Nx', 'Ny', 'Iterations', 'RMSAnalyticalErrorC', ...
    'MaximumAnalyticalErrorC', 'EnergyImbalance', ...
    'DirectSolutionDifferenceC'});
writetable(grid_refinement, fullfile(results_dir, 'grid_refinement.csv'));

finest = grid_solutions{end};

fig = figure('Color', 'w', 'Visible', 'off');
contourf(finest.X, finest.Y, finest.temperature, 24, 'LineColor', 'none');
hold on;
contour(finest.X, finest.Y, finest.temperature, 12, ...
    'LineColor', [0.2, 0.2, 0.2], 'LineWidth', 0.6);
colorbar; axis equal tight;
xlabel('x (m)'); ylabel('y (m)');
title('Steady temperature field (degC)');
exportgraphics(fig, fullfile(results_dir, 'temperature_field.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
bar(relaxation_factor, iterations, 0.65);
grid on; box on;
xlabel('SOR relaxation factor, \omega'); ylabel('Iterations to residual < 10^{-9}');
title('Effect of relaxation factor on SOR convergence');
exportgraphics(fig, fullfile(results_dir, 'relaxation_benchmark.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
loglog(nx_values.*ny_values, rms_analytical_error_C, ...
    'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on; box on;
xlabel('Number of control volumes'); ylabel('RMS analytical error (degC)');
title('Grid refinement against the analytical series solution');
exportgraphics(fig, fullfile(results_dir, 'grid_refinement.png'), ...
    'Resolution', 200);
close(fig);

assert(max(center_temperature_C) - min(center_temperature_C) < 1e-6, ...
    'Relaxation factor changed the converged center temperature.');
assert(max(heat_output_W_per_m) - min(heat_output_W_per_m) < 1e-5, ...
    'Relaxation factor changed the converged heat rate.');
assert(all(final_residual < 1e-11), 'A relaxation benchmark did not converge.');
assert(all(diff(rms_analytical_error_C) < 0), ...
    'Grid refinement did not reduce analytical error.');
assert(all(energy_imbalance < 1e-8), 'Global heat balance failed.');
assert(all(direct_solution_difference_C < 1e-6), ...
    'SOR does not agree with MATLAB sparse direct solution.');

expected_files = {
    'relaxation_benchmark.csv'
    'grid_refinement.csv'
    'temperature_field.png'
    'relaxation_benchmark.png'
    'grid_refinement.png'
};
for index = 1:numel(expected_files)
    assert(isfile(fullfile(results_dir, expected_files{index})), ...
        'Expected result file was not generated: %s', expected_files{index});
end

disp(relaxation_benchmark);
disp(grid_refinement);
fprintf('\nVerification passed. Results written to:\n%s\n', results_dir);
