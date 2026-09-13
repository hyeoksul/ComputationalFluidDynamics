%RUN_ALL Reproduce and verify the composite-wall finite-volume solution.

project_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(project_dir, 'results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end
addpath(project_dir);

cells_A = [5; 10; 20; 40];
cells_B = [2; 4; 8; 16];
number_of_meshes = numel(cells_A);

rms_error = zeros(number_of_meshes, 1);
maximum_error = zeros(number_of_meshes, 1);
energy_imbalance = zeros(number_of_meshes, 1);
interface_flux_relative_error = zeros(number_of_meshes, 1);
tdma_direct_difference = zeros(number_of_meshes, 1);
observed_order = nan(number_of_meshes, 1);

solutions = cell(number_of_meshes, 1);
for mesh = 1:number_of_meshes
    solutions{mesh} = solve_composite_wall(cells_A(mesh), cells_B(mesh));
    rms_error(mesh) = solutions{mesh}.rms_error;
    maximum_error(mesh) = solutions{mesh}.maximum_error;
    energy_imbalance(mesh) = solutions{mesh}.energy_imbalance;
    interface_flux_relative_error(mesh) = ...
        solutions{mesh}.interface_flux_relative_error;
    tdma_direct_difference(mesh) = solutions{mesh}.tdma_direct_difference;
end

mesh_spacing_A = 0.05./cells_A;
for mesh = 2:number_of_meshes
    observed_order(mesh) = log(rms_error(mesh-1)/rms_error(mesh)) ...
        /log(mesh_spacing_A(mesh-1)/mesh_spacing_A(mesh));
end

total_cells = cells_A + cells_B;
mesh_verification = table(cells_A, cells_B, total_cells, rms_error, ...
    maximum_error, energy_imbalance, interface_flux_relative_error, ...
    tdma_direct_difference, observed_order, 'VariableNames', ...
    {'CellsInA', 'CellsInB', 'TotalCells', 'RMSErrorK', ...
     'MaximumErrorK', 'EnergyImbalance', 'InterfaceFluxRelativeError', ...
     'TDMADirectDifferenceK', 'ObservedOrder'});
writetable(mesh_verification, fullfile(results_dir, 'mesh_verification.csv'));

finest = solutions{end};
conservation_summary = table(finest.generated_heat, finest.rejected_heat, ...
    finest.interface_heat_flux, finest.analytical_heat_flux, ...
    finest.interface_temperature, finest.surface_temperature, ...
    'VariableNames', {'GeneratedHeatWPerM2', 'RejectedHeatWPerM2', ...
    'InterfaceHeatFluxWPerM2', 'AnalyticalHeatFluxWPerM2', ...
    'InterfaceTemperatureK', 'SurfaceTemperatureK'});
writetable(conservation_summary, ...
    fullfile(results_dir, 'conservation_summary.csv'));

fig = figure('Color', 'w', 'Visible', 'off');
x_dense = linspace(0, 0.07, 600).';
analytical_dense = zeros(size(x_dense));
in_A = x_dense <= 0.05;
analytical_dense(in_A) = 413.15 - 10000*x_dense(in_A).^2;
analytical_dense(~in_A) = 413.15 - 500*x_dense(~in_A);
plot(x_dense*1000, analytical_dense, 'k-', 'LineWidth', 2.2); hold on;
plot(finest.centers*1000, finest.temperature, 'o', ...
    'MarkerSize', 4.5, 'LineWidth', 1.1);
xline(50, '--', 'Material interface', 'LabelVerticalAlignment', 'bottom');
grid on; box on;
xlabel('x (mm)'); ylabel('Temperature (K)');
title('Steady temperature distribution in the composite wall');
legend('Analytical', 'Finite volume', 'Location', 'southwest');
exportgraphics(fig, fullfile(results_dir, 'temperature_distribution.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
second_order_reference = rms_error(1)*(total_cells/total_cells(1)).^(-2);
loglog(total_cells, rms_error, 'o-', 'LineWidth', 1.8, ...
    'MarkerSize', 7); hold on;
loglog(total_cells, second_order_reference, '--', 'LineWidth', 1.6);
grid on; box on;
xlabel('Total number of control volumes'); ylabel('RMS temperature error (K)');
title('Mesh convergence of the composite-wall solution');
legend('Finite-volume error', 'Second-order reference', ...
    'Location', 'southwest');
exportgraphics(fig, fullfile(results_dir, 'mesh_convergence.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
values = [finest.generated_heat, finest.rejected_heat, ...
          finest.interface_heat_flux, finest.analytical_heat_flux]/1000;
bar(values);
xticklabels({'Generated', 'Rejected', 'At interface', 'Analytical'});
ylabel('Heat flux (kW/m^2)');
title('Global conservation and interface-flux verification');
grid on; box on;
ylim([0, 1.15*max(values)]);
exportgraphics(fig, fullfile(results_dir, 'conservation_verification.png'), ...
    'Resolution', 200);
close(fig);

assert(all(observed_order(2:end) > 1.99), ...
    'The temperature solution did not demonstrate second-order convergence.');
assert(maximum_error(end) < 0.005, ...
    'The finest-mesh temperature error exceeds the verification tolerance.');
assert(all(energy_imbalance < 1e-10), ...
    'Global energy-conservation verification failed.');
assert(all(interface_flux_relative_error < 1e-10), ...
    'Interface heat-flux verification failed.');
assert(all(tdma_direct_difference < 1e-10), ...
    'TDMA does not agree with MATLAB direct solution.');

expected_files = {
    'mesh_verification.csv'
    'conservation_summary.csv'
    'temperature_distribution.png'
    'mesh_convergence.png'
    'conservation_verification.png'
};
for index = 1:numel(expected_files)
    assert(isfile(fullfile(results_dir, expected_files{index})), ...
        'Expected result file was not generated: %s', expected_files{index});
end

disp(mesh_verification);
disp(conservation_summary);
fprintf('\nVerification passed. Results written to:\n%s\n', results_dir);
