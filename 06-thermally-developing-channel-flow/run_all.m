%RUN_ALL Reproduce and verify thermally developing channel-flow results.

project_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(project_dir,'results');
if ~isfolder(results_dir), mkdir(results_dir); end
addpath(project_dir);

%% Convection-scheme comparison on the reporting grid
schemes = ["Central";"Upwind";"Hybrid";"PowerLaw"];
reporting_cell_size = 0.0625;
scheme_solutions = cell(numel(schemes),1);
outlet_nusselt = zeros(numel(schemes),1);
correlation_rms_error = zeros(numel(schemes),1);
minimum_temperature = zeros(numel(schemes),1);
maximum_temperature = zeros(numel(schemes),1);

for index = 1:numel(schemes)
    current = solve_channel_temperature(reporting_cell_size, ...
        schemes(index),'Direct');
    scheme_solutions{index} = current;
    outlet_nusselt(index) = current.local_nusselt(end);
    comparison = current.x >= 0.5 & current.x <= 0.9*current.x(end);
    correlation_rms_error(index) = sqrt(mean((current.local_nusselt(comparison) ...
        -current.correlation_nusselt(comparison)).^2));
    minimum_temperature(index) = min(current.temperature,[],'all');
    maximum_temperature(index) = max(current.temperature,[],'all');
end

scheme_summary = table(schemes,outlet_nusselt,correlation_rms_error, ...
    minimum_temperature,maximum_temperature, ...
    'VariableNames',{'Scheme','OutletNusselt','CorrelationRMSError', ...
    'MinimumTemperatureK','MaximumTemperatureK'});
writetable(scheme_summary,fullfile(results_dir,'scheme_comparison.csv'));

%% Grid refinement using the bounded power-law scheme
cell_size = [0.5;0.25;0.125;0.0625];
monitor_temperature = zeros(size(cell_size));
outlet_nusselt_grid = zeros(size(cell_size));
grid_unknowns = zeros(size(cell_size));
for index = 1:numel(cell_size)
    current = solve_channel_temperature(cell_size(index),'PowerLaw','Direct');
    [~,ix] = min(abs(current.x-12.5));
    [~,jy] = min(abs(current.y-1.0));
    monitor_temperature(index) = current.temperature(jy,ix);
    outlet_nusselt_grid(index) = current.local_nusselt(end);
    grid_unknowns(index) = current.nx*current.ny;
end
difference_from_finest = abs(monitor_temperature-monitor_temperature(end));
grid_summary = table(cell_size,grid_unknowns,monitor_temperature, ...
    outlet_nusselt_grid,difference_from_finest, ...
    'VariableNames',{'CellSize','Unknowns','TemperatureAt12p5And1K', ...
    'OutletNusselt','DifferenceFromFinestK'});
writetable(grid_summary,fullfile(results_dir,'grid_refinement.csv'));

%% Linear-solver comparison on a common algebraic system
solver_names = ["Jacobi";"GaussSeidel";"SOR"];
relaxation = [1;1;1.1];
iterations = zeros(size(relaxation));
elapsed_time_s = zeros(size(relaxation));
final_scaled_residual = zeros(size(relaxation));
direct_difference_K = zeros(size(relaxation));
solver_solutions = cell(numel(solver_names),1);
for index = 1:numel(solver_names)
    current = solve_channel_temperature(0.25,'PowerLaw', ...
        solver_names(index),relaxation(index),1e-8,20000);
    assert(current.converged,'The %s solver did not converge.',solver_names(index));
    solver_solutions{index} = current;
    iterations(index) = current.iterations;
    elapsed_time_s(index) = current.elapsed_time;
    final_scaled_residual(index) = current.final_scaled_residual;
    direct_difference_K(index) = current.direct_difference;
end
solver_summary = table(solver_names,relaxation,iterations,elapsed_time_s, ...
    final_scaled_residual,direct_difference_K, ...
    'VariableNames',{'Solver','Relaxation','Iterations','ElapsedTimeS', ...
    'FinalScaledResidual','DirectDifferenceK'});
writetable(solver_summary,fullfile(results_dir,'solver_comparison.csv'));

%% Verification
primary = scheme_solutions{schemes=="PowerLaw"};
symmetry_error = max(abs(primary.temperature-flipud(primary.temperature)),[],'all');
assert(primary.final_scaled_residual < 1e-10,'Direct-solver residual is too large.');
assert(symmetry_error < 1e-9,'The temperature field is not symmetric.');
assert(all(minimum_temperature >= primary.inlet_temperature-1e-9));
assert(all(maximum_temperature <= primary.wall_temperature+1e-9));
assert(all(diff(difference_from_finest(1:end-1)) < 0), ...
    'Grid refinement is not monotonic relative to the finest result.');
assert(all(direct_difference_K < 2e-5), ...
    'An iterative solution differs excessively from the direct solution.');

verification = table(primary.reynolds_gap,primary.reynolds_hydraulic, ...
    primary.prandtl,primary.fully_developed_x,symmetry_error, ...
    primary.final_scaled_residual, ...
    'VariableNames',{'ReynoldsBasedOnGap','ReynoldsBasedOnHydraulicDiameter', ...
    'Prandtl','NumericalFullyDevelopedXOverL','SymmetryErrorK', ...
    'DirectScaledResidual'});
writetable(verification,fullfile(results_dir,'verification_summary.csv'));

%% Figures
fig = figure('Color','w','Visible','off');
contourf(primary.x/primary.half_gap,primary.y/primary.half_gap, ...
    primary.temperature,24,'LineColor','none');
colorbar; colormap(turbo); axis tight;
xlabel('x/l'); ylabel('y/l');
title('Thermally developing temperature field (K)');
exportgraphics(fig,fullfile(results_dir,'temperature_field.png'),'Resolution',200);
close(fig);

fig = figure('Color','w','Visible','off'); hold on;
line_styles = {'-','--','-.',':'};
for index = 1:numel(schemes)
    plot(scheme_solutions{index}.x(2:end), ...
        scheme_solutions{index}.local_nusselt(2:end),line_styles{index}, ...
        'LineWidth',1.5);
end
plot(primary.x(2:end),primary.correlation_nusselt(2:end),'k-', ...
    'LineWidth',2.2);
yline(7.541,'k:','Fully developed: 7.541');
grid on; box on; xlim([0,15]); ylim([6,35]);
xlabel('x/l'); ylabel('Local Nu_{D_h}');
title('Local Nusselt-number development');
legend([cellstr(schemes);{'Shah-London correlation'}],'Location','northeast');
exportgraphics(fig,fullfile(results_dir,'nusselt_development.png'),'Resolution',200);
close(fig);

fig = figure('Color','w','Visible','off');
loglog(cell_size(1:end-1),difference_from_finest(1:end-1),'o-', ...
    'LineWidth',1.8,'MarkerSize',7); grid on; box on;
set(gca,'XDir','reverse');
xlabel('Uniform cell size'); ylabel('|T-T_{finest}| at (12.5l,l) (K)');
title('Grid-refinement verification');
exportgraphics(fig,fullfile(results_dir,'grid_refinement.png'),'Resolution',200);
close(fig);

fig = figure('Color','w','Visible','off');
for index = 1:numel(solver_names)
    semilogy(1:solver_solutions{index}.iterations, ...
        solver_solutions{index}.residual_history,'LineWidth',1.5);
    hold on;
end
set(gca,'YScale','log');
grid on; box on;
xlabel('Iteration'); ylabel('Scaled maximum algebraic residual');
title('Stationary iterative-solver convergence');
legend(cellstr(solver_names),'Location','northeast');
exportgraphics(fig,fullfile(results_dir,'solver_convergence.png'),'Resolution',200);
close(fig);

disp(scheme_summary); disp(grid_summary); disp(solver_summary); disp(verification);
fprintf('\nAll verification checks passed. Results written to:\n%s\n',results_dir);
