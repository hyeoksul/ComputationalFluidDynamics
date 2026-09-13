function result = solve_composite_heating(cells_in_concrete, cells_in_steel, time_step)
%SOLVE_COMPOSITE_HEATING Implicit FVM for radiative composite-wall heating.

arguments
    cells_in_concrete (1,1) double {mustBeInteger, mustBePositive}
    cells_in_steel (1,1) double {mustBeInteger, mustBePositive}
    time_step (1,1) double {mustBePositive}
end

concrete_length = 0.075;
steel_length = 0.025;
final_time = 10000;

concrete_k = 1.4;
concrete_rho = 2300;
concrete_cp = 880;
steel_k = 55;
steel_rho = 7850;
steel_cp = 450;

ambient_temperature = 900 + 273.15;
initial_temperature = 27 + 273.15;
convection_coefficient = 40;
emissivity = 0.9;
stefan_boltzmann = 5.67e-8;

concrete_dx = concrete_length/cells_in_concrete;
steel_dx = steel_length/cells_in_steel;
assert(abs(concrete_dx - steel_dx) < 100*eps, ...
    'This refinement study expects equal cell widths in both materials.');

cell_width = [repmat(concrete_dx, cells_in_concrete, 1); ...
              repmat(steel_dx, cells_in_steel, 1)];
faces = [0; cumsum(cell_width)];
centers = 0.5*(faces(1:end-1) + faces(2:end));
number_of_cells = numel(cell_width);

conductivity = [repmat(concrete_k, cells_in_concrete, 1); ...
                repmat(steel_k, cells_in_steel, 1)];
density = [repmat(concrete_rho, cells_in_concrete, 1); ...
           repmat(steel_rho, cells_in_steel, 1)];
specific_heat = [repmat(concrete_cp, cells_in_concrete, 1); ...
                 repmat(steel_cp, cells_in_steel, 1)];
thermal_capacity = density.*specific_heat.*cell_width;

face_conductance = zeros(number_of_cells-1, 1);
for face = 1:number_of_cells-1
    resistance = cell_width(face)/(2*conductivity(face)) ...
               + cell_width(face+1)/(2*conductivity(face+1));
    face_conductance(face) = 1/resistance;
end

lower = zeros(number_of_cells, 1);
upper = zeros(number_of_cells, 1);
base_diagonal = thermal_capacity/time_step;
for cell_index = 1:number_of_cells
    if cell_index > 1
        lower(cell_index) = -face_conductance(cell_index-1);
        base_diagonal(cell_index) = base_diagonal(cell_index) ...
            + face_conductance(cell_index-1);
    end
    if cell_index < number_of_cells
        upper(cell_index) = -face_conductance(cell_index);
        base_diagonal(cell_index) = base_diagonal(cell_index) ...
            + face_conductance(cell_index);
    end
end

number_of_steps = round(final_time/time_step);
assert(abs(number_of_steps*time_step - final_time) < 1e-10, ...
    'final_time must be divisible by time_step.');

temperature = initial_temperature*ones(number_of_cells, 1);
surface_temperature = initial_temperature;

time_history = (0:number_of_steps).'*time_step;
surface_history = zeros(number_of_steps+1, 1);
mid_concrete_history = zeros(number_of_steps+1, 1);
core_history = zeros(number_of_steps+1, 1);
convective_flux_history = zeros(number_of_steps+1, 1);
radiative_flux_history = zeros(number_of_steps+1, 1);

surface_history(1) = surface_temperature;
mid_concrete_history(1) = interp1(centers, temperature, 0.05, 'linear');
core_history(1) = temperature(end);
convective_flux_history(1) = convection_coefficient ...
    *(ambient_temperature - surface_temperature);
radiative_flux_history(1) = emissivity*stefan_boltzmann ...
    *(ambient_temperature^4 - surface_temperature^4);
maximum_energy_residual = 0;
maximum_nonlinear_iterations = 0;
maximum_surface_iteration_error = 0;

for step = 1:number_of_steps
    previous_temperature = temperature;
    previous_energy = sum(thermal_capacity.*previous_temperature);
    surface_guess = surface_temperature;

    for nonlinear_iteration = 1:50
        radiation_coefficient = emissivity*stefan_boltzmann ...
            *(ambient_temperature + surface_guess) ...
            *(ambient_temperature^2 + surface_guess^2);
        external_coefficient = convection_coefficient + radiation_coefficient;
        boundary_conductance = 1/(cell_width(1)/(2*conductivity(1)) ...
            + 1/external_coefficient);

        diagonal = base_diagonal;
        diagonal(1) = diagonal(1) + boundary_conductance;
        rhs = thermal_capacity/time_step.*previous_temperature;
        rhs(1) = rhs(1) + boundary_conductance*ambient_temperature;
        candidate_temperature = tdma_solver(lower, diagonal, upper, rhs);

        total_boundary_flux = boundary_conductance ...
            *(ambient_temperature - candidate_temperature(1));
        candidate_surface_temperature = ambient_temperature ...
            - total_boundary_flux/external_coefficient;
        surface_iteration_error = ...
            abs(candidate_surface_temperature - surface_guess);
        surface_guess = candidate_surface_temperature;
        if surface_iteration_error < 1e-10
            break;
        end
    end
    assert(nonlinear_iteration < 50, ...
        'Surface radiation iteration failed to converge.');

    temperature = candidate_temperature;
    surface_temperature = candidate_surface_temperature;
    maximum_nonlinear_iterations = max(maximum_nonlinear_iterations, ...
        nonlinear_iteration);
    maximum_surface_iteration_error = max(maximum_surface_iteration_error, ...
        surface_iteration_error);

    convective_flux = convection_coefficient ...
        *(ambient_temperature - surface_temperature);
    radiative_flux = emissivity*stefan_boltzmann ...
        *(ambient_temperature^4 - surface_temperature^4);
    current_energy = sum(thermal_capacity.*temperature);
    energy_residual = abs((current_energy - previous_energy)/time_step ...
        - (convective_flux + radiative_flux)) ...
        /max(abs(convective_flux + radiative_flux), 1);
    maximum_energy_residual = max(maximum_energy_residual, energy_residual);

    surface_history(step+1) = surface_temperature;
    mid_concrete_history(step+1) = interp1(...
        centers, temperature, 0.05, 'linear');
    core_history(step+1) = temperature(end);
    convective_flux_history(step+1) = convective_flux;
    radiative_flux_history(step+1) = radiative_flux;
end

matrix = diag(diagonal) + diag(upper(1:end-1), 1) ...
       + diag(lower(2:end), -1);
direct_temperature = matrix\rhs;

interface_flux = face_conductance(cells_in_concrete) ...
    *(temperature(cells_in_concrete) - temperature(cells_in_concrete+1));
interface_temperature = temperature(cells_in_concrete) ...
    - interface_flux*cell_width(cells_in_concrete) ...
    /(2*conductivity(cells_in_concrete));

result = struct();
result.faces = faces;
result.centers = centers;
result.temperature = temperature;
result.surface_temperature = surface_temperature;
result.interface_temperature = interface_temperature;
result.core_temperature = temperature(end);
result.mid_concrete_temperature = interp1(centers, temperature, 0.05, 'linear');
result.interface_flux = interface_flux;
result.time_history = time_history;
result.surface_history = surface_history;
result.mid_concrete_history = mid_concrete_history;
result.core_history = core_history;
result.convective_flux_history = convective_flux_history;
result.radiative_flux_history = radiative_flux_history;
result.maximum_energy_residual = maximum_energy_residual;
result.maximum_nonlinear_iterations = maximum_nonlinear_iterations;
result.maximum_surface_iteration_error = maximum_surface_iteration_error;
result.tdma_direct_difference = max(abs(temperature - direct_temperature));
result.minimum_temperature_increment = min(diff(...
    [surface_history, mid_concrete_history, core_history]), [], 'all');
end
