function result = solve_slab_cooling(number_of_intervals, time_step)
%SOLVE_SLAB_COOLING Fully implicit FVM solution for a convectively cooled slab.

arguments
    number_of_intervals (1,1) double {mustBeInteger, mustBePositive}
    time_step (1,1) double {mustBePositive}
end
assert(mod(number_of_intervals, 2) == 0, ...
    'Use an even number of intervals so a node lies at the midplane.');

thickness = 0.2;
initial_temperature = 1400 + 273.15;
ambient_temperature = 50 + 273.15;
target_temperature = 300 + 273.15;
convection_coefficient = 5000;
conductivity = 30;
density = 7800;
specific_heat = 700;

number_of_nodes = number_of_intervals + 1;
dx = thickness/number_of_intervals;
x = linspace(0, thickness, number_of_nodes).';
control_volume_width = dx*ones(number_of_nodes, 1);
control_volume_width([1, end]) = dx/2;
thermal_capacity = density*specific_heat*control_volume_width;
conductance = conductivity/dx;

lower = -conductance*ones(number_of_nodes, 1);
diagonal = thermal_capacity/time_step + 2*conductance;
upper = -conductance*ones(number_of_nodes, 1);
lower(1) = 0;
upper(end) = 0;
diagonal(1) = thermal_capacity(1)/time_step ...
            + conductance + convection_coefficient;
diagonal(end) = thermal_capacity(end)/time_step ...
              + conductance + convection_coefficient;

matrix = diag(diagonal) + diag(upper(1:end-1), 1) ...
       + diag(lower(2:end), -1);

temperature = initial_temperature*ones(number_of_nodes, 1);
time = 0;
maximum_time = 1000;
midplane_index = number_of_intervals/2 + 1;

time_history = 0;
surface_history = initial_temperature;
midplane_history = initial_temperature;
maximum_energy_residual = 0;
tdma_direct_difference = 0;

while time < maximum_time
    previous_temperature = temperature;
    previous_energy = sum(thermal_capacity ...
        .*(previous_temperature - ambient_temperature));

    rhs = thermal_capacity/time_step.*previous_temperature;
    rhs([1, end]) = rhs([1, end]) ...
        + convection_coefficient*ambient_temperature;

    temperature = tdma_solver(lower, diagonal, upper, rhs);
    time = time + time_step;

    current_energy = sum(thermal_capacity ...
        .*(temperature - ambient_temperature));
    convection_loss = convection_coefficient ...
        *((temperature(1) - ambient_temperature) ...
        + (temperature(end) - ambient_temperature));
    energy_residual = abs((current_energy - previous_energy)/time_step ...
        + convection_loss)/max(convection_loss, eps);
    maximum_energy_residual = max(maximum_energy_residual, energy_residual);

    time_history(end+1, 1) = time; %#ok<AGROW>
    surface_history(end+1, 1) = 0.5*(temperature(1) + temperature(end)); %#ok<AGROW>
    midplane_history(end+1, 1) = temperature(midplane_index); %#ok<AGROW>

    if temperature(1) <= target_temperature
        previous_surface = previous_temperature(1);
        event_fraction = (previous_surface - target_temperature) ...
            /(previous_surface - temperature(1));
        event_time = time - time_step + event_fraction*time_step;
        event_temperature = previous_temperature ...
            + event_fraction*(temperature - previous_temperature);
        tdma_direct_difference = max(abs(temperature - matrix\rhs));
        break;
    end
end

assert(exist('event_time', 'var') == 1, ...
    'The target surface temperature was not reached before maximum_time.');

analytical_temperature = analytical_slab_temperature(x, event_time, 100);
temperature_error = event_temperature - analytical_temperature;

result = struct();
result.x = x;
result.time = event_time;
result.temperature = event_temperature;
result.analytical_temperature = analytical_temperature;
result.time_history = time_history;
result.surface_history = surface_history;
result.midplane_history = midplane_history;
result.target_temperature = target_temperature;
result.midplane_temperature = event_temperature(midplane_index);
result.rms_analytical_error = sqrt(mean(temperature_error.^2));
result.maximum_analytical_error = max(abs(temperature_error));
result.symmetry_error = max(abs(event_temperature - flipud(event_temperature)));
result.maximum_energy_residual = maximum_energy_residual;
result.tdma_direct_difference = tdma_direct_difference;
end
