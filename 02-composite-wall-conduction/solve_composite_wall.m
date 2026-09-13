function result = solve_composite_wall(cells_in_A, cells_in_B)
%SOLVE_COMPOSITE_WALL Solve a two-material wall with heat generation.
%   RESULT = SOLVE_COMPOSITE_WALL(NA, NB) uses NA and NB cell-centered
%   control volumes in materials A and B. The interface coincides with a
%   control-volume face.

arguments
    cells_in_A (1,1) double {mustBeInteger, mustBePositive}
    cells_in_B (1,1) double {mustBeInteger, mustBePositive}
end

% Geometry and properties (SI units).
length_A = 0.05;
length_B = 0.02;
conductivity_A = 75.0;
conductivity_B = 150.0;
generation_A = 1.5e6;
convection_coefficient = 1000.0;
ambient_temperature = 303.15;

cell_width = [repmat(length_A/cells_in_A, cells_in_A, 1); ...
              repmat(length_B/cells_in_B, cells_in_B, 1)];
faces = [0; cumsum(cell_width)];
centers = 0.5*(faces(1:end-1) + faces(2:end));
number_of_cells = numel(cell_width);

conductivity = [repmat(conductivity_A, cells_in_A, 1); ...
                repmat(conductivity_B, cells_in_B, 1)];
generation = [repmat(generation_A, cells_in_A, 1); ...
              zeros(cells_in_B, 1)];

% Face conductance per unit wall area. Writing the coefficient as a sum of
% half-cell resistances is valid for unequal cell widths and discontinuous k.
face_conductance = zeros(number_of_cells-1, 1);
for face = 1:number_of_cells-1
    resistance = cell_width(face)/(2*conductivity(face)) ...
               + cell_width(face+1)/(2*conductivity(face+1));
    face_conductance(face) = 1/resistance;
end

lower = zeros(number_of_cells, 1);
diagonal = zeros(number_of_cells, 1);
upper = zeros(number_of_cells, 1);
rhs = generation.*cell_width;

% Insulated left boundary: only the east-face heat flux remains.
diagonal(1) = face_conductance(1);
upper(1) = -face_conductance(1);

for cell_index = 2:number_of_cells-1
    lower(cell_index) = -face_conductance(cell_index-1);
    diagonal(cell_index) = face_conductance(cell_index-1) ...
                         + face_conductance(cell_index);
    upper(cell_index) = -face_conductance(cell_index);
end

% Combine conduction across the final half cell with external convection.
right_resistance = cell_width(end)/(2*conductivity(end)) ...
                 + 1/convection_coefficient;
right_conductance = 1/right_resistance;
lower(end) = -face_conductance(end);
diagonal(end) = face_conductance(end) + right_conductance;
rhs(end) = rhs(end) + right_conductance*ambient_temperature;

temperature = tdma_solver(lower, diagonal, upper, rhs);

matrix = diag(diagonal) + diag(upper(1:end-1), 1) ...
       + diag(lower(2:end), -1);
direct_temperature = matrix\rhs;

total_length = length_A + length_B;
analytical_heat_flux = generation_A*length_A;
surface_temperature = ambient_temperature ...
                    + analytical_heat_flux/convection_coefficient;
interface_temperature = surface_temperature ...
                      + analytical_heat_flux*length_B/conductivity_B;

analytical_temperature = zeros(number_of_cells, 1);
in_A = centers <= length_A;
analytical_temperature(in_A) = interface_temperature ...
    + generation_A/(2*conductivity_A)*(length_A^2 - centers(in_A).^2);
analytical_temperature(~in_A) = surface_temperature ...
    + analytical_heat_flux/conductivity_B*(total_length - centers(~in_A));

temperature_error = temperature - analytical_temperature;
generated_heat = sum(generation.*cell_width);
rejected_heat = right_conductance*(temperature(end) - ambient_temperature);
interface_heat_flux = face_conductance(cells_in_A) ...
                    *(temperature(cells_in_A) - temperature(cells_in_A+1));
numerical_surface_temperature = ambient_temperature ...
                              + rejected_heat/convection_coefficient;

result = struct();
result.faces = faces;
result.centers = centers;
result.temperature = temperature;
result.analytical_temperature = analytical_temperature;
result.surface_temperature = numerical_surface_temperature;
result.interface_temperature = interface_temperature;
result.generated_heat = generated_heat;
result.rejected_heat = rejected_heat;
result.interface_heat_flux = interface_heat_flux;
result.analytical_heat_flux = analytical_heat_flux;
result.rms_error = sqrt(mean(temperature_error.^2));
result.maximum_error = max(abs(temperature_error));
result.energy_imbalance = abs(generated_heat - rejected_heat)/generated_heat;
result.interface_flux_relative_error = ...
    abs(interface_heat_flux - analytical_heat_flux)/analytical_heat_flux;
result.tdma_direct_difference = max(abs(temperature - direct_temperature));
end
