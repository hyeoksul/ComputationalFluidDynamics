function temperature = analytical_slab_temperature(x, time, number_of_terms)
%ANALYTICAL_SLAB_TEMPERATURE Eigenfunction solution for convective cooling.

if nargin < 3
    number_of_terms = 100;
end

thickness = 0.2;
half_thickness = thickness/2;
initial_temperature = 1400 + 273.15;
ambient_temperature = 50 + 273.15;
convection_coefficient = 5000;
conductivity = 30;
density = 7800;
specific_heat = 700;
thermal_diffusivity = conductivity/(density*specific_heat);

biot_number = convection_coefficient*half_thickness/conductivity;
fourier_number = thermal_diffusivity*time/half_thickness^2;
centered_coordinate = abs(x - half_thickness)/half_thickness;

dimensionless_temperature = zeros(size(x));
root_offset = 1e-10;
for mode = 1:number_of_terms
    lower_bound = (mode-1)*pi + root_offset;
    upper_bound = (mode-0.5)*pi - root_offset;
    eigenvalue = fzero(@(z) z*tan(z) - biot_number, ...
        [lower_bound, upper_bound]);
    coefficient = 4*sin(eigenvalue) ...
        /(2*eigenvalue + sin(2*eigenvalue));
    dimensionless_temperature = dimensionless_temperature ...
        + coefficient*cos(eigenvalue*centered_coordinate) ...
        .*exp(-eigenvalue^2*fourier_number);
end

temperature = ambient_temperature ...
    + (initial_temperature - ambient_temperature)*dimensionless_temperature;
end
