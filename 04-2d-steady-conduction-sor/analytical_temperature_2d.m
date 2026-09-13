function temperature = analytical_temperature_2d(x, y, number_of_terms)
%ANALYTICAL_TEMPERATURE_2D Series solution for the rectangular domain.
%   X and Y must be arrays of equal size containing coordinates at which
%   the solution is evaluated.

if nargin < 3
    number_of_terms = 200;
end

width = 0.3;
height = 0.6;
conductivity = 2.0;
convection_coefficient = 20.0;
fixed_temperature = 200.0;
ambient_temperature = 25.0;
temperature_difference = fixed_temperature - ambient_temperature;

temperature = fixed_temperature*ones(size(x));
for mode = 1:2:number_of_terms
    eigenvalue = mode*pi/height;
    constant_coefficient = 4*temperature_difference/(mode*pi);
    denominator = cosh(eigenvalue*width) ...
        *(conductivity*eigenvalue*tanh(eigenvalue*width) ...
        + convection_coefficient);
    x_factor = cosh(eigenvalue*x)/denominator;
    temperature = temperature ...
        - convection_coefficient*constant_coefficient*x_factor ...
        .*sin(eigenvalue*y);
end
end
