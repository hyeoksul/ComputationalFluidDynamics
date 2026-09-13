function result = solve_steady_conduction_2d(nx, ny, omega, tolerance)
%SOLVE_STEADY_CONDUCTION_2D Cell-centered FVM and point-SOR solution.

arguments
    nx (1,1) double {mustBeInteger, mustBePositive}
    ny (1,1) double {mustBeInteger, mustBePositive}
    omega (1,1) double {mustBeGreaterThan(omega,0), mustBeLessThan(omega,2)}
    tolerance (1,1) double {mustBePositive} = 1e-9
end

width = 0.3;
height = 0.6;
conductivity = 2.0;
convection_coefficient = 20.0;
ambient_temperature = 25.0;
fixed_temperature = 200.0;
maximum_iterations = 200000;

dx = width/nx;
dy = height/ny;
x = ((1:nx) - 0.5)*dx;
y = ((1:ny) - 0.5)*dy;
[X, Y] = meshgrid(x, y);

east_west_conductance = conductivity*dy/dx;
north_south_conductance = conductivity*dx/dy;
dirichlet_conductance = conductivity*dx/(dy/2);
convection_conductance = dy/(dx/(2*conductivity) ...
    + 1/convection_coefficient);

aP = zeros(ny, nx);
aW = zeros(ny, nx);
aE = zeros(ny, nx);
aS = zeros(ny, nx);
aN = zeros(ny, nx);
b = zeros(ny, nx);

for row = 1:ny
    for column = 1:nx
        if column > 1
            aW(row, column) = east_west_conductance;
            aP(row, column) = aP(row, column) + east_west_conductance;
        end
        if column < nx
            aE(row, column) = east_west_conductance;
            aP(row, column) = aP(row, column) + east_west_conductance;
        else
            aP(row, column) = aP(row, column) + convection_conductance;
            b(row, column) = b(row, column) ...
                + convection_conductance*ambient_temperature;
        end
        if row > 1
            aS(row, column) = north_south_conductance;
            aP(row, column) = aP(row, column) + north_south_conductance;
        else
            aP(row, column) = aP(row, column) + dirichlet_conductance;
            b(row, column) = b(row, column) ...
                + dirichlet_conductance*fixed_temperature;
        end
        if row < ny
            aN(row, column) = north_south_conductance;
            aP(row, column) = aP(row, column) + north_south_conductance;
        else
            aP(row, column) = aP(row, column) + dirichlet_conductance;
            b(row, column) = b(row, column) ...
                + dirichlet_conductance*fixed_temperature;
        end
    end
end

temperature = 0.5*(fixed_temperature + ambient_temperature)*ones(ny, nx);
residual_history = zeros(maximum_iterations, 1);

for iteration = 1:maximum_iterations
    for row = 1:ny
        for column = 1:nx
            neighbor_sum = b(row, column);
            if column > 1
                neighbor_sum = neighbor_sum ...
                    + aW(row, column)*temperature(row, column-1);
            end
            if column < nx
                neighbor_sum = neighbor_sum ...
                    + aE(row, column)*temperature(row, column+1);
            end
            if row > 1
                neighbor_sum = neighbor_sum ...
                    + aS(row, column)*temperature(row-1, column);
            end
            if row < ny
                neighbor_sum = neighbor_sum ...
                    + aN(row, column)*temperature(row+1, column);
            end
            predicted_temperature = neighbor_sum/aP(row, column);
            temperature(row, column) = temperature(row, column) ...
                + omega*(predicted_temperature - temperature(row, column));
        end
    end

    equation_residual = compute_residual(temperature, aP, aW, aE, aS, aN, b);
    normalized_residual = max(abs(equation_residual), [], 'all') ...
        /max(abs(b), [], 'all');
    residual_history(iteration) = normalized_residual;
    if normalized_residual < tolerance
        break;
    end
end
assert(iteration < maximum_iterations, 'SOR failed to converge.');
residual_history = residual_history(1:iteration);

number_of_unknowns = nx*ny;
matrix = spalloc(number_of_unknowns, number_of_unknowns, 5*number_of_unknowns);
rhs = zeros(number_of_unknowns, 1);
for row = 1:ny
    for column = 1:nx
        index = sub2ind([ny, nx], row, column);
        matrix(index, index) = aP(row, column);
        rhs(index) = b(row, column);
        if column > 1
            matrix(index, sub2ind([ny, nx], row, column-1)) = -aW(row, column);
        end
        if column < nx
            matrix(index, sub2ind([ny, nx], row, column+1)) = -aE(row, column);
        end
        if row > 1
            matrix(index, sub2ind([ny, nx], row-1, column)) = -aS(row, column);
        end
        if row < ny
            matrix(index, sub2ind([ny, nx], row+1, column)) = -aN(row, column);
        end
    end
end
direct_temperature = reshape(matrix\rhs, ny, nx);

analytical_temperature = analytical_temperature_2d(X, Y, 300);
analytical_error = temperature - analytical_temperature;

heat_input_bottom = sum(dirichlet_conductance ...
    *(fixed_temperature - temperature(1, :)));
heat_input_top = sum(dirichlet_conductance ...
    *(fixed_temperature - temperature(end, :)));
heat_output_right = sum(convection_conductance ...
    *(temperature(:, end) - ambient_temperature));
energy_imbalance = abs(heat_input_bottom + heat_input_top - heat_output_right) ...
    /heat_output_right;

result = struct();
result.x = x;
result.y = y;
result.X = X;
result.Y = Y;
result.temperature = temperature;
result.analytical_temperature = analytical_temperature;
result.iterations = iteration;
result.residual_history = residual_history;
result.final_residual = residual_history(end);
result.rms_analytical_error = sqrt(mean(analytical_error.^2, 'all'));
result.maximum_analytical_error = max(abs(analytical_error), [], 'all');
result.energy_imbalance = energy_imbalance;
result.heat_input = heat_input_bottom + heat_input_top;
result.heat_output = heat_output_right;
result.direct_solution_difference = max(abs(temperature - direct_temperature), [], 'all');
result.center_temperature = temperature(round(ny/2), round(nx/2));
end

function residual = compute_residual(temperature, aP, aW, aE, aS, aN, b)
[ny, nx] = size(temperature);
residual = b - aP.*temperature;
for row = 1:ny
    for column = 1:nx
        if column > 1
            residual(row, column) = residual(row, column) ...
                + aW(row, column)*temperature(row, column-1);
        end
        if column < nx
            residual(row, column) = residual(row, column) ...
                + aE(row, column)*temperature(row, column+1);
        end
        if row > 1
            residual(row, column) = residual(row, column) ...
                + aS(row, column)*temperature(row-1, column);
        end
        if row < ny
            residual(row, column) = residual(row, column) ...
                + aN(row, column)*temperature(row+1, column);
        end
    end
end
end
