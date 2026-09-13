function result = solve_channel_temperature(cell_size, scheme, solver, ...
    relaxation, tolerance, maximum_iterations)
%SOLVE_CHANNEL_TEMPERATURE Solve the channel energy equation.

if nargin < 4, relaxation = 1; end
if nargin < 5, tolerance = 1e-8; end
if nargin < 6, maximum_iterations = 20000; end

model = build_channel_system(cell_size, scheme);
A = model.A;
b = model.b;
temperature_scale = model.wall_temperature-model.inlet_temperature;
diagonal_values = diag(A);
scaled_residual = @(vector) max(abs(A*vector-b) ...
    ./max(abs(diagonal_values), eps))/temperature_scale;

tic;
if strcmpi(solver, 'Direct')
    vector = A\b;
    iterations = 1;
    residual_history = scaled_residual(vector);
    converged = true;
    direct_difference = 0;
else
    vector = model.inlet_temperature*ones(size(b));
    vector(model.boundary_indices) = b(model.boundary_indices);
    diagonal = spdiags(diagonal_values, 0, size(A,1), size(A,2));
    lower_matrix = tril(A,-1);
    upper_matrix = triu(A,1);
    residual_history = nan(maximum_iterations,1);

    for iterations = 1:maximum_iterations
        switch lower(string(solver))
            case "jacobi"
                vector = diagonal\(b-(lower_matrix+upper_matrix)*vector);
            case "gaussseidel"
                vector = (diagonal+lower_matrix)\(b-upper_matrix*vector);
            case "sor"
                vector = (diagonal+relaxation*lower_matrix)\(...
                    relaxation*b-((relaxation-1)*diagonal ...
                    +relaxation*upper_matrix)*vector);
            otherwise
                error('Unknown linear solver: %s', solver);
        end
        residual_history(iterations) = scaled_residual(vector);
        if residual_history(iterations) < tolerance
            break;
        end
    end
    residual_history = residual_history(1:iterations);
    converged = residual_history(end) < tolerance;
    direct_solution = A\b;
    direct_difference = max(abs(vector-direct_solution));
end
elapsed_time = toc;

temperature = reshape(vector, model.ny, model.nx);
heat_transfer = compute_nusselt_number(model, temperature);

result = model;
result.temperature = temperature;
result.solver = string(solver);
result.relaxation = relaxation;
result.iterations = iterations;
result.residual_history = residual_history;
result.final_scaled_residual = scaled_residual(vector);
result.converged = converged;
result.elapsed_time = elapsed_time;
result.direct_difference = direct_difference;
result.bulk_temperature = heat_transfer.bulk_temperature;
result.local_nusselt = heat_transfer.local_nusselt;
result.correlation_nusselt = heat_transfer.correlation_nusselt;
result.fully_developed_x = heat_transfer.fully_developed_x;
end
