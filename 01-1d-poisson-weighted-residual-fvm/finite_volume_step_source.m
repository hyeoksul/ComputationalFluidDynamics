function metrics = finite_volume_step_source(output_dir)
%FINITE_VOLUME_STEP_SOURCE Solve a 1D discontinuous-source Poisson problem.
%   Uses a uniform node-centered finite-volume discretization with
%   T(0)=0 and T(1)=50, then compares against the analytical solution.

if nargin < 1
    output_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

L = 1.0;
number_of_nodes = 41;
x = linspace(0, L, number_of_nodes).';
dx = x(2) - x(1);
% Use control-volume averages at the discontinuity. The node at x=L/2
% represents a volume containing equal portions of S=1 and S=0.
source = double(x < L/2);
interface_node = abs(x - L/2) < 10*eps(L);
source(interface_node) = 0.5;

number_of_interior_nodes = number_of_nodes - 2;
A = diag(-2*ones(number_of_interior_nodes, 1)) ...
  + diag(ones(number_of_interior_nodes - 1, 1), 1) ...
  + diag(ones(number_of_interior_nodes - 1, 1), -1);
b = -source(2:end-1)*dx^2;

left_boundary_value = 0.0;
right_boundary_value = 50.0;
b(1) = b(1) - left_boundary_value;
b(end) = b(end) - right_boundary_value;

temperature = [left_boundary_value; A\b; right_boundary_value];

analytical = zeros(size(x));
left = x <= L/2;
analytical(left) = -0.5*x(left).^2 + 50.375*x(left);
analytical(~left) = 49.875*x(~left) + 0.125;

error = temperature - analytical;
metrics = table("Finite volume", sqrt(mean(error.^2)), max(abs(error)), ...
    'VariableNames', {'Method', 'RMSError', 'MaxAbsoluteError'});

fig = figure('Color', 'w', 'Visible', 'off');
plot(x, analytical, 'k-', 'LineWidth', 2.2); hold on;
plot(x, temperature, 'o', 'MarkerSize', 4.5, 'LineWidth', 1.2);
grid on; box on;
xlabel('x'); ylabel('T(x)');
title('Finite-volume solution of a discontinuous-source Poisson equation');
legend('Analytical', 'Finite volume', 'Location', 'southeast');
exportgraphics(fig, fullfile(output_dir, 'finite_volume_step_source.png'), ...
    'Resolution', 200);
close(fig);

fig = figure('Color', 'w', 'Visible', 'off');
plot(x, abs(error), 'LineWidth', 1.7);
grid on; box on;
xlabel('x'); ylabel('|T_{FVM}-T_{analytical}|');
title('Finite-volume absolute error');
exportgraphics(fig, fullfile(output_dir, 'finite_volume_absolute_error.png'), ...
    'Resolution', 200);
close(fig);
end
