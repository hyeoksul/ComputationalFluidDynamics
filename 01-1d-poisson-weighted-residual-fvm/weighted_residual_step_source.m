function metrics = weighted_residual_step_source(output_dir)
%WEIGHTED_RESIDUAL_STEP_SOURCE Compare three weighted-residual approximations.
%   Solves phi'' + S(x) = 0 on [0,L], where S is unity on the first
%   half of the domain and zero on the second half. Homogeneous Dirichlet
%   conditions are prescribed at both ends.

if nargin < 1
    output_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

L = 1.0;
x = linspace(0, L, 1001);

exact = zeros(size(x));
left = x <= L/2;
exact(left) = -0.5*x(left).^2 + (3*L/8)*x(left);
exact(~left) = -(L/8)*x(~left) + L^2/8;

galerkin = (2*L^2/pi^3)*sin(pi*x/L) ...
         + (L^2/(2*pi^3))*sin(2*pi*x/L) ...
         + (2*L^2/(27*pi^3))*sin(3*pi*x/L);

point_collocation = ((2 + sqrt(2))*L^2/(4*pi^2))*sin(pi*x/L) ...
                  + (L^2/(8*pi^2))*sin(2*pi*x/L) ...
                  + ((sqrt(2) - 2)*L^2/(36*pi^2))*sin(3*pi*x/L);

subdomain_collocation = (2*L^2/(9*pi))*sin(pi*x/L) ...
                      + (L^2/(18*pi))*sin(2*pi*x/L) ...
                      + (L^2/(108*pi))*sin(3*pi*x/L);

names = ["Galerkin"; "Point collocation"; "Subdomain collocation"];
solutions = [galerkin; point_collocation; subdomain_collocation];
errors = solutions - exact;
metrics = table(names, sqrt(mean(errors.^2, 2)), max(abs(errors), [], 2), ...
    'VariableNames', {'Method', 'RMSError', 'MaxAbsoluteError'});

fig = figure('Color', 'w', 'Visible', 'off');
plot(x/L, exact, 'k-', 'LineWidth', 2.2); hold on;
plot(x/L, galerkin, '--', 'LineWidth', 1.7);
plot(x/L, point_collocation, '-.', 'LineWidth', 1.7);
plot(x/L, subdomain_collocation, ':', 'LineWidth', 2.0);
grid on; box on;
xlabel('x/L'); ylabel('\phi(x)');
title('Discontinuous source: three-term weighted-residual approximations');
legend('Analytical', 'Galerkin', 'Point collocation', ...
    'Subdomain collocation', 'Location', 'best');
exportgraphics(fig, fullfile(output_dir, 'weighted_residual_step_source.png'), ...
    'Resolution', 200);
close(fig);
end
