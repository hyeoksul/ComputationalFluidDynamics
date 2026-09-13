function metrics = weighted_residual_quadratic_source(output_dir)
%WEIGHTED_RESIDUAL_QUADRATIC_SOURCE Compare two-term approximations.
%   Solves phi'' + 10*x^2 = 0 on [0,1] with phi(0)=0 and phi(1)=1.

if nargin < 1
    output_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

x = linspace(0, 1, 1001);
exact = -(5/6)*x.^4 + (11/6)*x;

a1_g = (20/pi^3)*(1 - 4/pi^2);
a2_g = -5/(2*pi^3);
galerkin = x + a1_g*sin(pi*x) + a2_g*sin(2*pi*x);

a1_p = 50/(9*pi^2*sqrt(3));
a2_p = -5/(6*pi^2*sqrt(3));
point_collocation = x + a1_p*sin(pi*x) + a2_p*sin(2*pi*x);

a1_s = 5/(3*pi);
a2_s = -5/(16*pi);
subdomain_collocation = x + a1_s*sin(pi*x) + a2_s*sin(2*pi*x);

names = ["Galerkin"; "Point collocation"; "Subdomain collocation"];
solutions = [galerkin; point_collocation; subdomain_collocation];
errors = solutions - exact;
metrics = table(names, sqrt(mean(errors.^2, 2)), max(abs(errors), [], 2), ...
    'VariableNames', {'Method', 'RMSError', 'MaxAbsoluteError'});

fig = figure('Color', 'w', 'Visible', 'off');
plot(x, exact, 'k-', 'LineWidth', 2.2); hold on;
plot(x, galerkin, '--', 'LineWidth', 1.7);
plot(x, point_collocation, '-.', 'LineWidth', 1.7);
plot(x, subdomain_collocation, ':', 'LineWidth', 2.0);
grid on; box on;
xlabel('x'); ylabel('\phi(x)');
title('Quadratic source: two-term weighted-residual approximations');
legend('Analytical', 'Galerkin', 'Point collocation', ...
    'Subdomain collocation', 'Location', 'best');
exportgraphics(fig, fullfile(output_dir, 'weighted_residual_quadratic_source.png'), ...
    'Resolution', 200);
close(fig);
end
