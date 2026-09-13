%RUN_ALL Reproduce all results for the 1D Poisson project.

project_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(project_dir, 'results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end
addpath(project_dir);

step_metrics = weighted_residual_step_source(results_dir);
step_metrics.Problem = repmat("Discontinuous source", height(step_metrics), 1);

quadratic_metrics = weighted_residual_quadratic_source(results_dir);
quadratic_metrics.Problem = repmat("Quadratic source", height(quadratic_metrics), 1);

fvm_metrics = finite_volume_step_source(results_dir);
fvm_metrics.Problem = repmat("Discontinuous source with nonzero boundary", ...
    height(fvm_metrics), 1);

summary = [step_metrics; quadratic_metrics; fvm_metrics];
summary = movevars(summary, 'Problem', 'Before', 'Method');
writetable(summary, fullfile(results_dir, 'error_summary.csv'));

assert(all(isfinite(summary.RMSError)), 'An RMS error is not finite.');
assert(all(isfinite(summary.MaxAbsoluteError)), ...
    'A maximum absolute error is not finite.');
assert(fvm_metrics.RMSError < 1e-10, ...
    'The finite-volume RMS error exceeds the verification tolerance.');

expected_files = {
    'weighted_residual_step_source.png'
    'weighted_residual_quadratic_source.png'
    'finite_volume_step_source.png'
    'finite_volume_absolute_error.png'
    'error_summary.csv'
};
for index = 1:numel(expected_files)
    assert(isfile(fullfile(results_dir, expected_files{index})), ...
        'Expected result file was not generated: %s', expected_files{index});
end

disp(summary);
fprintf('\nVerification passed. Results written to:\n%s\n', results_dir);
