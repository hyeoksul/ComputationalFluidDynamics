function solution = tdma_solver(lower, diagonal, upper, right_hand_side)
%TDMA_SOLVER Solve a tridiagonal system using the Thomas algorithm.

number_of_equations = numel(diagonal);
modified_upper = zeros(number_of_equations, 1);
modified_rhs = zeros(number_of_equations, 1);

pivot = diagonal(1);
assert(abs(pivot) > eps, 'Zero pivot encountered in TDMA row 1.');
modified_upper(1) = upper(1)/pivot;
modified_rhs(1) = right_hand_side(1)/pivot;

for row = 2:number_of_equations
    pivot = diagonal(row) - lower(row)*modified_upper(row-1);
    assert(abs(pivot) > eps, 'Zero pivot encountered in TDMA row %d.', row);
    if row < number_of_equations
        modified_upper(row) = upper(row)/pivot;
    end
    modified_rhs(row) = ...
        (right_hand_side(row) - lower(row)*modified_rhs(row-1))/pivot;
end

solution = zeros(number_of_equations, 1);
solution(end) = modified_rhs(end);
for row = number_of_equations-1:-1:1
    solution(row) = modified_rhs(row) ...
                  - modified_upper(row)*solution(row+1);
end
end
