function predictions = rwRtoODE_matrix(parameters, trialGrid, initialValues, numberOfAnimals)
%RWRTOODE_MATRIX Simulate the R-W model for all animals.
%   PREDICTIONS = RWRTOODE_MATRIX(PARAMETERS, TRIALGRID, INITIALVALUES,
%   NUMBEROFANIMALS) returns a number-of-trials-by-number-of-animals matrix.
%
%   PARAMETERS(1) is lambda, the asymptotic freezing level.
%   PARAMETERS(2) is alpha, the learning rate.
%   TRIALGRID is the interpolated model coordinate vector.
%   INITIALVALUES contains one first-trial freezing value per animal.

    numberOfTimePoints = numel(trialGrid);
    predictions = zeros(numberOfTimePoints, numberOfAnimals);

    for animalIndex = 1:numberOfAnimals
        solution = ode45(@(~, y) rwdiffun(y, parameters), ...
            trialGrid, initialValues(animalIndex));
        predictions(:, animalIndex) = deval(solution, trialGrid);
    end
end
