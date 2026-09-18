function predictions = RtoODE_matrix(parameters, trialGrid, initialValues, numberOfAnimals)
%RTOODE_MATRIX Simulate the tdRL model for all animals.
%   PREDICTIONS = RTOODE_MATRIX(PARAMETERS, TRIALGRID, INITIALVALUES,
%   NUMBEROFANIMALS) returns a number-of-trials-by-number-of-animals matrix.
%
%   PARAMETERS = [lambda, alpha, mu, sigma].
%   TRIALGRID is the interpolated model coordinate vector.
%   INITIALVALUES contains one first-trial freezing value per animal.

    numberOfTimePoints = numel(trialGrid);
    predictions = zeros(numberOfTimePoints, numberOfAnimals);

    for animalIndex = 1:numberOfAnimals
        solution = ode45(@(t, y) drldiffun(t, y, parameters), ...
            trialGrid, initialValues(animalIndex));
        predictions(:, animalIndex) = deval(solution, trialGrid);
    end
end
