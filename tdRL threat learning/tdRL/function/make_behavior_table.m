function outputTable = make_behavior_table(trialCoordinates, values)
%MAKE_BEHAVIOR_TABLE Add explicit trial and animal labels to an output matrix.

    numberOfAnimals = size(values, 2);
    animalNames = cellstr(compose('animal_%02d', 1:numberOfAnimals));
    variableNames = [{'trial_coordinate'}, animalNames];
    outputTable = array2table([trialCoordinates, values], ...
        'VariableNames', variableNames);
end
