function [originalTrials, observedOriginal, trialGrid, observedInterpolated] = ...
    prepare_behavior_data(inputFile, interpolationStep)
%PREPARE_BEHAVIOR_DATA Read, validate, and interpolate freezing data.
%   The input must contain no header row. Column 1 contains strictly
%   increasing trial coordinates; columns 2 onward contain one animal per
%   column, in freezing percent. All cells must be finite numeric values.

    if ~isfile(inputFile)
        error('BehaviorData:FileNotFound', ...
            'Input data file not found: %s', inputFile);
    end

    data = readmatrix(inputFile);
    if ~ismatrix(data) || size(data, 1) < 2 || size(data, 2) < 2
        error('BehaviorData:InvalidShape', ...
            ['Input must have at least two rows and two columns: ', ...
             'trial coordinate plus at least one animal.']);
    end
    if any(~isfinite(data), 'all')
        error('BehaviorData:NonNumericOrMissing', ...
            ['All input cells must be finite numbers. Remove headers, ', ...
             'blank cells, text, NaN, and Inf values.']);
    end
    if ~isscalar(interpolationStep) || ~isfinite(interpolationStep) || ...
            interpolationStep <= 0
        error('BehaviorData:InvalidStep', ...
            'interpolationStep must be a finite positive scalar.');
    end

    originalTrials = data(:, 1);
    observedOriginal = data(:, 2:end);

    if any(diff(originalTrials) <= 0)
        error('BehaviorData:TrialOrder', ...
            'Trial coordinates in column 1 must be unique and increasing.');
    end
    if any(observedOriginal < 0 | observedOriginal > 100, 'all')
        error('BehaviorData:FreezingRange', ...
            'Freezing values must be percentages in the range [0, 100].');
    end

    numberOfSteps = round( ...
        (originalTrials(end) - originalTrials(1)) / interpolationStep);
    reconstructedEnd = originalTrials(1) + ...
        numberOfSteps * interpolationStep;
    if abs(reconstructedEnd - originalTrials(end)) > 1e-10
        error('BehaviorData:InterpolationGrid', ...
            ['The interval between first and last trial coordinates must ', ...
             'be an integer multiple of %.6g.'], interpolationStep);
    end

    trialGrid = originalTrials(1) + ...
        (0:numberOfSteps)' * interpolationStep;
    observedInterpolated = interp1(originalTrials, observedOriginal, ...
        trialGrid, 'linear');
end
