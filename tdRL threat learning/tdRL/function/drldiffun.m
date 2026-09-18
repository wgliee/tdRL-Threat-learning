function derivative = drldiffun(trialCoordinate, state, parameters)
%DRLDIFFUN Trial-wise discount reinforcement-learning update rule.
%   dV/dt = (1 - phi(t)) * alpha * (lambda - V), where
%   phi(t) = exp(-(t - mu)^2 / (2 * sigma^2)).
%
%   PARAMETERS(1) is lambda, the asymptotic freezing level.
%   PARAMETERS(2) is alpha, the learning rate.
%   PARAMETERS(3) is mu, the center of the discount function.
%   PARAMETERS(4) is sigma, its positive width.

    lambda = parameters(1);
    alpha = parameters(2);
    mu = parameters(3);
    sigma = parameters(4);

    discountWeight = exp(-((trialCoordinate - mu).^2) ./ (2 * sigma^2));
    derivative = (1 - discountWeight) .* alpha .* (lambda - state);
end
