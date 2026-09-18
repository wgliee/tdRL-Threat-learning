function derivative = rwdiffun(state, parameters)
%RWDIFFUN Rescorla-Wagner learning rule in continuous trial coordinates.
%   dV/dt = alpha * (lambda - V)
%
%   PARAMETERS(1) is lambda, the asymptotic freezing level.
%   PARAMETERS(2) is alpha, the learning rate.

    lambda = parameters(1);
    alpha = parameters(2);
    derivative = alpha .* (lambda - state);
end
