function density = myNormalPDF(x, mu, sigma)
%MYNORMALPDF Evaluate a normal probability density without a toolbox call.
%   X may be a scalar or array. MU is the center and SIGMA must be a
%   positive scalar width.

    if ~isscalar(sigma) || ~isfinite(sigma) || sigma <= 0
        error('myNormalPDF:InvalidSigma', ...
            'sigma must be a finite positive scalar.');
    end
    density = exp(-((x - mu).^2) ./ (2 * sigma^2)) ...
        ./ (sigma * sqrt(2 * pi));
end
