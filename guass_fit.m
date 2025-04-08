function [FWHM_x, FWHM_y, fitback] = guass_fit(data)
    % Gaussian function fitting for a 2D image.

    [X, Y] = meshgrid(1:size(data,2), 1:size(data,1));
    initial_guess = [max(data(:)), mean(X(:)), std(X(:)), mean(Y(:)), std(Y(:))];
    
    gaussian = fittype(@(a, x0, sigmax, y0, sigmay, x, y) ...
        a * exp(-((x-x0).^2/(2*sigmax^2) + (y-y0).^2/(2*sigmay^2))), ...
        'independent', {'x', 'y'}, 'dependent', 'z');
    fittedmodel = fit([X(:), Y(:)], data(:), gaussian, 'StartPoint', initial_guess);
    fitback = fittedmodel(X, Y);
    coeffs = coeffvalues(fittedmodel);
    A = coeffs(1);
    x0 = coeffs(2);
    sigmax = coeffs(3);
    y0 = coeffs(4);
    sigmay = coeffs(5);

    FWHM_x = abs(2 * sqrt(2 * log(2)) * sigmax);
    FWHM_y = abs(2 * sqrt(2 * log(2)) * sigmay);

end

