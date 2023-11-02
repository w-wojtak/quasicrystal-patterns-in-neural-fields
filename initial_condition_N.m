function u = initial_condition_N(X, Y, N)
    % Calculate the radial distance
    r = sqrt(X.^2 + Y.^2);

    % Define angular step size
    step = round(360 / N);

    % Generate equally spaced angles in degrees and radians
    thetas_deg = step:step:360;
    thetas_rad = deg2rad(thetas_deg);

    % Initialize the u_0 array
    u_0 = zeros(size(r, 1), size(r, 1), step);

    % Determine q_n based on N
    switch N
        case 8
            q_n = 8;
        case 10
            q_n = 5;
        case 12
            q_n = 12;
        otherwise
            error('N must be equal to 8, 10, or 12.');
    end

    for j = 1:length(thetas_deg)
        % Calculate the current angle in radians
        theta = thetas_rad(j);

        % Calculate the dot product for all points
        dot_product = X * cos(theta) + Y * sin(theta);
        k_j = reshape(dot_product, size(r, 1), size(r, 1));

        % Calculate q
        q = k_j * 2 * cos(pi / q_n);

        % Calculate u_0 for this angle
        u_0(:, :, j) = exp(1i * k_j) + exp(1i * q);
    end

    % Compute the resulting pattern
    u = real(sum(u_0, 3));
end

