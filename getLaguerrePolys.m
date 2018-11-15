function laguerreFuncs = getLaguerrePolys(length,maxOrder,a)
    
    x = zeros(length,1);
    x(1) = 1;

    laguerreFuncs = zeros(length, maxOrder);
    dt = 1;
    % Apply the recursion relations for the zeroth-order Laguerre function at time 0
    laguerreFuncs(1, 1) = dt * sqrt(1-a) * x(1);

    % Apply the recursion relations for higher-order Laguerre functions at time 0
    for j = 2:maxOrder
        laguerreFuncs(1,j) = sqrt(a) * laguerreFuncs(1, j-1);
    end

    % Iterate over timepoints
    for t = 2:length
        % Apply the recursion relations for the zeroth-order Laguerre function
        laguerreFuncs(t, 1) = sqrt(a) * laguerreFuncs(t-1, 1) + dt * sqrt(1-a) * x(t);

        % Apply the recursion relations for higher-order Laguerre functions
        for j = 2:maxOrder
            laguerreFuncs(t, j) = sqrt(a) * laguerreFuncs(t-1, j) + sqrt(a) * laguerreFuncs(t, j-1) - laguerreFuncs(t-1, j-1);
        end
    end
