function yi = interpolateNan(y);

% function replaces nan values in y with linearly interpolated values from
% non-nan values

x = find(~isnan(y));
xi = [1:length(y)];
yi = interp1(x,y(x),xi,'linear',NaN);

