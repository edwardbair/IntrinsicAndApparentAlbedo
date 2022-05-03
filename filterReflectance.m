function [outputRefl,varargout] = filterReflectance(dt,refl,maxValue)
% [outputRefl] = filterReflectance(dt,refl,maxValue)
%filter reflectance to eliminate unlikely values
% dt - datetime of input values, probably not continuous
% refl - input reflectance
% maxValue - maximum possible value
%
% outputRefl - filtered for outliers and smoothed

% sort by time, probably okay but just in case
[dt,ia] = sort(dt);
refl = refl(ia);
% if less than N values, just filter, no interpolation or smoothing
N = 10;
if length(dt)<N
    t = refl<0 | refl>maxValue | isnan(refl);
    outputRefl = refl(~t);
    dt = dt(~t);
else % filtering and smoothing
    % datenums for the input datetimes
    dval = datenum(dt);
    % continuous dval with no gaps
    diffval = sort(diff(dval));
    xval = dval(1):diffval(1):dval(end);
    % eliminate outliers
    t = refl<0 | refl>1 | isnan(refl);
    % smoothing spline through the input
    F = fit(dval(~t),refl(~t),'smoothingspline');
    % SLM fit across continuous values
    slm = slmengine(xval,F(xval),'knots',round(length(dt)/2),...
        'minvalue',max(0,min(refl(~t))),'maxvalue',min(maxValue,max(refl(~t))));
    % output
    outputRefl = slmeval(dval,slm);
end
if nargout>1
    varargout{1} = dt;
end

end