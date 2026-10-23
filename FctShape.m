function Z = FctShape(X,Y,h,EpsX,EpsY)

if nargin == 4, EpsY = 0; end

Zx = ones(size(X(:)));
Xn = X*2*pi/(max(X(:))-min(X(:)));  % Angular coordinate spanning the full x-bandwidth

Per = 2*EpsX*pi; X0 = (1-EpsX)*pi;   % Roll-off period and flat/tapered boundary (angular units)

% Positive-side (highest x-orders): taper weight down from 1 to 0
Px = find(Xn>(1-EpsX)*max(Xn(:)));
Zx(Px) = .5+.5*cos(2*pi/Per*(Xn(Px)-X0));
% Negative-side (lowest/most negative x-orders): taper weight down from 1 to 0
Px = find(Xn<-(1-EpsX)*max(Xn(:)));
Zx(Px) = .5+.5*cos(2*pi/Per*(Xn(Px)+X0));
%
% --- 1-D raised-cosine (Tukey) window along y (same construction) ---
Zy = ones(size(Y(:)));
Yn = Y*2*pi/(max(Y(:))-min(Y(:)));  % Angular coordinate spanning the full y-bandwidth

Per = 2*EpsY*pi; Y0 = (1-EpsY)*pi;   % Roll-off period and flat/tapered boundary (angular units)
% Positive-side (highest y-orders): taper weight down from 1 to 0
Py = find(Yn>(1-EpsY)*max(Yn(:)));
Zy(Py) = .5+.5*cos(2*pi/Per*(Yn(Py)-Y0));
% Negative-side (lowest/most negative y-orders): taper weight down from 1 to 0
Py = find(Yn<-(1-EpsY)*max(Yn(:)));
Zy(Py) = .5+.5*cos(2*pi/Per*(Yn(Py)+Y0));


% Separable 2-D apodization weight for each harmonic: product of the x-
% and y-direction taper weights
Z = h*reshape(Zx.*Zy,size(X));

end
