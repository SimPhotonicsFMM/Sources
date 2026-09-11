function TfApod = Apod(BetaX,BetaY,Eps)

% APOD - Build a smooth (raised-cosine / Tukey-type) apodization window in
%        the Fourier/modal domain, applied to the truncated set of Fourier
%        harmonics used in a modal (Fourier Modal Method / RCWA-like)
%        field expansion.
%
% DESCRIPTION
%   A field is here represented by a finite (truncated) set of Fourier
%   harmonics, each harmonic k having transverse wavevector components
%   (BetaX(k),BetaY(k)). Truncating the Fourier series introduces spurious
%   oscillations (Gibbs phenomenon) when the field is reconstructed in
%   real space, especially near material discontinuities. APOD builds a
%   set of real weights TfApod(k), one per harmonic, that are equal to 1
%   for the "low-order" harmonics and smoothly roll off to 0 for the
%   highest-order harmonics near the truncation limit. Multiplying the
%   modal coefficients by TfApod before the inverse Fourier/modal
%   synthesis (see e.g. FieldD2E) damps the highest spatial frequencies
%   and therefore reduces ringing artifacts in the reconstructed field.
%
%   The window is SEPARABLE: a 1-D raised-cosine (Tukey) window is built
%   independently along x (from BetaX) and along y (from BetaY), and the
%   2-D weight applied to harmonic k is simply the product
%   TfApod(k) = TfApodX(k) * TfApodY(k).
%
%   Tukey/raised-cosine shape (per direction):
%     - FLAT region (weight = 1) for the central fraction (1-Eps) of the
%       spectrum,
%     - COSINE ROLL-OFF (from 1 down to 0) over the outer fraction Eps of
%       the spectrum on each side (highest positive and negative orders).
%
% SYNTAX
%   TfApod = Apod(BetaX,BetaY,Eps)
%
% INPUT PARAMETERS
%   BetaX  - Transverse wavevector (or Fourier order) component along x
%            for each retained harmonic [vector, one value per harmonic]
%   BetaY  - Transverse wavevector (or Fourier order) component along y
%            for each retained harmonic [vector, one value per harmonic,
%            paired index-by-index with BetaX]
%   Eps    - Dimensionless taper fraction (Tukey parameter), in [0,1]:
%            the fraction of the spectrum, on each side, over which the
%            weight is smoothly tapered from 1 to 0.
%              Eps -> 0 : almost no tapering (rectangular window,
%                         TfApod ~= 1 everywhere)
%              Eps = 1  : tapering applied across the whole spectrum
%                         (Hann-type window)
%            (a typical default value of 0.75 is kept commented below)
%
% OUTPUT PARAMETERS
%   TfApod - Apodization weights, one per harmonic (same length/shape as
%            BetaX(:)), to be multiplied element-wise onto the modal
%            field coefficients before real-space synthesis.
%
% ALGORITHM
%   1. Degenerate case: if only a single harmonic exists in each direction
%      (BetaX and BetaY both scalars), there is nothing to taper -> return
%      TfApod = 1.
%   2. For each direction (x, then y):
%      a. Normalize the wavevector/order to an angular variable spanning
%         the full bandwidth: alpha = Beta * 2*pi / range(Beta).
%      b. The flat/tapered boundary sits at alpha0 = (1-Eps)*pi, and the
%         cosine roll-off has period Per = 2*Eps*pi.
%      c. Harmonics with |alpha| beyond (1-Eps)*max(alpha) (i.e. in the
%         outer Eps fraction of the spectrum, positive or negative side)
%         get their unit weight replaced by a raised-cosine (half-Hann)
%         taper that smoothly decreases from 1 (at the flat/tapered
%         boundary) to 0 (at the highest/lowest retained order).
%   3. Combine the two 1-D windows multiplicatively:
%      TfApod = TfApodX .* TfApodY.
%
%
% VERSION HISTORY
%   Author: M. Besbes (LCF/CNRS/IOGS) 31 August 2026

% Degenerate case: a single harmonic in each direction (no truncated
% spectrum to smooth) -> no apodization needed
if isscalar(BetaX(:)) && isscalar(BetaY(:)), TfApod = 1; return; end

% Typical default taper fraction if not supplied by the caller
if nargin == 2, Eps = 0.75; end

%% --- 1-D raised-cosine (Tukey) window along x ---
TfApodX = ones(size(BetaX(:)));
alphaX = BetaX*2*pi/(max(BetaX(:))-min(BetaX(:)));  % Angular coordinate spanning the full x-bandwidth

Per = 2*Eps*pi; alpha0 = (1-Eps)*pi;   % Roll-off period and flat/tapered boundary (angular units)

% Positive-side (highest x-orders): taper weight down from 1 to 0
Px = find(alphaX>(1-Eps)*max(alphaX(:)));
TfApodX(Px) = .5+.5*cos(2*pi/Per*(alphaX(Px)-alpha0));
% Negative-side (lowest/most negative x-orders): taper weight down from 1 to 0
Px = find(alphaX<-(1-Eps)*max(alphaX(:)));
TfApodX(Px) = .5+.5*cos(2*pi/Per*(alphaX(Px)+alpha0));
%
%% --- 1-D raised-cosine (Tukey) window along y (same construction) ---
TfApodY = ones(size(BetaY(:)));
alphaY = BetaY'*2*pi/(max(BetaY(:))-min(BetaY(:)));  % Angular coordinate spanning the full y-bandwidth

% Positive-side (highest y-orders): taper weight down from 1 to 0
Py = find(alphaY>(1-Eps)*max(alphaY(:)));
TfApodY(Py) = .5+.5*cos(2*pi/Per*(alphaY(Py)-alpha0));
% Negative-side (lowest/most negative y-orders): taper weight down from 1 to 0
Py = find(alphaY<-(1-Eps)*max(alphaY(:)));
TfApodY(Py) = .5+.5*cos(2*pi/Per*(alphaY(Py)+alpha0));


% Separable 2-D apodization weight for each harmonic: product of the x-
% and y-direction taper weights
TfApod = TfApodX.*TfApodY;


end
