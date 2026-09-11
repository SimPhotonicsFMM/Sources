function Exy = FieldD2E(Em,Epxy,Mesh,Phase,Phys,TfApod,x,y,c)

% FIELDD2E - Reconstruct one Cartesian component (Ex or Ey) of the electric
%            field from its modal/Fourier representation, correctly
%            handling the discontinuity of the field component normal to a
%            material interface.
%
% DESCRIPTION
%   A direct inverse Fourier/modal synthesis of a field component from its
%   modal coefficients (Em) is only accurate where that component behaves
%   as a locally tangential field, i.e. where the medium is homogeneous
%   along the direction of that component. Wherever the component is
%   normal to an interface (the medium changes along that direction),
%   synthesizing it directly from Em is physically wrong and produces poor
%   (Gibbs-type) convergence, since Em represents a genuinely discontinuous
%   quantity there.
%
%   This function therefore proceeds strip by strip along the direction of
%   the requested component:
%     1. It scans the material grid lines perpendicular to that direction
%        and, for each resulting strip, checks whether the diagonal
%        permittivity components eps_xx and eps_yy are both uniform over
%        the strip (i.e. a single homogeneous material).
%     2. In a HOMOGENEOUS strip, the requested component is locally
%        tangential (continuous): it is reconstructed directly by inverse
%        Fourier/modal synthesis, Exy = Phase*(TfApod.*Em).
%     3. In a strip that straddles a material interface (eps_xx and/or
%        eps_yy not uniform along the scanning direction), the
%        requested component is normal to that interface and therefore
%        discontinuous. The function instead reconstructs the physically
%        continuous displacement field D = Phase*(TfApod.*(Epxy*Em)), then
%        recovers E = D / eps_local separately for each element, using the
%        LOCAL permittivity of the medium the evaluation point actually
%        lies in. This correctly reproduces the jump of the normal
%        component on each side of the interface instead of smoothing it
%        out.
%
% SYNTAX
%   Exy = FieldD2E(Em,Epxy,Mesh,Phase,Phys,TfApod,x,y,c)
%
% INPUT PARAMETERS
%   Em      - Modal/Fourier coefficients of the field to synthesize.
%             Rows 1:size(Em,1)/2 hold the coefficients used for the Ex
%             problem, rows size(Em,1)/2+1:end those used for the Ey
%             problem (selected below via ni:nf depending on c). Columns
%             correspond to the TM and TE polarisations
%   Epxy    - Fourier-factorized permittivity operator : Operator
% 	      relating the modal coefficients of E to those of
%             D = eps*E, used to reconstruct the interface-continuous
% 	      quantity D instead of the discontinuous E directly.
%   Mesh    - Mesh structure:
%             . CoorN  : node coordinates [Nnodes x 2]
%             . Cn     : element-to-node connectivity
%             . CoorV  : coordinates of the material-grid lines used to
%                        detect homogeneous strips
%                        (column 1 = x-lines, column 2 = y-lines)
%             . Nsd    : subdomain (material) index of each element
%   Phase   - Inverse Fourier/modal synthesis matrix: reconstructs field
%             values at the requested (x,y) points from modal coefficients
%   Phys    - Physical constants/material structure:
%             . CaractEps(c,c,kd) : (c,c) component of the relative
%                                    permittivity tensor of subdomain kd
%             . K0                : vacuum wavenumber
%   TfApod  - Apodization/filtering window applied to the modal
%             coefficients before synthesis (limits Gibbs ringing)
%   x, y    - Coordinates of all evaluation points; need not lie within the
%             fundamental period cell, they are folded back into it
%             (modulo the cell size) before use, per periodic/Floquet
%             boundary conditions
%   c       - Requested field component: 1 = Ex (x-component),
%                                        2 = Ey (y-component)
%
% OUTPUT PARAMETERS
%   Exy     - [length(x) x 2] array containing the requested field
%             component (Ex if c=1, Ey if c=2).
%
% ALGORITHM SUMMARY
%   0. Fold the evaluation points (x,y) back into the fundamental periodic
%      cell (modulo the mesh extent Dx, Dy), so periodic repetitions of
%      the structure are handled transparently.
%   For each material-grid line xyi(ki) perpendicular to the requested
%   component:
%     - Gather the elements attached to that line (Pe) and the evaluation
%       points geometrically located in the corresponding strip (in).
%     - If eps_xx and eps_yy are both uniform over the strip: the
%       component is locally tangential -> direct synthesis from Em.
%     - Otherwise (eps_xx and/or eps_yy vary along the strip, i.e. the
%       requested component is normal to at least one interface inside
%       it): synthesize the continuous D field, then recover
%       E = D / (eps_cc * K0) element by element, using each element's own
%       subdomain permittivity, so that the correct discontinuity of the
%       normal component is reproduced across the interface.
%
% SEE ALSO
%   CalculFieldFMM, CalculFieldFD_FMM, Field
%
% VERSION HISTORY
%   Author: M. Besbes (LCF/CNRS/IOGS) 28 August 2026

Exy = zeros(length(x),2);
epsz = max(abs(Mesh.CoorN(:)))/1e6;  % Small geometric tolerance for coordinate matching

%
% Select the scanning direction according to the requested component:
%   - Ex (c=1) can be discontinuous across interfaces whose normal is
%     along x, so we scan strips of constant y and test whether the
%     material is homogeneous along x inside each strip.
%   - Ey (c=2) can be discontinuous across interfaces whose normal is
%     along y, so we scan strips of constant x and test whether the
%     material is homogeneous along y inside each strip.
%
if c == 1
    xyi = unique(Mesh.CoorV(:,2));      % y-coordinates of the horizontal strips
    CoorVxy = Mesh.CoorV(:,2);
    ni = 1; nf = size(Em,1)/2;          % Rows of Em holding the Ex modal coefficients
elseif c == 2
    xyi = unique(Mesh.CoorV(:,1));      % x-coordinates of the vertical strips
    CoorVxy = Mesh.CoorV(:,1);
    ni = 1+size(Em,1)/2; nf = size(Em,1); % Rows of Em holding the Ey modal coefficients
else
    error('c is equal to 1 for x-component or 2 for y-component');
end
%
% Fold the evaluation points back into the fundamental periodic cell
% (Floquet/periodic boundary conditions): the modal/Fourier representation
% is only defined over one period, so any (x,y) point located outside
% [-Dx/2,Dx/2] x [-Dy/2,Dy/2] must be wrapped modulo the cell size before
% being located in the mesh.
Dx = max(Mesh.CoorN(:,1))-min(Mesh.CoorN(:,1)); % x-Period
Dy = max(Mesh.CoorN(:,2))-min(Mesh.CoorN(:,2)); % y-Period
%
P = abs(x) > Dx/2;
x(P) = x(P) - floor((x(P)+Dx/2)/Dx)*Dx;
P = abs(y) > Dy/2;
y(P) = y(P) - floor((y(P)+Dy/2)/Dy)*Dy;
%
for ki = 1:length(xyi)
    % Elements attached to the current material-grid line xyi(ki)
    Pe = find(abs(CoorVxy-xyi(ki)) < epsz);
    Pen = unique(Mesh.Cn(Pe,:));
    %
    % Bounding box of the strip, used to select which of the evaluation
    % points (x,y) geometrically fall inside it
    xd = Mesh.CoorN(Pen,1); yd = Mesh.CoorN(Pen,2);
    xd = [min(xd) max(xd) max(xd) min(xd) min(xd)];
    yd = [min(yd) min(yd) max(yd) max(yd) min(yd)];
    in = inpolygon(x(:),y(:),xd(:),yd(:));
    %
    NumSd = unique(Mesh.Nsd(Pe));
    %
    if isscalar(unique(Phys.CaractEps(1,1,NumSd))) && ...
       isscalar(unique(Phys.CaractEps(2,2,NumSd)))
	%
        % Strip is a single homogeneous material along the scanning direction
        Exy(in,:) = Phase(in,:)*(TfApod.*Em(ni:nf,:));
    else
        % the requested component is normal to that interface
        % and therefore discontinuous. Reconstruct the physically
        % continuous displacement field D first...
	%
        Dxy = Phase*(TfApod.*(Epxy*Em(ni:nf,:)));
        for ke = 1:length(Pe)
            ie = Pe(ke);
            kd = Mesh.Nsd(ie);   % Local subdomain (material) of this element
            Pen = unique(Mesh.Cn(ie,:));
            %
            xd = Mesh.CoorN(Pen,1); yd = Mesh.CoorN(Pen,2);
            xd = [min(xd) max(xd) max(xd) min(xd) min(xd)];
            yd = [min(yd) min(yd) max(yd) max(yd) min(yd)];
            in = inpolygon(x(:),y(:),xd(:),yd(:));
            %
            % ...then recover E = D / eps_local on each side of the
            % interface, using this element's own subdomain permittivity
	    %
            Exy(in,:) = 1/Phys.CaractEps(c,c,kd)*Dxy(in,:)/Phys.K0;
        end
    end
end

end
