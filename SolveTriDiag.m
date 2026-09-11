function [Ek,Hk,Db,Dh] = SolveTriDiag(Sb,MatS,Sh,Fb,Fh,Ib,Ih)

% SOLVETRIDIAG - Solve the block-tridiagonal system arising from the
%                z-discretization of a multilayer stack (modified block
%                Thomas algorithm), returning the Fourier-harmonic
%                coefficients of the electric and magnetic fields at
%                every z-discretization node across the Ns layers,
%                together with the diffraction-order amplitude
%                coefficients at the bottom and top boundaries.
%
% DESCRIPTION
%   The structure is a stack of Ns layers along z, sandwiched between a
%   homogeneous bottom medium and a homogeneous top medium. Each layer k
%   is uniformly subdivided into Nsub(k) sub-layers of thickness
%   dz(k) = hc(k)/Nsub(k), giving Np = sum(Nsub) sub-intervals and
%   Np+1 z-nodes overall (including the two outer boundaries and every
%   interior layer interface). At each node the modal/Fourier-harmonic
%   field coefficient vector (length m, one column per polarisation: TE
%   and TM, as in FieldD2E) is coupled to its two neighbours through a
%   three-point (tridiagonal) finite-difference relation built from each
%   layer's differential operators A_k and B_k. Stacked over all nodes
%   this yields one large block-tridiagonal linear system, which
%   SOLVETRIDIAG solves with a MODIFIED BLOCK THOMAS ALGORITHM (block
%   Gaussian elimination specialised to a tridiagonal block structure):
%   a forward sweep reduces every block-row to an equivalent one-sided
%   (upper block-bidiagonal) system, and a backward sweep recovers the
%   electric field Ek at every node from the reduced system.
%
%   Boundary conditions are enforced by folding a radiation/impedance-
%   type boundary operator - the Sb12 block extracted from Sb at the
%   bottom, the Sh21 block extracted from Sh at the top (see
%   IMPLEMENTATION NOTES) - together with a source/excitation term (Fb
%   at the bottom, Fh at the top) into the first and last block-rows of
%   the system before elimination.
%
%   Once Ek is known, the function performs two further steps:
%     1. MAGNETIC FIELD: Hk is reconstructed node by node from Ek, using
%        a central finite-difference derivative through each layer's
%        admittance-like operator ak{.} = inv(A_k) (interior nodes), and
%        a one-sided/averaged O(dz) formula at every layer interface;
%        at the two outer boundaries Hk is obtained algebraically from
%        Ek via the bottom/top boundary operators.
%     2. DIFFRACTION EFFICIENCIES: Db and Dh, the Fourier-harmonic
%        amplitude coefficients of the outgoing (diffracted/reflected)
%        orders at the bottom and top boundaries respectively, are
%        obtained from Ek at the corresponding boundary node together
%        with the known incident-order amplitudes Ib (from below) and
%        Ih (from above), via the S-matrix blocks of Sb and Sh. Squaring
%        and normalising these (by the caller) yields the actual
%        bottom/top diffraction efficiencies.
%
% SYNTAX
%   [Ek,Hk,Db,Dh] = SolveTriDiag(Sb,MatS,Sh,Fb,Fh,Ib,Ih)
%
% INPUT PARAMETERS
%   Sb    - Cell array describing the bottom (z=0 side) homogeneous
%           medium/boundary, in S-matrix form. Only the following
%           elements are used by this function:
%           • Sb{1} : partitioned scattering operator of size
%                     (m+length(Pd)) x (length(Pdi)+m), split below into
%                     4 blocks Sb11,Sb12,Sb21,Sb22 (see IMPLEMENTATION
%                     NOTES for what each block couples).
%           • Sb{6} : reference array whose row count gives m, the
%                     block size (number of Fourier harmonics per
%                     polarisation) of the tangential-field unknown at
%                     the bottom boundary.
%           • Sb{7} = Pdi : indices of the retained incident plane-wave
%                     orders at the bottom (rows of Ib to be used).
%           • Sb{8} = Pd  : indices of the retained diffracted/reflected
%                     plane-wave orders at the bottom (defines the size
%                     of Db).
%           (Sb{2}..Sb{5} are not used inside SolveTriDiag.)
%
%   MatS  - Cell array [Ns x 3], one row per layer k = 1..Ns:
%           • MatS{k,1} = A_k : layer-k differential operator that
%             multiplies the second-difference (curvature) term of the
%             z-discretization. May be given as a plain [m x m] matrix
%             or as a cell array of sub-blocks (normalised to a plain
%             matrix via cell2mat before use).
%           • MatS{k,2} = B_k : layer-k differential operator that
%             multiplies the zero-order (local/potential-like) term of
%             the z-discretization (e.g. built from the transverse
%             Fourier/modal operator and the layer's permittivity).
%             Same [m x m] or cell-array form as A_k.
%           • MatS{k,3} = Data(k) : struct describing layer k:
%               .hc    : physical thickness of layer k
%               .Nsub  : number of uniform z-sub-layers used to
%                        discretize layer k (finite-difference
%                        resolution within the layer)
%
%   Sh    - Cell array describing the top (z=zmax side) homogeneous
%           medium/boundary, mirroring Sb: Sh{1} partitioned scattering
%           operator (split below into Sh11,Sh12,Sh21,Sh22), Sh{6}
%           (unused here), Sh{7} = Pui (incident orders from the top,
%           rows of Ih to be used), Sh{8} = Pu (diffracted/transmitted
%           orders at the top, defines the size of Dh).
%
%   Fb    - Source/excitation term folded into the bottom boundary
%           block-row of the tridiagonal system, [m x 2].
%
%   Fh    - Source/excitation term folded into the top boundary
%           block-row of the tridiagonal system, [m x 2].
%
%   Ib    - Incident-field amplitude spectrum impinging from below,
%           indexed by Pdi (only Ib(Pdi,:) is used). Pass Ib such that
%           Ib(Pdi,:) is empty (e.g. Ib = []) to indicate that there is
%           no incident field from below (purely outgoing/radiating
%           bottom boundary); Db is then computed from Ek{1} alone.
%
%   Ih    - Incident-field amplitude spectrum impinging from above,
%           indexed by Pui (only Ih(Pui,:) is used). Same empty
%           convention as Ib, for "no incident field from above"; Dh is
%           then computed from Ek{Np+1} alone.
%
% OUTPUT PARAMETERS
%   Ek    - [Np+1 x 1] cell array of electric-field Fourier-harmonic
%           coefficients, Np = sum(Nsub) over the Ns layers. Ek{k} is an
%           [m x 2] array holding the harmonic coefficients at the k-th
%           z-node (k=1 at the bottom boundary, k=Np+1 at the top
%           boundary); columns 1 and 2 correspond to the TE and TM
%           polarisations respectively (same convention as FieldD2E).
%
%   Hk    - [Np+1 x 1] cell array of magnetic-field Fourier-harmonic
%           coefficients, same node/column convention as Ek.
%
%   Db    - [length(Pd) x 2] diffracted/reflected-order Fourier-harmonic
%           amplitude coefficients at the bottom boundary; the bottom
%           diffraction efficiencies are obtained from Db by the caller
%           (typically |Db|^2 normalised by the relevant propagation-
%           constant/power ratio of each order).
%
%   Dh    - [length(Pu) x 2] diffracted/transmitted-order Fourier-
%           harmonic amplitude coefficients at the top boundary; the top
%           diffraction efficiencies are obtained from Dh the same way
%           as for Db.
%
% ALGORITHM SUMMARY (modified block Thomas algorithm)
%   1. Unpack per-layer data (Data(k) <- MatS{k,3}) and compute the
%      sub-layer thicknesses dz(k) = hc(k)/Nsub(k) and the total number
%      of sub-intervals Np = sum(Nsub). Split the bottom/top S-matrices
%      Sb{1}/Sh{1} into their 4 constituent blocks (see IMPLEMENTATION
%      NOTES).
%   2. FORWARD SWEEP (bottom -> top), building at every node P a
%      reduced pivot ck{P} = inv(M1) and reduced right-hand side
%      fk{P}, exactly as in the scalar Thomas algorithm but with block
%      (matrix) pivots and right-hand sides:
%        a. Bottom boundary (node 1): combine the boundary operator
%           Sb12 with layer 1's interior operator D1 into Cb = Db+D1,
%           and the source term Fb into fk{1}.
%        b. Interior sub-layer nodes of a layer: three-point recursion
%           using D = Id + dz^2/2*(A_k*B_k), i.e. a second-order
%           (O(dz^2)) accurate finite-difference discretization of the
%           interior 1-D propagation equation within the layer.
%        c. Interface nodes between two consecutive layers k-1 and k:
%           combine both layers' operators (D1 from layer k-1, D2 from
%           layer k) into C1 = D1+D2, enforcing continuity of the field
%           and of its associated "flux" term across the material
%           interface.
%        d. Top boundary (node Np+1): combine the top boundary operator
%           Sh21 with the last layer's interior operator D2 into Ch
%           (for a single-layer stack, Ns==1, D2 is set equal to D1
%           since the layers-2..Ns loop never runs), and the source
%           term Fh into the final reduced right-hand side.
%   3. BACKWARD SUBSTITUTION (top -> bottom): starting from
%      Ek{Np+1} = -M1\fk{end}, recover the electric field at every node
%      by sweeping back down through each layer,
%        Ek{P} = -ck{P}*(fk{P} - Ek{P+1})
%      (with an extra InvA2 = ak{k}/dz factor applied exactly at each
%      layer's bottom-most node, to account for the change of local
%      operator A across the interface).
%   4. MAGNETIC FIELD RECONSTRUCTION: given Ek, compute Hk at every
%      node - central difference through ak{k} for interior nodes of a
%      layer, and an averaged one-sided O(dz) formula (combining both
%      sides' ak/B operators) at every interface node - plus the two
%      boundary values Hk{1}, Hk{Np+1} obtained algebraically from Ek
%      via the boundary S-matrix blocks.
%   5. DIFFRACTION-ORDER AMPLITUDES: assemble Db from Ek{1} (plus the
%      incident term Ib(Pdi,:) if not empty) via Sb21/Sb22/Sb12/Sb11,
%      and Dh from Ek{Np+1} (plus Ih(Pui,:) if not empty) via
%      Sh11/Sh12.
%
% IMPLEMENTATION NOTES
%   • S-matrix block partitioning. Sb{1} is partitioned as
%       Sb1 = [ Sb11  Sb12 ]      (rows 1:m        , cols 1:length(Pdi) and length(Pdi)+1:end)
%             [ Sb21  Sb22 ]      (rows m+1:end    ,  "        "                "        "  )
%     where the m columns/rows correspond to the internal tangential
%     field at the bottom node and the length(Pdi)/length(Pd) columns/
%     rows correspond to the incident/diffracted plane-wave bases:
%       - Sb12 [m x m]            : the bottom boundary "impedance"-type
%                                    operator relating internal H to
%                                    internal E (Hk{1} = Sb12\Ek{1} when
%                                    there is no incident field from
%                                    below); this is the operator folded
%                                    into the forward sweep as the
%                                    bottom boundary condition.
%       - Sb11 [m x length(Pdi)]  : couples the incident amplitudes
%                                    Ib(Pdi,:) into the internal E/H
%                                    relation.
%       - Sb21 [length(Pd) x length(Pdi)] : direct/specular coupling
%                                    from incident to diffracted orders.
%       - Sb22 [length(Pd) x m]   : couples the internal field into the
%                                    diffracted-order amplitudes.
%     Sh{1} is partitioned the same way but with the block order
%     mirrored (incident/diffracted vs. internal columns and rows are
%     swapped relative to Sb1, since the "outgoing" direction is
%     upward at the top instead of downward at the bottom):
%       - Sh21 [m x m]             : top-boundary counterpart of Sb12
%                                     (folded into the forward sweep as
%                                     the top boundary condition).
%       - Sh11 [length(Pu) x m], Sh12 [length(Pu) x length(Pui)],
%         Sh22 [m x length(Pui)]  : counterparts of Sb22, Sb21, Sb11.
%   • Hk{Np+1} (top boundary magnetic field) is computed as
%     Sh21\Ek{Np+1} (+ Sh12/Sh22 incident terms where relevant), the top
%     counterpart of Hk{1} = Sb12\Ek{1} at the bottom - both branches
%     (isempty(Ih(Pui)) or not) now consistently use the top operator
%     Sh21. An earlier version used the bottom operator Sb12 here by
%     mistake; this was caught during documentation review and fixed.
%   • MatS{k,1} and MatS{k,2} may be plain matrices or cell arrays of
%     sub-blocks; both are normalised to plain matrices via cell2mat
%     before use.
%   • Data is pre-allocated by assigning its last element first
%     (Data(Ns) = MatS{Ns,3}) and then filling Data(1:Ns-1) in a loop -
%     a common MATLAB idiom that allocates the whole struct array (with
%     its final field set) in one step instead of growing it one field
%     at a time.
%
% SEE ALSO
%   CalculFieldFD_FMM, FieldD2E
%
% VERSION HISTORY
%   Author: Mondher Besbes (LCF/CNRS/IOGS) 2026-09-10.

Ns = size(MatS,1); % Slice number
%
[Pdi,Pd,Pui,Pu] = deal(Sb{7},Sb{8},Sh{7},Sh{8});
m = size(Sb{6},1);
%
% Split the top S-matrix Sh{1} into its 4 blocks: the internal-internal
% "impedance" block Sh21 [m x m] used as the top boundary operator, and
% Sh11/Sh12/Sh22 coupling the internal field to the incident (Pui) and
% diffracted (Pu) plane-wave amplitudes - see IMPLEMENTATION NOTES.
Sh1 = Sh{1};
[n1,n2,n3,n4] = deal(m,length(Pui),length(Pu),m);
[Sh11,Sh12,Sh21,Sh22] = deal(Sh1(1:n3,1:n1),Sh1(1:n3,n1+1:end),...
                             Sh1(n3+1:end,1:n1),Sh1(n3+1:end,n1+1:end));
%
% Split the bottom S-matrix Sb{1} into its 4 blocks: the internal-
% internal "impedance" block Sb12 [m x m] used as the bottom boundary
% operator, and Sb11/Sb21/Sb22 coupling the internal field to the
% incident (Pdi) and diffracted (Pd) plane-wave amplitudes.
Sb1 = Sb{1};
[n1,n2,n3,n4] = deal(length(Pdi),m,m,length(Pd));
[Sb11,Sb12,Sb21,Sb22] = deal(Sb1(1:n3,1:n1),Sb1(1:n3,n1+1:end),...
                             Sb1(n3+1:end,1:n1),Sb1(n3+1:end,n1+1:end));
%
% Pre-allocate the per-layer struct array in one step (last element
% first), then unpack the per-layer data (thickness hc, sub-layer count
% Nsub) stored in the 3rd column of MatS into Data(1:Ns).
Data(Ns) = MatS{Ns,3};
%
for k = 1:Ns-1, Data(k) = MatS{k,3}; end

%
hc1 = cell2mat({Data.hc});     % Physical thickness of each layer, hc1(k)
Np1 = cell2mat({Data.Nsub});   % Number of z-sub-layers per layer, Np1(k)
Np = sum(Np1);                 % Total number of z-sub-intervals over the stack
%
ak = cell(Ns,1);       % Cached inv(A_k) for each layer (reused in the backward sweep and in the Hk reconstruction)
ck = cell(Np+1,1);     % Reduced (Thomas-eliminated) block pivots, one per z-node
fk = cell(Np+1,1);     % Reduced (Thomas-eliminated) block right-hand sides, one per z-node
%
dz1 = hc1(1)/Np1(1);    % Sub-layer thickness within layer 1

%% --- Bottom boundary (node 1): layer 1 operators + bottom boundary operator Sb12 ---
A1 = MatS{1,1};
B1 = MatS{1,2};
%
if iscell(A1), A1 = cell2mat(A1); end
if iscell(B1), B1 = cell2mat(B1); end
%
Id = speye(size(A1));
%
InvSb12 = Sb12\Id; %
Db = InvSb12;  % Local working variable for the bottom boundary contribution to Cb; NOT yet the output Db (which is overwritten below once Ek{1} is known)
%
ak{1} = inv(A1);
InvA1 = ak{1}/dz1;
D1 = InvA1 + dz1/2*B1;   % Layer-1 interior operator (O(dz^2) finite-difference term)
%
Cb = Db+D1;
%
M1 = Cb;
fk{1} = Fb;
ck{1} = inv(M1);
%
if Np1(1)>1, D = Id + dz1^2/2*(A1*B1);end   % Interior finite-difference operator for layer 1
%
%% --- Forward sweep: interior sub-layer nodes of layer 1 ---
P = 1;
for ks = 1:Np1(1)-1
    P = P+1;
    if ks == 1
        M1 = 2*D - ck{P-1}*InvA1;
    else
        M1 = 2*D - ck{P-1};
    end
    fk{P} = ck{P-1}*fk{P-1};
    ck{P} = inv(M1);
end

%% --- Forward sweep: layers 2..Ns (interface node + interior sub-layer nodes) ---
for k = 2:Ns
    P = P+1;
    %
    dz2 = hc1(k)/Np1(k);    % Sub-layer thickness within layer k
    %
    A2 = MatS{k,1};
    B2 = MatS{k,2};
    %
    if iscell(A2), A2 = cell2mat(A2); end
    if iscell(B2), B2 = cell2mat(B2); end
    %
    ak{k} = inv(A2);
    InvA2 = ak{k}/dz2;
    D2 = InvA2 + dz2/2*B2;   % Layer-k interior operator (O(dz^2) finite-difference term)
    C1 = D1+D2;     % Interface pivot: continuity between layer k-1 and layer k
    %
    if Np1(k)>1
        M1 = C1 - InvA1*ck{P-1};
    else
        M1 = C1 - InvA1*ck{P-1}*InvA1;
    end
    fk{P} = InvA1*(ck{P-1}*fk{P-1});
    ck{P} = inv(M1);
    %
    InvA1 = InvA2;
    D1 = D2;
    %
    if Np1(k)>1, D = Id + dz2^2/2*(A2*B2); end   % Interior finite-difference operator for layer k
    %
    % Interior sub-layer nodes of layer k
    for ks = 1:Np1(k)-1
        P = P+1;
        if ks == 1
            M1 = 2*D - ck{P-1}*InvA1;
        else
            M1 = 2*D - ck{P-1};
        end
        fk{P} = ck{P-1}*fk{P-1};
        ck{P} = inv(M1);
    end
end

%% --- Top boundary (node Np+1): last layer operators + top boundary operator Sh21 ---
Dh = -Sh21;     % Local working variable for the top boundary contribution to Ch; NOT yet the output Dh (which is overwritten below once Ek{Np+1} is known)
if Ns == 1, D2 = D1; end   % Single-layer stack: the layers-2..Ns loop above never ran, so D2 must be taken equal to D1 here
Ch = Dh+D2;     % Combined top-boundary block-row pivot
%
if Np1(Ns)>1
    M1 = Ch - InvA1*ck{P};
else
    M1 = Ch - InvA1*ck{P}*InvA1;
end

fk{P+1} = Fh + InvA1*(ck{P}*fk{P});
%
Ek = cell(Np+1,1);
Ek{Np+1} = -M1\fk{end};    % Field at the top boundary node, from the fully reduced system
%
P = Np+1;
%
%% --- Backward substitution: recover the electric field at every node, top -> bottom ---
for k = Ns:-1:1
    dz2 = hc1(k)/Np1(k);
    InvA2 = ak{k}/dz2; %inv(A2);
    for ks = Np1(k):-1:2
        P = P-1;
        Ek{P} = -ck{P}*(fk{P} - Ek{P+1}); ck{P} = [];
    end
    P = P-1;
    Ek{P} = -ck{P}*(fk{P} - InvA2*Ek{P+1}); ck{P} = [];   % Layer's bottom-most node: extra InvA2 factor across the interface
end
%

%%
% For magnetic field calculation

Hk = cell(Np+1,1);

% Bottom boundary: diffracted-order amplitudes Db and tangential Hk{1},
% from Ek{1} alone (no incident field from below), or including the
% incident contribution Ib(Pdi,:) otherwise.
if isempty(Ib(Pdi))
    Db = Sb22*(Sb12\Ek{1});
    Hk{1} = Sb12\Ek{1};
else
    Db = Sb21*Ib(Pdi,:)-Sb22*(Sb12\(Sb11*Ib(Pdi,:)))+Sb22*(Sb12\Ek{1});
    Hk{1} = Sb12\Ek{1} - Sb12\(Sb11*Ib(Pdi,:));
end
%

% Top boundary: diffracted/transmitted-order amplitudes Dh and
% tangential Hk{Np+1}, from Ek{Np+1} alone (no incident field from
% above), or including the incident contribution Ih(Pui,:) otherwise.
if isempty(Ih(Pui))
    Dh = Sh11*Ek{Np+1};
    Hk{Np+1} = Sh21\Ek{Np+1};
else
    Dh = Sh11*Ek{Np+1} + Sh12*Ih(Pui,:);
    Hk{Np+1} = Sh21\Ek{Np+1}  + Sh22*Ih(Pui,:);
end

%
% Interior nodes: central-difference Hk from Ek within each layer, and
% an averaged one-sided O(dz) formula at every layer interface (where
% the local operator ak/B switches from one layer to the next).
P = 1;
for k = 1:Ns
    dz1 = hc1(k)/Np1(k);
    %
    for ks = 1:Np1(k)-1
        P = P+1;
        Hk{P} = ak{k}*(Ek{P+1}-Ek{P-1})/(2*dz1);
    end
    %
    if P <= Np-1
        dz2 = hc1(k+1)/Np1(k+1);
        %
        B1 = MatS{k,2};
        if iscell(B1), B1 = cell2mat(B1); end
        B2 = MatS{k+1,2};
        if iscell(B2), B2 = cell2mat(B2); end
        %
        P = P+1;
        Hk{P} = 0.5*(ak{k+1}*(Ek{P+1}-Ek{P})/dz2-dz2/2*B2*Ek{P} ...
                    -ak{k}*(Ek{P-1}-Ek{P})/dz1+dz1/2*B1*Ek{P});
    end

end
%

end

