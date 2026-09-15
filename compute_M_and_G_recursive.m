function [M, G] = compute_M_and_G_recursive(q, params)
N  = params.N;
nd = 2*N;

L  = params.L;
r  = params.r;
g  = params.g(:);

M = zeros(nd, nd);
G = zeros(nd, 1);

% Cumulative pose at base of current section: T0 = [R0 P0; 0 1]
R0 = eye(3);
P0 = zeros(3,1);

% Store cumulative derivatives of R0 and P0 w.r.t all previous DoFs (1..2(i-1))
% We store dR0/dqk as a cell array of 3x3 blocks.
dR0 = cell(1, nd);   % only entries 1..2(i-1) are used at step i
dP0 = zeros(3, nd);  % 3 x nd, only previous columns nonzero

for i = 1:N
    % local DoF indices in global q
    idx = (2*i-1):(2*i);
    li2 = q(idx(1));
    li3 = q(idx(2));

    % NOTE: You must map your 2-DoF representation into the 3x1 l vector.
    % If your model uses l = [l1; l2; l3] with l1 constrained, adjust here.
    % Current placeholder assumes l1 = 0.
    lvec = [0; li2; li3];

    % ---- Local kinematics at CoG (xi = cog_xi(i)) ----
    xi_c = params.cog_xi(i);
    [~, Rloc_c, ploc_c] = HTM_nume(lvec, xi_c, L, r);
    [PosJ_c, RotJ_c, ~, ~] = LocalJacob_nume(lvec, xi_c, L, r); % PosJ_c: 3x2, RotJ_c: 3x6 (two 3x3 blocks)

    % ---- Local kinematics at tip (xi=1) for recursion update ----
    [~, Rloc_tip, ploc_tip] = HTM_nume(lvec, 1, L, r);
    [PosJ_tip, RotJ_tip, ~, ~] = LocalJacob_nume(lvec, 1, L, r);

    % ============================================================
    % A) CoG global Jacobian Jv_cog (3 x nd)
    % ============================================================
    % Global CoG position: Pcog = P0 + R0*ploc_c
    % Derivative wrt previous DoFs k: dP0(:,k) + dR0{k}*ploc_c
    % Derivative wrt local DoFs: R0*PosJ_c(:,local_col)
    Jv_cog = zeros(3, nd);

    % previous columns
    if i > 1
        prev = 1:(2*(i-1));
        Jv_cog(:, prev) = dP0(:, prev);
        for k = prev
            if ~isempty(dR0{k})
                Jv_cog(:, k) = Jv_cog(:, k) + dR0{k} * ploc_c;
            end
        end
    end

    % local columns
    Jv_cog(:, idx) = R0 * PosJ_c;

    % ============================================================
    % B) Optional angular Jacobian Jw_cog (3 x nd) from dR/dq
    % ============================================================
    % CoG global rotation: Rcog = R0 * Rloc_c
    Rcog = R0 * Rloc_c;

    % Build dRcog/dqk for all k (previous + local), then map to omega via vee(R^T dR/dqk).
    Jw_cog = zeros(3, nd);

    % previous DoFs: dRcog/dqk = dR0{k} * Rloc_c
    if i > 1
        prev = 1:(2*(i-1));
        for k = prev
            dRc = dR0{k} * Rloc_c;           % 3x3
            Jw_cog(:, k) = vee( Rcog.' * dRc ); % 3x1
        end
    end

    % local DoFs: dRcog/dq_local = R0 * (dRloc_c/dq_local)
    dRloc_c_1 = RotJblock(RotJ_c, 1); % 3x3 = dR/dq_local1
    dRloc_c_2 = RotJblock(RotJ_c, 2); % 3x3 = dR/dq_local2

    dRc1 = R0 * dRloc_c_1;
    dRc2 = R0 * dRloc_c_2;

    Jw_cog(:, idx(1)) = vee( Rcog.' * dRc1 );
    Jw_cog(:, idx(2)) = vee( Rcog.' * dRc2 );

    % ============================================================
    % C) Accumulate M and G for this section mass mi
    % ============================================================
    mi = params.mi(i);

    % Translational kinetic energy term
    M = M + mi * (Jv_cog.' * Jv_cog);

    % Rotational kinetic energy term (optional)
    % Provide params.Ii as either:
    %   - a cell array: params.Ii{i} is 3x3 inertia of section i about CoG (expressed in CoG frame),
    %   - or omit Ii, then rotational term is ignored.
    if isfield(params, 'Ii') && ~isempty(params.Ii)
        if iscell(params.Ii)
            Ii = params.Ii{i};
        else
            Ii = params.Ii(:,:,i);
        end
        M = M + (Jw_cog.' * Ii * Jw_cog);
    end

    % Gravity generalized force: G += Jv^T * (m g)
    G = G + Jv_cog.' * (mi * g);

    % ============================================================
    % D) Update recursion to the tip of section i: new (R0,P0,dR0,dP0)
    % ============================================================
    % Tip global pose:
    % Ptip = P0 + R0*ploc_tip
    % Rtip = R0*Rloc_tip
    Ptip = P0 + R0 * ploc_tip;
    Rtip = R0 * Rloc_tip;

    % Update dP0 and dR0 to represent derivatives of tip pose wrt global q
    % For previous doFs k: 
    %   dPtip/dqk = dP0(:,k) + dR0{k}*ploc_tip
    %   dRtip/dqk = dR0{k}*Rloc_tip
    % For local doFs:
    %   dPtip/dq_local = R0*PosJ_tip(:,local)
    %   dRtip/dq_local = R0*(dRloc_tip/dq_local)
    dPtip = dP0; % start with previous
    dRtip = dR0; % cell array, copy

    if i > 1
        prev = 1:(2*(i-1));
        for k = prev
            dPtip(:, k) = dP0(:, k) + dR0{k} * ploc_tip;
            dRtip{k}    = dR0{k} * Rloc_tip;
        end
    end

    % local dP
    dPtip(:, idx) = R0 * PosJ_tip;

    % local dR
    dRloc_tip_1 = RotJblock(RotJ_tip, 1);
    dRloc_tip_2 = RotJblock(RotJ_tip, 2);

    dRtip{idx(1)} = R0 * dRloc_tip_1;
    dRtip{idx(2)} = R0 * dRloc_tip_2;

    % advance recursion
    P0  = Ptip;
    R0  = Rtip;
    dP0 = dPtip;
    dR0 = dRtip;
end

% enforce symmetry in M
M = 0.5*(M + M.');
end
