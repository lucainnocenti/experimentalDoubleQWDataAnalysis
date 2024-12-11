BeginPackage["doubleQWExperimentalModel`"];


experimentalDoubleQwEvolution;


Begin["`Private`"];

(* WE ONLY USE THE LR BASIS IN THIS HOUSE *)

walkerSpaceDimension = 5;
maxOccupationNumber = (walkerSpaceDimension - 1)/2;
(* these are the indices corresponding to the walker starting in position "0" *)
inputStateIndices = {walkerSpaceDimension, walkerSpaceDimension + 1};

base[R, -2] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
base[R, -1] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
base[R, 0] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
base[R, 1] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
base[R, 2] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
base[R, 3] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

base[L, -3] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
base[L, -2] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
base[L, -1] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
base[L, 0] = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
base[L, 1] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
base[L, 2] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1};


projectionOperator[vec1_List] := KroneckerProduct[vec1, Conjugate@vec1];
projectionOperator[vec1_List, vec2_List] := KroneckerProduct[vec1, Conjugate@vec2];

QP[alpha_, delta_] := Plus[
    Sum[
        Cos[delta/2] (projectionOperator@base[L, n] + projectionOperator@base[R, n]),
        {n, -maxOccupationNumber, maxOccupationNumber}
    ],
    Sum[
        I * Sin[delta / 2] (E^(2 I alpha) projectionOperator[base[L, n], base[R, n + 1]] + E^(-2 I alpha) projectionOperator[base[R, n], base[L, n - 1]]),
        {n, -maxOccupationNumber, maxOccupationNumber}
    ]
];

HWP[theta_] := KroneckerProduct[
    IdentityMatrix @ walkerSpaceDimension,
    {
        {0, Exp[2 I theta]},
        {Exp[-2 I theta], 0}
    }
];

QWP[varphi_] := KroneckerProduct[
    IdentityMatrix @ walkerSpaceDimension, 
    {
        {1 + I, (1 - I) Exp[2 I varphi]},
        {(1 - I) Exp[-2 I varphi], 1 + I}
    } / 2
];

(* RLtoHV = KroneckerProduct[
    IdentityMatrix @ walkerSpaceDimension, 
    {
        {1, 1},
        {-I, I}
    } / Sqrt @ 2
]; *)

qwUnitary[{alpha2_, delta2_ : Pi}, {varphi2_, theta2_, zeta2_}, {alpha1_, delta1_ : Pi/2}] := Dot[
    QP[alpha2, delta2]  (* second qplate) *),
    QWP[varphi2] . HWP[theta2] . QWP[zeta2], (* coin operation *)
    QP[alpha1, delta1] (* first qplate *)
];


projectUnitaryOnState[unitary_, polarizationState_] := Dot[
    (* first space is walker, second space is polarization *)
    KroneckerProduct[IdentityMatrix @ walkerSpaceDimension, ConjugateTranspose @ polarizationState],
    unitary
];

(* project the polarization on the specified polarization state (input space is left unchanged) *)
qwUnitaryAfterPolarizationProjection[{varphip_, thetap_}, {alpha2_, delta2_ : Pi}, {varphi2_, theta2_, zeta2_}, {alpha1_, delta1_ : Pi/2}] := projectUnitaryOnState[
    qwUnitary[{alpha2, delta2}, {varphi2, theta2, zeta2}, {alpha1, delta1}],
    {Cos[thetap], Sin[thetap] * Exp[I * varphip]}
];

(* project also on the possible input states *)
qwProjectedIsometry[{varphip_, thetap_}, {alpha2_, delta2_ : Pi}, {varphi2_, theta2_, zeta2_}, {alpha1_, delta1_ : Pi/2}] := qwUnitaryAfterPolarizationProjection[
    {varphip, thetap},
    {alpha2, delta2},
    {varphi2, theta2, zeta2},
    {alpha1, delta1}
][[All, inputStateIndices]];

(* qwIsometry[{varphip_, thetap_}, {alpha2_, delta2_ : Pi}, {varphi2_, theta2_, zeta2_}, {alpha1_, delta1_ : Pi/2}] := qwUnitary[
    {varphip, thetap}, (* projection angles *)
    {alpha2, delta2}, (* second qplate *)
    {varphi2, theta2, zeta2}, (* first (and only) coin operation *)
    {alpha1, delta1} (* first qplate *)
][[ ;; ;; 2, inputStateIndices]]; *)

experimentalDoubleQwEvolution = KroneckerProduct[
    With[{
            alpha1 = 19*Pi/180, delta1 = Pi/2,
            alpha2 = 77*Pi/180, delta2 = Pi,
            zeta2 = 2.908160255954825,
            varphi2 = 1.0070977275570985,
            theta2 = 1.5570696348907915,
            thetap = 0.7853981644358207, varphip = 0.6904379401618504
        },
        qwProjectedIsometry[{varphip, thetap}, {alpha2, delta2}, {varphi2, theta2, zeta2}, {alpha1, delta1}]
    ],
    With[{
            alpha1 = 336*Pi/180, delta1 = Pi/2,
            alpha2 = 163*Pi/180, delta2 = Pi,
            zeta2 = 2.89289820797498,
            varphi2 = 1.095755783171672,
            theta2 = 1.5937676381596888,
            thetap = 0.7853981643272621, varphip = 0.6567152340642829
        },
        qwProjectedIsometry[{varphip, thetap}, {alpha2, delta2}, {varphi2, theta2, zeta2}, {alpha1, delta1}]
    ]
];


End[];  (* End `Private` *)

EndPackage[];
