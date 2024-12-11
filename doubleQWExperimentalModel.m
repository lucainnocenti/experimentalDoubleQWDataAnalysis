BeginPackage["doubleQWExperimentalModel`"];


experimentalDoubleQwEvolution;


Begin["`Private`"];

(* EVERYTHING'S WRITTEN IN THE LR BASIS HERE *)

outputDimension = 5;
maxOccupationNumber = (outputDimension - 1)/2;

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
    IdentityMatrix @ outputDimension,
    {
        {0, Exp[2 I theta]},
        {Exp[-2 I theta], 0}
    }
];

QWP[varphi_] := KroneckerProduct[
    IdentityMatrix @ outputDimension, 
    {
        {1 + I, (1 - I) Exp[2 I varphi]},
        {(1 - I) Exp[-2 I varphi], 1 + I}
    } / 2
];

RLtoHV = KroneckerProduct[
    IdentityMatrix @ outputDimension, 
    {
        {1, 1},
        {-I, I}
    } / Sqrt @ 2
];

qwUnitary[{varphip_, thetap_}, {alpha2_, delta2_ : Pi}, {varphi2_, theta2_, zeta2_}, {alpha1_, delta1_ : Pi/2}] := Dot[
    RLtoHV,
    QWP[varphip] . HWP[thetap],
    QP[alpha2, delta2],
    QWP[varphi2] . HWP[theta2] . QWP[zeta2],
    QP[alpha1, delta1]
];

qwIsometry[{varphip_, thetap_}, {alpha2_, delta2_ : Pi}, {varphi2_, theta2_, zeta2_}, {alpha1_, delta1_ : Pi/2}] := Plus[
    qwUnitary[{varphip, thetap}, {alpha2, delta2}, {varphi2, theta2, zeta2}, {alpha1, delta1}][[1 ;; ;; 2, {outputDimension, outputDimension + 1}]],
    qwUnitary[{varphip, thetap}, {alpha2, delta2}, {varphi2, theta2, zeta2}, {alpha1, delta1}][[2 ;; ;; 2, {outputDimension, outputDimension + 1}]]
] / Sqrt @ 2;

experimentalDoubleQwEvolution = KroneckerProduct[
    With[{
            alpha1 = 19*Pi/180, delta1 = Pi/2,
            alpha2 = 77*Pi/180, delta2 = Pi,
            zeta2 = 2.908160255954825,
            varphi2 = 1.0070977275570985,
            theta2 = 1.5570696348907915,
            thetap = 0.7853981644358207, varphip = 0.6904379401618504
        },
        qwIsometry[{varphip, thetap}, {alpha2, delta2}, {varphi2, theta2, zeta2}, {alpha1, delta1}]
    ],
    With[{
            alpha1 = 336*Pi/180, delta1 = Pi/2,
            alpha2 = 163*Pi/180, delta2 = Pi,
            zeta2 = 2.89289820797498,
            varphi2 = 1.095755783171672,
            theta2 = 1.5937676381596888,
            thetap = 0.7853981643272621, varphip = 0.6567152340642829
        },
        qwIsometry[{varphip, thetap}, {alpha2, delta2}, {varphi2, theta2, zeta2}, {alpha1, delta1}]
    ]
];


End[];  (* End `Private` *)

EndPackage[];
