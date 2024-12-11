BeginPackage["doubleQWExperimentalModel`"];

(* Unprotect and clear all package symbols *)
Unprotect @@ Names["doubleQWExperimentalModel`*"];
ClearAll @@ Names["doubleQWExperimentalModel`*"];

(* define public symbols *)
experimentalDoubleQwEvolution;


Begin["`Private`"];

(* WE ONLY USE THE LR BASIS IN THIS HOUSE *)

walkerSpaceDimension = 5;
maxOccupationNumber = (walkerSpaceDimension - 1)/2;
(* these are the indices corresponding to the walker starting in position "0" *)
inputStateIndices = {walkerSpaceDimension, walkerSpaceDimension + 1};

(* the weird way to chose the element with 1 is to convert between notation with (-2, -1, 0, 1, 2) as indices to the natural Mathematica one with (1, 2, 3, 4, 5) *)
(* polarization convention is here R->{1, 0}, L->{0, 1} *)
base["R", n_Integer] := Normal @ Flatten @ KroneckerProduct[SparseArray[{n + maxOccupationNumber + 1 -> 1}, walkerSpaceDimension], {1, 0}];
base["L", n_Integer] := Normal @ Flatten @ KroneckerProduct[SparseArray[{n + maxOccupationNumber + 1 -> 1}, walkerSpaceDimension], {0, 1}];


(* define qplate unitary *)
QP[alpha_, delta_] := Plus[
    Cos[delta / 2] * KroneckerProduct[IdentityMatrix @ walkerSpaceDimension, IdentityMatrix @ 2],
    I * Sin[delta / 2] * Sum[
        (# + ConjugateTranspose @ #) &@ (Exp[2 I alpha] KroneckerProduct[base["L", n], Conjugate @ base["R", n + 1] ]),
        {n, -maxOccupationNumber, maxOccupationNumber - 1}
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

(* define the unitary describing A SINGLE quantum walk (without specifying input state or final projection) *)
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
] // Chop;


(* Protect all package symbols *)
With[{syms = Names["doubleQWExperimentalModel`*"]},
  SetAttributes[syms, {Protected, ReadProtected}]
];

End[];  (* End `Private` *)

EndPackage[];
