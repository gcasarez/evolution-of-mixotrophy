#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Parameters*)


pW = 3;
gammaP = 2.9/3;
sigmaP = 0;
kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;
hW = 50;
hP = hP2 = 50;
lW = 0.5;
lH = lH2 = 0.1;
lP = lP2 = 0.4;
dMax = 0.3;
sigmaD = 1; (* New! *)
a = a2 = 0.15;
gammaA = 0; (* We want this to be 0 consistently, because we've decided the phototroph doesn't eat *)
sigmaA = 0; (* 0 for no attack rate trade-off, 1 for attack rate trade-off *)
c = 0.1;
gammaC = 1; (* 1 for no conversion efficiency trade-off, 0.7 for conversion efficiency trade-off *)
sigmaC = 0; (* Added sigmaC because it looks like it's not defined here *)
iIn = 25; (* Incident light *)


(* ::Section:: *)
(*Equations*)


(* ::Text:: *)
(*Prey*)


wDot[w_, hW_, lW_, a_, pW_, kappa_, iIn_, gammaA_, eta_, phi_, eta2_, phi2_, sigmaA_, p_, h_, a2_, h2_, p2_]:= 
(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-(a(1-sigmaA (phi eta))) w h-(gammaA a(1-sigmaA (phi eta))) w p-(a2(1-sigmaA (phi2 eta2))) w h2 - (gammaA a2(1-sigmaA (phi2 eta2))) w p2 


(* ::Text:: *)
(*Resident heterotroph*)


hDot[eta_, c_, gammaC_, kappa_, gammaP_, pW_, sigmaP_, p_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, lH_, dMax_, sigmaD_]:=
(1-eta) (c/(gammaC c))((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi (a(1-sigmaA (phi eta))) w h + (1-phi)c (a(1-sigmaA (phi eta))) w h+ (1-eta) c (gammaA a(1-sigmaA (phi eta))) w p -lH h + (dMax (1-sigmaD eta))p


(* ::Text:: *)
(*Resident phototroph*)


pDot[eta_, gammaP_, sigmaP_, pW_, p_, kappa_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, gammaC_, sigmaC_, c_, dMax_, sigmaD_, lP_] :=
eta ((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi (a(1-sigmaA (phi eta))) w h +eta (gammaA a(1-sigmaA (phi eta))) (gammaC c) p w - (dMax (1-sigmaD eta)) p - lP p


(* ::Text:: *)
(*Mutant heterotroph*)


h2Dot[eta2_, c_, gammaC_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, lH2_, dMax_, sigmaD_] :=
(1-eta2) (c/(gammaC c))((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 (a2(1-sigmaA (phi2 eta2))) w h2 + (1-phi2) c (a2(1-sigmaA (phi2 eta2))) w h2 +(1-eta2) c (gammaA a2(1-sigmaA (phi2 eta2))) w p2 - lH2 h2+ (dMax(1-sigmaD eta2)) p2


(* ::Text:: *)
(*Mutant phototroph*)


p2Dot[eta2_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, c_, gammaC_, sigmaC_, dMax_, sigmaD_, lP2_] :=
eta2 ((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 (a2(1-sigmaA (phi2 eta2))) w h2+eta2 (gammaA a2(1-sigmaA (phi2 eta2))) (gammaC c) p2 w-(dMax(1-sigmaD eta2)) p2-lP2 p2


(* ::Text:: *)
(*Light absorbance*)


kappa[kW_, kH_, kP_, kH2_, kP2_, w_, h_, p_, h2_, p2_]:=
kW * w + kH * h + kP * p + kH2 * h2 + kP2 * p2 


(* ::Section:: *)
(*Finding Equilibria*)


equilibrium[light_, retention_, replication_, guess_:{6, 1, 1}] := 
FindRoot[{0 == wDot[w, hW,lW, a, pW, k, light, gammaA, replication, retention, 0, 0, sigmaA, p, h, a2, 0, 0],
	0 == hDot[replication, c, gammaC, k, gammaP, pW, sigmaP, p, hP, light, retention, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	0 == pDot[replication, gammaP, sigmaP, pW, p, k, hP, light, retention, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}/.
	{k -> kappa[kW, kH, kP, 0, 0, w, h, p, 0, 0]}, 
	{w, guess[[1]]}, {h, guess[[2]]}, {p, guess[[3]]}, AccuracyGoal -> 6, PrecisionGoal -> 6]


(* Do not use! Slooooooooooow, maybe not something Mathematica can solve numerically *)
(*equilibria[light_, retention_, replication_] := 
NSolve[{0 \[Equal] wDot[w, hW,lW, a, pW, k, light, gammaA, replication, retention, 0, 0, sigmaA, p, h, a2, 0, 0],
	0 \[Equal] hDot[replication, c, gammaC, k, gammaP, pW, sigmaP, p, hP, light, retention, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	0 \[Equal] pDot[replication, gammaP, sigmaP, pW, p, k, hP, light, retention, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}/.
	{k \[Rule] kappa[kW, kH, kP, 0, 0, w, h, p, 0, 0]}]*)


jacobian[light_, retention_, replication_] := 
D[{wDot[w, hW,lW, a, pW, k, light, gammaA, replication, retention, 0, 0, sigmaA, p, h, a2, 0, 0],
	hDot[replication, c, gammaC, k, gammaP, pW, sigmaP, p, hP, light, retention, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	pDot[replication, gammaP, sigmaP, pW, p, k, hP, light, retention, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}/.
	{k -> kappa[kW, kH, kP, 0, 0, w, h, p, 0, 0]}, 
	{{w, h, p}}]


equilibrium[0, 0, 0, {6, 0, 0}]


Table[equilibrium[light, 0, 0, {2, 2, 0}], {light, 0, 30}]


ListLinePlot[{w/.%54, h/.%54}, PlotLegends -> {"Prey", "Heterotroph"}]


equilibria = Table[{light, equilibrium[light, 0, 0, {4, 4, 0}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]]]


(* ::Text:: *)
(*Get the Jacobian for each light level, evaluated at equilibrium*)


jacobianVsLight = Map[(jacobian[#[[1]], 0, 0] /. #[[2]])&, equilibria];


(* ::Text:: *)
(*Find the eigenvalues of the Jacobian at equilibrium*)


eigenvaluesVsLight = Map[Eigenvalues, jacobianVsLight];


(* ::Text:: *)
(*Check if the real parts of eigenvalues are negative*)


stableEquil = Map[AllTrue[#, Negative]&, Re[eigenvaluesVsLight]];


?AccuracyGoal


(* ::Subsection:: *)
(*Phototroph*)


equilibria = Table[{light, equilibrium[light, 0, 1, {1, 1, 1}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]], PlotRange -> All]


(* ::Text:: *)
(*Get the Jacobian for each light level, evaluated at equilibrium*)


jacobianVsLight = Map[(jacobian[#[[1]], 0, 1] /. #[[2]])&, equilibria];


(* ::Text:: *)
(*Find the eigenvalues of the Jacobian at equilibrium*)


eigenvaluesVsLight = Map[Eigenvalues, jacobianVsLight];


(* ::Text:: *)
(*Check if the real parts of eigenvalues are negative*)


stableEquil = Map[AllTrue[#, Negative]&, Re[eigenvaluesVsLight]];


(* ::Subsubsection:: *)
(*Prey only equilibrium*)


equilibria = Table[{light, equilibrium[light, 0, 1, {50, 0, 0}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]], PlotRange -> All]


(* ::Text:: *)
(*Get the Jacobian for each light level, evaluated at equilibrium*)


jacobianVsLight = Map[(jacobian[#[[1]], 0, 1] /. #[[2]])&, equilibria];


(* ::Text:: *)
(*Find the eigenvalues of the Jacobian at equilibrium*)


eigenvaluesVsLight = Map[Eigenvalues, jacobianVsLight];


(* ::Text:: *)
(*Check if the real parts of eigenvalues are negative*)


stableEquil = Map[AllTrue[#, Negative]&, Re[eigenvaluesVsLight]];


Table[{i, stableEquil[[i + 1]]}, {i, 0, 50}] // Transpose // MatrixForm


(* ::Subsubsection:: *)
(*Prey and heterotroph equilibrium*)


equilibria = Table[{light, equilibrium[light, 0, 1, {4, 4, 0}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]], PlotRange -> All]


(* ::Text:: *)
(*Get the Jacobian for each light level, evaluated at equilibrium*)


jacobianVsLight = Map[(jacobian[#[[1]], 0, 1] /. #[[2]])&, equilibria];


(* ::Text:: *)
(*Find the eigenvalues of the Jacobian at equilibrium*)


eigenvaluesVsLight = Map[Eigenvalues, jacobianVsLight];


(* ::Text:: *)
(*Check if the real parts of eigenvalues are negative*)


stableEquil = Map[AllTrue[#, Negative]&, Re[eigenvaluesVsLight]];


Table[{i, stableEquil[[i + 1]]}, {i, 0, 50}] // Transpose // MatrixForm


(* ::Subsubsection:: *)
(*Phototroph equilibrium*)


(* ::Text:: *)
(*Heterotroph present due to plastid decay*)


equilibria = Table[{light, equilibrium[light, 0, 0, {5, 1, 2}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]], PlotRange -> All]


(* ::Text:: *)
(*Get the Jacobian for each light level, evaluated at equilibrium*)


jacobianVsLight = Map[(jacobian[#[[1]], 0, 1] /. #[[2]])&, equilibria];


(* ::Text:: *)
(*Find the eigenvalues of the Jacobian at equilibrium*)


eigenvaluesVsLight = Map[Eigenvalues, jacobianVsLight];


(* ::Text:: *)
(*Check if the real parts of eigenvalues are negative*)


stableEquil = Map[AllTrue[#, Negative]&, Re[eigenvaluesVsLight]];


equilibria[[23]]


Table[{i, stableEquil[[i + 1]]}, {i, 0, 50}] // Transpose // MatrixForm


equilibria = Table[{light, equilibrium[light, 0, 1, {1, 4, 2}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]], PlotRange -> All]


sigmaD=1;


(* ::Subsubsection:: *)
(*No organisms present*)


equilibria = Table[{light, {w -> 0, h -> 0, p -> 0}}, {light, 0, 50}]


jacobianVsLight = Map[(jacobian[#[[1]], 0, 1] /. #[[2]])&, equilibria];


jacobianVsLight


D[wDot[w,hW,lW,a,pW,kW * w,iIn,gammaA,1,0,0,0,sigmaA,0,0, a2,0,0], {w, 1}]/. w-> 0


wDot[w,hW,lW,a,pW,kW * w,iIn,gammaA,1,0,0,0,sigmaA,0,0, a2,0,0]


(* ::Subsection:: *)
(*Obligate kleptoplast*)


sigmaD = 0;


equilibria = Table[{light, equilibrium[light, 1, 0, {1, 4, 1}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]], PlotRange -> All]


(* ::Text:: *)
(*Get the Jacobian for each light level, evaluated at equilibrium*)


jacobianVsLight = Map[(jacobian[#[[1]], 0, 1] /. #[[2]])&, equilibria];


(* ::Text:: *)
(*Find the eigenvalues of the Jacobian at equilibrium*)


eigenvaluesVsLight = Map[Eigenvalues, jacobianVsLight];


(* ::Text:: *)
(*Check if the real parts of eigenvalues are negative*)


stableEquil = Map[AllTrue[#, Negative]&, Re[eigenvaluesVsLight]];


Table[{i, stableEquil[[i + 1]]}, {i, 0, 50}] //Transpose // MatrixForm


(* ::Subsubsection:: *)
(*Prey alone*)


equilibria = Table[{light, equilibrium[light, 1, 0, {20, 0, 0}]}, {light, 0, 50}];
ListLinePlot[{w /. equilibria[[All, 2]], h /. equilibria[[All, 2]], p /. equilibria[[All, 2]]}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}, DataRange -> MinMax[equilibria[[All, 1]]], PlotRange -> All]


(* ::Text:: *)
(*Get the Jacobian for each light level, evaluated at equilibrium*)


jacobianVsLight = Map[(jacobian[#[[1]], 0, 1] /. #[[2]])&, equilibria];


(* ::Text:: *)
(*Find the eigenvalues of the Jacobian at equilibrium*)


eigenvaluesVsLight = Map[Eigenvalues, jacobianVsLight];


(* ::Text:: *)
(*Check if the real parts of eigenvalues are negative*)


stableEquil = Map[AllTrue[#, Negative]&, Re[eigenvaluesVsLight]];


Table[{i, stableEquil[[i + 1]]}, {i, 0, 50}] //Transpose // MatrixForm


Table[residentEcoDynamics[1, 0, equilibria[[i, 1]], (w /. equilibria[[i, 2]]) + 0.001, (h /. equilibria[[i, 2]]) + 0.001, (p /. equilibria[[i, 2]]) + 0.001], {i, {17}}][[1]]


Plot[{w[t] /. %[[1]], h[t] /. %[[1]], p[t] /. %[[1]]}, {t, 0, 4000}]


(* ::Subsection:: *)
(*Kleptoplastidic phototroph*)


(* ::Subsection:: *)
(*Function to find the resident ecological dynamics, with no mutant present*)


eigenvaluesVsLight // MatrixForm
