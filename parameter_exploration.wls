#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*What happens if we use the parameters measured for real mixotrophs?*)


(* ::Section::Closed:: *)
(*Parameters*)


pW = 1/3;(*1/2; from moeller et al. 2016 *)
gammaP = 2/5;
sigmaP = 0;
kW = 10^-5;
kH = kH2 = 5*10^-6;
kP = kP2 = 15*10^-6;
hW = 50;
hP = hP2 = 50;
lW = 2/10;
lH = lH2 = 3/100;
lP = lP2 = 1/100;
a = a2 = 32*10^-6; (* 32*10^-7 is not from moeller et al. 2016 but seems to help? *)
dMax=1/10; (* 1/100 is not from moeller et al. 2016 but seems to help? *)
iIn = 2000;


(* ::Section::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Function for simulating ecological dynamics for arbitrary amount of time*)


residentEcoDynamics[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), initialPrey_ (* W[0] *), initialHeterotroph_ (* H[0] *), initialPhototroph_ (* P[0] *), tmax_:4000, opts:OptionsPattern[]] :=

(* Set up a module to hold some temporary parameter values *)
Module[{
(* This variable will hold the solution to the differential equations describing the ecological dynamics *)
ecoDynamics,

(* Match input with the parameters used in equations *)
phi = retention, eta = replication, iIn = surfaceLight, 

(* Implement trade-offs *)
sigmaA = If[OptionValue[attackTradeOff], 1, 0], (* 1 for attack rate trade-off, 0 for no attack rate trade-off *)
gammaC = If[OptionValue[conversionTradeOff], 0.7, 1], (* 0.7 for conversion efficiency trade-off, 1 for no conversion efficiency trade-off *)

(* Set up a variable that will become true if the mixotroph population cannot persist *)
extinction = False},

(* Now onto the equation solving!*)
ecoDynamics = NDSolve[
	{w'[t] == wDot[w[t], hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], 0, 0], iIn, gammaA, eta, phi, 0, 0, sigmaA, p[t], h[t], a2, 0, 0],
	h'[t] == hDot[eta, c, gammaC, kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], 0, 0], gammaP, pW, sigmaP, p[t], hP, iIn, phi, a, w[t], h[t], gammaA, sigmaA, lH, dMax, sigmaD],
	p'[t] == pDot[eta, gammaP, sigmaP, pW, p[t], kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], 0, 0], hP, iIn, phi, a, w[t], h[t], gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP],
	
	(* Initial conditions for the prey, heterotroph, and phototroph population sizes *)
	w[0] == initialPrey, h[0] == initialHeterotroph, p[0] == initialPhototroph, 
	
	(* If the entire mixotroph population goes to 0, stop solving the differential equations and set the variable "extinction" to true *)
	WhenEvent[(h[t] <= 0) && (p[t] <= 0), extinction = True; "StopIntegration"],
	
	(* If any single population goes below 0, set it to 0. This is just to be safe; technically the equations shouldn't allow this, but sometimes small dips below 0 happen when solving numerically *)
	WhenEvent[w[t] < 0, w[t] -> 0], WhenEvent[h[t] < 0, h[t] -> 0], WhenEvent[p[t] < 0, p[t] -> 0]},
	
	(* Find the population sizes over tmax days *)
	{w, h, p}, {t, 0, tmax},
	
	(* Optional input to adjust the differential equation solver *)
	FilterRules[{opts}, Options[NDSolve]]]; 
	
(* Clean up the formatting of ecoDynamics *)
ecoDynamics = Flatten[ecoDynamics, 1];

(* Check that extinction did not occur on the last time step *)
If[extinction == False,
 extinction = ((h[tmax] <= 0) && (p[tmax] <= 0)) /. ecoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{ecoDynamics, extinction}]]


(* ::Section::Closed:: *)
(*Investigating ecological equilibria*)


(* ::Text:: *)
(*Which equilibria are stable? What initial state variables are a good choice for getting to nontrivial equilibria?*)


Table[{exampleEcoDynamics, extinction} = residentEcoDynamics[retention, replication, 500, 10000, 100, 1, 10000, PrecisionGoal -> 10];
{retention, replication, extinction}
(*Plot[{w[t] /. exampleEcoDynamics, h[t] /. exampleEcoDynamics, p[t] /. exampleEcoDynamics}, {t, 0, 10000}, PlotLegends \[Rule] {"Prey", "Heterotroph", "Phototroph"}]*), {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}]


ListDensityPlot[Map[{#[[1]], #[[2]], Boole[#[[3]]]}&, Flatten[%506,1]], InterpolationOrder -> 0, PlotLegends -> Automatic, FrameLabel -> {"Retention", "Replication"}]


{exampleEcoDynamics, extinction} = residentEcoDynamics[1, 0, 500, 10000, 100, 1, 10000, PrecisionGoal -> 10]


Plot[{w[t] /. exampleEcoDynamics, h[t] /. exampleEcoDynamics, p[t] /. exampleEcoDynamics}, {t, 0, 2000}, PlotLegends -> {"Prey", "Heterotroph", "Phototroph"}]


(* ::Subsection::Closed:: *)
(*Prey alone*)


FindRoot[0 == wDot[w, hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w, 0, 0, 0, 0], iIn, gammaA, 0, 0, 0, 0, sigmaA, 0, 0, a2, 0, 0], {w, 10^6}]


(* ::Subsection::Closed:: *)
(*Prey and heterotroph (no retention or replication ability)*)


(* ::Text:: *)
(*(10^5 for both will do it for initial values)*)


Table[FindRoot[{0 == wDot[w, hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w, h, 0, 0, 0], iIn, gammaA, 0, 0, 0, 0, sigmaA, 0, h, a2, 0, 0],
	0 == hDot[0, c, gammaC, kappa[kW, kH, kP, kH2, kP2, w, h, 0, 0, 0], gammaP, pW, sigmaP, 0, hP, iIn, 0, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD]}, {w, 10^i}, {h, 10^j}],
	{i, 6}, {j, 6}] //Chop // MatrixForm


(* ::Subsection::Closed:: *)
(*Prey and phototroph (retention = replication = 1)*)


(* ::Text:: *)
(*Some heterotrophs will be produced due to plastid decay*)


Table[FindRoot[{0 == wDot[w, hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], iIn, gammaA, 1, 1, 0, 0, sigmaA, p, h, a2, 0, 0],
	0 == hDot[1, c, gammaC, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], gammaP, pW, sigmaP, p, hP, iIn, 1, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	0 == pDot[1, gammaP, sigmaP, pW, p, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], hP, iIn, 1, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}
	, {w, 10^i}, {h, 0}, {p, 10^j}],
	{i, 6}, {j, 6}] //Chop // MatrixForm


(* ::Text:: *)
(*With initial h = 0, we see either only prey or only mixotroph (or extinction). Initialize one to 10^6 and the other to some other values < 10^6 (both at 10^6 produces negative numbers)*)


Table[FindRoot[{0 == wDot[w, hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], iIn, gammaA, 1, 1, 0, 0, sigmaA, p, h, a2, 0, 0],
	0 == hDot[1, c, gammaC, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], gammaP, pW, sigmaP, p, hP, iIn, 1, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	0 == pDot[1, gammaP, sigmaP, pW, p, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], hP, iIn, 1, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}
	, {w, 10^i}, {h, 10^5}, {p, 10^j}],
	{i, {5}}, {j, {6}}] //Chop // MatrixForm


(* ::Text:: *)
(*Initializing with h < 0 seems to open up some equilibria where h is negative (if  w or p < 10^6, or both are 10^6) *)


(* ::Text:: *)
(*Based on this, I propose the following initialization: w = 10^5, h = 10^5, p=10^6*)


mySolns = ParallelTable[FindRoot[{0 == wDot[w, hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], iIn, gammaA, eta, phi, 0, 0, sigmaA, p, h, a2, 0, 0],
	0 == hDot[eta, c, gammaC, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], gammaP, pW, sigmaP, p, hP, iIn, phi, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	0 == pDot[eta, gammaP, sigmaP, pW, p, kappa[kW, kH, kP, kH2, kP2, w, h, p, 0, 0], hP, iIn, phi, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}, {{w, 10^5}, {h, 10^5}, {p, p0}}], 
	{phi, 0, 1, 1/10}, {eta, 0, 1, 1/10}, {p0, {0, 10^6}}];


Chop[mySolns]//MatrixForm


MatrixPlot[Chop[h/.mySolns], PlotLegends -> Automatic]
