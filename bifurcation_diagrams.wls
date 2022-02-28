#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Parameters*)


pW = 3;
gammaP = 2/3;(*2/3;*)(*2.9/3;*)
sigmaP = 0;
kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;
hW = 50;(*50;*)
hP = hP2 = 50;(*10;*)(*50;*)
lW = 0.5;
lH = lH2 = 0.25;(*0.5;*)(*0.1;*)
lP = lP2 = 0.25;(*0.1;*)(*0.4;*)
dMax = 0.3;(*0.3;*)
sigmaD = 0; (* No more decay trade-off *)
a = a2 = 0.15;
gammaA = 0; (* We want this to be 0 consistently, because we've decided the phototroph doesn't eat *)
(*sigmaA = 0;*) (* 0 for no attack rate trade-off, 1 for attack rate trade-off *)
c = 0.1;
(*gammaC = 1;*) (* 1 for no conversion efficiency trade-off, 0.7 for conversion efficiency trade-off *)
sigmaC = 0; (* Added sigmaC because it looks like it's not defined here *)
(*iIn = 25; *)(* Incident light *)


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


Options[equilibrium] = Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[FindRoot]];
Options[jacobian] = {attackTradeOff -> False, conversionTradeOff -> False};
Options[residentEcoDynamics] = Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[NDSolve]];


(* ::Text:: *)
(*Find an equilibrium (unknown stability, finds just one equilibrium), starting from a guess (prey, heterotroph, phototroph)*)


equilibrium[light_, retention_, replication_, guess_, opts:OptionsPattern[]] := 

Module[{
(* Implement trade-offs *)
sigmaA = If[OptionValue[attackTradeOff], 1, 0], (* 1 for attack rate trade-off, 0 for no attack rate trade-off *)
gammaC = If[OptionValue[conversionTradeOff], 0.7, 1], (* 0.7 for conversion efficiency trade-off, 1 for no conversion efficiency trade-off *)},

(* Numerically find an equilibrium, using guess as a starting point *)
FindRoot[{0 == wDot[w, hW,lW, a, pW, k, light, gammaA, replication, retention, 0, 0, sigmaA, p, h, a2, 0, 0],
	0 == hDot[replication, c, gammaC, k, gammaP, pW, sigmaP, p, hP, light, retention, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	0 == pDot[replication, gammaP, sigmaP, pW, p, k, hP, light, retention, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}/.
	{k -> kappa[kW, kH, kP, 0, 0, w, h, p, 0, 0]}, 
	
	{{w, guess[[1]]}, {h, guess[[2]]}, {p, guess[[3]]}}, 
	
	(* Optional inputs to FindRoot *)
	Evaluate[FilterRules[{opts}, Options[FindRoot]]]]]


(* Do not use! Slooooooooooow, maybe not something Mathematica can solve numerically *)
(*equilibria[light_, retention_, replication_] := 
NSolve[{0 \[Equal] wDot[w, hW,lW, a, pW, k, light, gammaA, replication, retention, 0, 0, sigmaA, p, h, a2, 0, 0],
	0 \[Equal] hDot[replication, c, gammaC, k, gammaP, pW, sigmaP, p, hP, light, retention, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	0 \[Equal] pDot[replication, gammaP, sigmaP, pW, p, k, hP, light, retention, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}/.
	{k \[Rule] kappa[kW, kH, kP, 0, 0, w, h, p, 0, 0]}]*)


(* ::Text:: *)
(*Get an expression for the Jacobian (it returns an expression in terms of prey, heterotroph, and phototroph population sizes that can be evaluated at different values later)*)


jacobian[light_, retention_, replication_, opts:OptionsPattern[]] := 
Module[{
(* Implement trade-offs *)
sigmaA = If[OptionValue[attackTradeOff], 1, 0], (* 1 for attack rate trade-off, 0 for no attack rate trade-off *)
gammaC = If[OptionValue[conversionTradeOff], 0.7, 1], (* 0.7 for conversion efficiency trade-off, 1 for no conversion efficiency trade-off *)},

D[{wDot[w, hW,lW, a, pW, k, light, gammaA, replication, retention, 0, 0, sigmaA, p, h, a2, 0, 0],
	hDot[replication, c, gammaC, k, gammaP, pW, sigmaP, p, hP, light, retention, a, w, h, gammaA, sigmaA, lH, dMax, sigmaD],
	pDot[replication, gammaP, sigmaP, pW, p, k, hP, light, retention, a, w, h, gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP]}/.
	{k -> kappa[kW, kH, kP, 0, 0, w, h, p, 0, 0]}, 
	{{w, h, p}}]]


(* ::Subsection:: *)
(*Solving the system (when we can't find a stable equilibrium*)


residentEcoDynamics[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), initialPrey_ (* W[0] *), initialHeterotroph_ (* H[0] *), initialPhototroph_ (* P[0] *), opts:OptionsPattern[]] :=

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
	
	(* Find the population sizes over 4000 days *)
	{w, h, p}, {t, 0, 4000},
	
	(* Optional input to adjust the differential equation solver *)
	FilterRules[{opts}, Options[NDSolve]]]; 
	
(* Clean up the formatting of ecoDynamics *)
ecoDynamics = Flatten[ecoDynamics, 1];

(* Check that extinction did not occur on the last time step *)
If[extinction == False,
 extinction = ((h[4000] <= 0) && (p[4000] <= 0)) /. ecoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{ecoDynamics, extinction}]]


(* ::Subsubsection:: *)
(*Heterotroph!*)


candidateStableEquils = Table[equilibrium[light, 0, 0, {guessW, guessH, guessP}, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6],
{light, 0, 60}, {guessW, 1/1000, 50, 5}, {guessH, 0, 50, 5}, {guessP, 0, 50, 5}];


candidateStableEquils = Flatten[candidateStableEquils, {{1}, {2, 3, 4}, {5}}];


duplicatesRemoved = Table[DeleteDuplicates[candidateStableEquils[[light+1]], EuclideanDistance[({w, h, p} /. #1), ({w, h, p} /. #2)] < 1/1000&], {light, 0, 60}];


jacobiansList = Table[jacobian[light, 0, 0, conversionTradeOff -> True] /. duplicatesRemoved[[light + 1]], {light, 0, 60}];
stability = Map[AllTrue[Re[Eigenvalues[#]], Negative]&, jacobiansList, {2}];
stabilityColors = Map[If[#, Directive[Blue, Opacity[0.5]], Directive[Red, Opacity[0.5]]]&, stability, {2}];


stableEquils = Pick[duplicatesRemoved, stability, True];


Dimensions[stableEquils]


heterotrophPlot = ListLinePlot[Table[{light, stateVar /. stableEquils[[light + 1, 1]]}, {stateVar, {w, h, p}}, {light, 0, 60}], PlotRange -> All,
FrameLabel -> {"Light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s)", "Equilibrium density (cells/mL)"}, Frame -> {{True, False}, {True, False}}, PlotRangePadding -> {{0, 0}, {0.05, 0.2}},
PlotLabel -> "(A) Strict heterotroph, no trade-off", PlotLegends -> {"Prey", "Heteroph", "Phototroph"}]


ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {w /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {h /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {p /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]


(* ::Subsubsection:: *)
(*Phototroph!*)


phototrophEquilibria = Table[equilibrium[light, 0, 1, {guessW, guessH, guessP}, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6, WorkingPrecision -> 18],
	{light, 0, 60}, {guessW, 1/1000, 50, 5}, {guessH, 0, 50, 5}, {guessP, 0, 50, 5}];

phototrophEquilibria = Flatten[phototrophEquilibria, {{1}, {2, 3, 4}, {5}}];


duplicatesRemoved = Table[DeleteDuplicates[phototrophEquilibria[[light+1]], EuclideanDistance[({w, h, p} /. #1), ({w, h, p} /. #2)] < 1/1000&], {light, 0, 60}];


jacobiansList = Table[jacobian[light, 0, 1, conversionTradeOff -> True] /. duplicatesRemoved[[light + 1]], {light, 0, 60}];
stability = Map[AllTrue[Re[Eigenvalues[#]], Negative]&, jacobiansList, {2}];
stabilityColors = Map[If[#, Directive[Blue, Opacity[0.5]], Directive[Red, Opacity[0.5]]]&, stability, {2}];


ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {w /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {h /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {p /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]


stableEquils = Pick[duplicatesRemoved, stability, True];
posStableEquils = Map[Select[(w >= 0) && (h >= 0) && (p >= 0) /. Chop[#]&], stableEquils];


posStableEquils == stableEquils


Map[Length, posStableEquils]


(* ::Text:: *)
(*Check whether the system goes to 0 (as expected by neighbors) at light = 8 (where solving for stable equilibria hasn't produced any equilibria that look stable)*)


residentEcoDynamics[1, 1, 8, 1, 1, 1, conversionTradeOff -> True]


posStableEquils // MatrixForm


Pick[Range[0, 60], Map[Length[#] > 1&, posStableEquils], True]


ListPlot[Table[{light, stateVar /. Last[posStableEquils[[light + 1]]]}, {stateVar, {w, h, p}}, {light, Pick[Range[0, 60], Map[Length[#] > 1&, posStableEquils], True]}]]
ListPlot[Table[{light, stateVar /. First[posStableEquils[[light + 1]]]}, {stateVar, {w, h, p}}, {light, Pick[Range[0, 60], Map[Length[#] > 1&, posStableEquils], True]}]]


(* ::Text:: *)
(*Okay, so just looking at things, I want to have kind of 3 plots?*)
(*Before the *)


phototrophPlot = Show[ListLinePlot[Table[{light, stateVar /. Last[posStableEquils[[light + 1]]]}, {stateVar, {w, h, p}}, {light, Pick[Range[0, 43], Map[Length[#] > 0&, posStableEquils[[1;;44]]], True]}], PlotRange -> {{0, 60}, {0, 20}},
FrameLabel -> {"Light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s)", "Equilibrium density (cells/mL)"}, Frame -> {{True, False}, {True, False}}, PlotRangePadding -> {{0, 0}, {0.05, 0.2}},
PlotLabel -> "(D) Strict phototroph, no trade-off", PlotLegends -> {"Prey", "Heteroph", "Phototroph"}],

ListLinePlot[Table[{light, Chop[stateVar /. posStableEquils[[light + 1, 1]]]}, {stateVar, {w, h, p}}, {light, Pick[Range[23, 60], Map[Length[#] > 0&, posStableEquils[[24;;61]]], True]}], PlotStyle -> Thin]]


Column[{Show[heterotrophPlot, ImageSize -> 350], Show[phototrophPlot, ImageSize -> 350]}]
Export["bifurcation_diagrams_noAT_panelsAD.pdf", %]


(* ::Subsubsection::Closed:: *)
(*Kleptoplastidic phototroph!*)


kleptoPhototrophEquilibria = Table[Check[equilibrium[light, 1, 1, {guessW, guessH, guessP}, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6], {"Did not find equil", w -> guessW, h -> guessH, p -> guessP}],
	{light, 0, 30}, {guessW, 1/1000, 50, 5}, {guessH, 0, 50, 5}, {guessP, 0, 50, 5}];

kleptoPhototrophEquilibria = Flatten[kleptoPhototrophEquilibria, {{1}, {2, 3, 4}, {5}}];


Table[Length[Select[kleptoPhototrophEquilibria[[i]], Length[#] == 4&]], {i, 1, 31}]


nonEquilsRemoved = Table[Select[kleptoPhototrophEquilibria[[i]], Length[#] == 3&], {i, 1, 31}];


duplicatesRemoved = Table[DeleteDuplicates[nonEquilsRemoved[[light+1]], EuclideanDistance[({w, h, p} /. #1), ({w, h, p} /. #2)] < 1/1000&], {light, 0, 30}];


jacobiansList = Table[jacobian[light, 1, 1, conversionTradeOff -> True] /. duplicatesRemoved[[light + 1]], {light, 0, 30}];
stability = Map[AllTrue[Re[Eigenvalues[#]], Negative]&, jacobiansList, {2}];
stabilityColors = Map[If[#, Directive[Blue, Opacity[0.5]], Directive[Red, Opacity[0.5]]]&, stability, {2}];


ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {w /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {h /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {p /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 30}], PlotRange -> All]


stableEquils = Pick[duplicatesRemoved, stability, True];
posStableEquils = Map[Select[(w >= 0) && (h >= 0) && (p >= 0) /. Chop[#]&], stableEquils];


posStableEquils // MatrixForm


Pick[Range[0, 30], Map[Length[#] > 0&, posStableEquils], False]


ListPlot[Table[{light, stateVar /. posStableEquils[[light + 1, 1]]}, {stateVar, {w, h, p}}, {light, Pick[Range[0, 30], Map[Length[#] > 0&, posStableEquils], True]}], PlotRange -> All,
FrameLabel -> {"Light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s)", "Equilibrium density (cells/mL)"}, Frame -> {{True, False}, {True, False}}, PlotRangePadding -> {{0, 0}, {0.05, 0.2}},
PlotLabel -> "Kleptoplastidic phototroph, no attack trade-off", PlotLegends -> {"Prey", "Heteroph", "Phototroph"}]


(* ::Text:: *)
(*We need to check what's going on with these light levels. Is there cycling here? An equilibrium we just didn't pick up?*)


Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True]


noEquilDynamics = Table[residentEcoDynamics[1, 1, light, 1, 1, 1, conversionTradeOff -> True], {light, Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True]}]


posStableEquils


(* ::Text:: *)
(*Extinction at light = 6*)


plotPoints = posStableEquils;
plotPoints[[6 + 1]] = {{w -> 0, h -> 0, p -> 0}};


Map[Plot[{w[t] /. #[[1]], h[t] /. #[[1]], p[t] /. #[[1]]}, {t, 0, 4000}]&, Select[noEquilDynamics, #[[2]] == False&]]


(* ::Text:: *)
(*Cycling at light = 13*)


Plot[{w[t] /. %[[1]], h[t] /. %[[1]], p[t] /. %[[1]]}, {t, 0, 4000}]


plotPoints[[Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True][[2]] + 1]] = 
{{w -> NMinimize[{w[t] /. noEquilDynamics[[2, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMinimize[{h[t] /. noEquilDynamics[[2, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMinimize[{p[t] /. noEquilDynamics[[2, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]},
{w -> NMaximize[{w[t] /. noEquilDynamics[[2, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMaximize[{h[t] /. noEquilDynamics[[2, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMaximize[{p[t] /. noEquilDynamics[[2, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]}}


plotPoints[[Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True][[3]] + 1]] = 
{{w -> NMinimize[{w[t] /. noEquilDynamics[[3, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMinimize[{h[t] /. noEquilDynamics[[3, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMinimize[{p[t] /. noEquilDynamics[[3, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]},
{w -> NMaximize[{w[t] /. noEquilDynamics[[3, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMaximize[{h[t] /. noEquilDynamics[[3, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMaximize[{p[t] /. noEquilDynamics[[3, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]}}


plotPoints[[Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True][[4]] + 1]] = 
{{w -> NMinimize[{w[t] /. noEquilDynamics[[4, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMinimize[{h[t] /. noEquilDynamics[[4, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMinimize[{p[t] /. noEquilDynamics[[4, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]},
{w -> NMaximize[{w[t] /. noEquilDynamics[[4, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMaximize[{h[t] /. noEquilDynamics[[4, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMaximize[{p[t] /. noEquilDynamics[[4, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]}}


plotPoints[[Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True][[5]] + 1]] = 
{{w -> NMinimize[{w[t] /. noEquilDynamics[[5, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMinimize[{h[t] /. noEquilDynamics[[5, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMinimize[{p[t] /. noEquilDynamics[[5, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]},
{w -> NMaximize[{w[t] /. noEquilDynamics[[5, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMaximize[{h[t] /. noEquilDynamics[[5, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMaximize[{p[t] /. noEquilDynamics[[5, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]}}


kleptoPhototrophPlot = Show[ListLinePlot[Table[{light, stateVar /. plotPoints[[light + 1, 1]]}, {stateVar, {w, h, p}}, {light, 0, 12}], PlotRange -> {{0, 30}, {0, 15}},
FrameLabel -> {"Light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s)", "Equilibrium density (cells/mL)"}, Frame -> {{True, False}, {True, False}}, PlotRangePadding -> {{0, 0}, {0.05, 0.2}},
PlotLabel -> "(c) Kleptoplastidic phototroph, no attack trade-off", PlotLegends -> {"Prey", "Heteroph", "Phototroph"}],

ListLinePlot[Table[{light, Chop[stateVar /. plotPoints[[light + 1, 1]]]}, {stateVar, {w, h, p}}, {light, 17, 30}]],

ListLinePlot[Table[{light, Chop[stateVar /. plotPoints[[light + 1, 1]]]}, {stateVar, {w, h, p}}, {light, 12, 17}], PlotStyle -> Dashed],
ListLinePlot[Table[{light, Chop[stateVar /. Last[plotPoints[[light + 1]]]]}, {stateVar, {w, h, p}}, {light, 12, 17}], PlotStyle -> Dashed]]


(* ::Subsubsection:: *)
(*Obligate kleptoplast!*)


obligateKleptoEquilibria = Table[Check[equilibrium[light, 1, 0, {guessW, guessH, guessP}, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6], {"Did not find equil", w -> guessW, h -> guessH, p -> guessP}],
	{light, 0, 60}, {guessW, 1/1000, 50, 5}, {guessH, 0, 50, 5}, {guessP, 0, 50, 5}];

obligateKleptoEquilibria = Flatten[obligateKleptoEquilibria, {{1}, {2, 3, 4}, {5}}];


Table[Length[Select[obligateKleptoEquilibria[[i]], Length[#] == 4&]], {i, 1, 61}]


nonEquilsRemoved = Table[Select[obligateKleptoEquilibria[[i]], Length[#] == 3&], {i, 1, 61}];


duplicatesRemoved = Table[DeleteDuplicates[nonEquilsRemoved[[light+1]], EuclideanDistance[({w, h, p} /. #1), ({w, h, p} /. #2)] < 1/1000&], {light, 0, 60}];


jacobiansList = Table[jacobian[light, 1, 0, conversionTradeOff -> True] /. duplicatesRemoved[[light + 1]], {light, 0, 60}];
stability = Map[AllTrue[Re[Eigenvalues[#]], Negative]&, jacobiansList, {2}];
stabilityColors = Map[If[#, Directive[Blue, Opacity[0.5]], Directive[Red, Opacity[0.5]]]&, stability, {2}];


ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {w /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 60}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {h /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 60}], PlotRange -> All]
ListPlot[Table[MapThread[Style[{light, #1}, #2]&, {p /. duplicatesRemoved[[light + 1]], stabilityColors[[light + 1]]}], {light, 0, 60}], PlotRange -> All]


stableEquils = Pick[duplicatesRemoved, stability, True];
posStableEquils = Map[Select[(w >= 0) && (h >= 0) && (p >= 0) /. Chop[#]&], stableEquils];


Pick[Range[0, 60], Map[Length[#] > 0&, posStableEquils], False]


ListLinePlot[Table[{light, stateVar /. posStableEquils[[light + 1, 1]]}, {stateVar, {w, h, p}}, {light, Pick[Range[0, 60], Map[Length[#] > 0&, posStableEquils], True]}], PlotRange -> All,
FrameLabel -> {"Light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s)", "Equilibrium density (cells/mL)"}, Frame -> {{True, False}, {True, False}}, PlotRangePadding -> {{0, 0}, {0.05, 0.2}},
PlotLabel -> "(B) Obligately kleptoplastidic mixotroph, no trade-off", PlotLegends -> {"Prey", "Heteroph", "Phototroph"}]


(* ::Text:: *)
(*Everything above light = 14 seems to not have a stable equilibrium. Let's investigate further!*)


noEquilDynamics = Table[residentEcoDynamics[1, 0, light, 1, 1, 1, conversionTradeOff -> True, PrecisionGoal -> 10], {light, Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True]}]


Map[Plot[{w[t] /. #[[1]], h[t] /. #[[1]], p[t] /. #[[1]]}, {t, 0, 4000}]&, Select[noEquilDynamics, #[[2]] == False&]]


plotPoints = posStableEquils;


(* ::Text:: *)
(*Okay, let's find the max's and mins for our oscillating friends*)


plotPoints[[Pick[Range[0, 30], Map[Length[#] == 0&, posStableEquils], True] + 1]] = Table[
{{w -> NMinimize[{w[t] /. noEquilDynamics[[i, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMinimize[{h[t] /. noEquilDynamics[[i, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMinimize[{p[t] /. noEquilDynamics[[i, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]},
{w -> NMaximize[{w[t] /. noEquilDynamics[[i, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	h -> NMaximize[{h[t] /. noEquilDynamics[[i, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]], 
	p -> NMaximize[{p[t] /. noEquilDynamics[[i, 1]], 3000 <= t <= 4000}, t, Method -> "SimulatedAnnealing"][[1]]}},
{i, Length[noEquilDynamics]}]


obligateKleptoPlot = Show[ListLinePlot[Table[{light, stateVar /. plotPoints[[light + 1, 1]]}, {stateVar, {w, h, p}}, {light, 0, 14}], PlotRange -> {{0, 30}, {0, 20}},
FrameLabel -> {"Light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s)", "Equilibrium density (cells/mL)"}, Frame -> {{True, False}, {True, False}}, PlotRangePadding -> {{0, 0}, {0.05, 0.2}},
PlotLabel -> "(b) Obligate kleptoplast, no attack trade-off", PlotLegends -> {"Prey", "Heteroph", "Phototroph"}],

ListLinePlot[Table[{light, Chop[stateVar /. plotPoints[[light + 1, 1]]]}, {stateVar, {w, h, p}}, {light, 14, 30}], PlotStyle -> Dashed],
ListLinePlot[Table[{light, Chop[stateVar /. Last[plotPoints[[light + 1]]]]}, {stateVar, {w, h, p}}, {light, 14, 30}], PlotStyle -> Dashed]]


(* ::Subsection:: *)
(*Combine into plot*)


GraphicsGrid[{{heterotrophPlot[[1]], obligateKleptoPlot[[1]]}, {kleptoPhototrophPlot[[1]], phototrophPlot[[1]]}, {LineLegend[{Blue, Orange, Green}, {"Prey", "Heterotroph", "Phototroph"}], SpanFromLeft}}, ImageSize -> 800]
Export["bifurcation_diagrams_noAT.pdf", %]


(* ::Subsection:: *)
(*Phototroph*)


jacobian[50, 1, 0, conversionTradeOff -> True]/. {w -> 6, h -> 0, p -> 0}


equilibria = Table[{light, equilibrium[light, 0, 1, {6, 0, 0}, conversionTradeOff -> True]}, {light, 0, 50}];
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


stableEquil


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
