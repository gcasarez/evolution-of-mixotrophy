#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
(*Parameters*)


(* ::Text:: *)
(*These parameters produce evolution to different ESSs at different light levels (25 vs. 1000) with conversion trade +/- attack trade off (same ESSs with either)*)
(*Removing conversion trade-off means we lose this, partly I think because of changes in ability of the mixotroph to survive/maintain decent population sizes at high retention rates*)


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


(* ::Text:: *)
(*Things I've changed from Holly's paper: photosynthetic rate of phototroph (2 -> 0.5), half-saturation light of phototroph (10 -> 50), plastid decay rate (0.3 -> 0.1)*)
(*I'm kind of wondering if the half-saturation light difference is key?*)
(*Things I want to do next:*)
(** The other plots: bifurcation diagrams, evolution simulations*)
(** Try adjusting the parameters to get closer to Holly's paper's parameters*)


(* ::Text:: *)
(*Going back to Grace's original parameters!*)


(*pW = 3;
gammaP = 2.9/3;
sigmaP = 0;
kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;
hW = 50;
hP = hP2 = 50;(*10;*)(*50;*)
lW = 0.5;
lH = lH2 = 0.5;(*0.5;*)(*0.1;*)
lP = lP2 = 0.4;(*0.1;*)(*0.4;*)
dMax = 0.3;(*0.3;*)
sigmaD = 0; (* No more decay trade-off *)
a = a2 = 0.15;
gammaA = 0; (* We want this to be 0 consistently, because we've decided the phototroph doesn't eat *)
(*sigmaA = 0;*) (* 0 for no attack rate trade-off, 1 for attack rate trade-off *)
c = 0.1;
(*gammaC = 1;*) (* 1 for no conversion efficiency trade-off, 0.7 for conversion efficiency trade-off *)
sigmaC = 0; (* Added sigmaC because it looks like it's not defined here *)
(*iIn = 25; *)(* Incident light *)*)


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


(* ::Section::Closed:: *)
(*Resident and mutant dynamics*)


(* ::Text:: *)
(*Set the default options : both trade - offs turned off. Turn these on or off when calling residentEcoDynamics, not here!*)
(*The options for NDSolve are also included, in case the differential equation solving run by residentEcoDynamics has to be adjusted*)


Options[residentEcoDynamics] = Options[mutantAndResidentEcoDynamics] = Options[averageEcoState] = Options[invasionDerivatives] = Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[NDSolve], Options[NIntegrate]];


(* ::Subsection:: *)
(*Function to find the resident ecological dynamics, with no mutant present*)


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


(* ::Subsection::Closed:: *)
(*Function to find the ecological dynamics when a resident and a mutant mixotroph are both present*)


mutantAndResidentEcoDynamics[residentRetention_ (* phi *), residentReplication_ (* eta *), mutantRetention_ (* phi2 *), mutantReplication_ (* eta2 *), surfaceLight_ (* iIn *), 
initialPrey_ (* W[0] *), initialResidentHeterotroph_ (* H[0] *), initialResidentPhototroph_ (* P[0] *),  initialMutantHeterotroph_ (* H2[0] *), initialMutantPhototroph_ (* P2[0] *), opts:OptionsPattern[]] :=

(* Set up a module to hold some temporary parameter values *)
Module[{
(* This variable will hold the solution to the differential equations describing the ecological dynamics *)
ecoDynamics,

(* Match input with the parameters used in equations *)
phi = residentRetention, eta = residentReplication, 
phi2 = mutantRetention, eta2 = mutantReplication,
iIn = surfaceLight, 

(* Implement trade-offs *)
sigmaA = If[OptionValue[attackTradeOff], 1, 0], (* 1 for attack rate trade-off, 0 for no attack rate trade-off *)
gammaC = If[OptionValue[conversionTradeOff], 0.7, 1], (* 0.7 for conversion efficiency trade-off, 1 for no conversion efficiency trade-off *)

(* Set up a variable that will become true if both the resident and mutant mixotroph populations cannot persist *)
extinction = False},

(* Quick check to make sure we're not starting with the resident and mutant extinct *)
If[(initialResidentHeterotroph == 0) && (initialResidentPhototroph == 0) && (initialMutantHeterotroph == 0) && (initialMutantPhototroph == 0),
extinction = True];

(* Now onto the equation solving!*)
ecoDynamics = NDSolve[
	{w'[t] == wDot[w[t], hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], h2[t], p2[t]], iIn, gammaA, eta, phi, eta2, phi2, sigmaA, p[t], h[t], a2, h2[t], p2[t]],
	h'[t] == hDot[eta, c, gammaC, kappa[kW, kH, kP, kH2, kP2 ,w[t], h[t], p[t], h2[t], p2[t]], gammaP, pW, sigmaP, p[t], hP, iIn, phi, a, w[t], h[t], gammaA, sigmaA, lH, dMax, sigmaD],
	p'[t] == pDot[eta, gammaP, sigmaP, pW, p[t], kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], h2[t], p2[t]], hP, iIn, phi, a, w[t], h[t], gammaA, sigmaA, gammaC, sigmaC, c, dMax, sigmaD, lP],
	h2'[t] == h2Dot[eta2, c, gammaC, gammaP, sigmaP, pW, p2[t], kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], h2[t], p2[t]], hP2, iIn, phi2, a2, w[t], h2[t], gammaA, sigmaA, lH2, dMax, sigmaD],
	p2'[t] == p2Dot[eta2, gammaP, sigmaP, pW, p2[t], kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], h2[t], p2[t]], hP2, iIn, phi2, a2, w[t], h2[t], gammaA, sigmaA, c, gammaC, sigmaC, dMax, sigmaD, lP2],
	
	(* Initial conditions for the prey, heterotroph, and phototroph population sizes *)
	w[0] == initialPrey, h[0] == initialResidentHeterotroph, p[0] == initialResidentPhototroph, h2[0] == initialMutantHeterotroph, p2[0] == initialMutantPhototroph,
	
	(* If the entire mixotroph population goes to 0, stop solving the differential equations and set the variable "extinction" to true *)
	WhenEvent[((h[t] <= 0) && (p[t] <= 0)) && ((h2[t] <= 0) && (p2[t] <= 0)), extinction = True; "StopIntegration"],
	
	(* If any single population goes below 0, set it to 0. This is just to be safe; technically the equations shouldn't allow this, but sometimes small dips below 0 happen when solving numerically *)
	WhenEvent[w[t] < 0, w[t] -> 0], WhenEvent[h[t] < 0, h[t] -> 0], WhenEvent[p[t] < 0, p[t] -> 0], WhenEvent[h2[t] < 0, h2[t] -> 0], WhenEvent[p2[t] < 0, p2[t] -> 0]},
	
	(* Find the population sizes over 4000 days *)
	{w, h, p, h2, p2}, {t, 0, 4000},
	
	(* Optional input to adjust the differential equation solver *)
	FilterRules[{opts}, Options[NDSolve]]]; 
	
(* Clean up the formatting of ecoDynamics *)
ecoDynamics = Flatten[ecoDynamics, 1];

(* Check that extinction did not occur on the last time step *)
If[extinction == False,
 extinction = ((h[4000] <= 0) && (p[4000] <= 0) && (h2[4000] <= 0) && (p2[4000] <= 0)) /. ecoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{ecoDynamics, extinction}]]


(* ::Subsection:: *)
(*Function to find prey ecological dynamics, no mixotrophs present*)


preyAloneEcoDynamics[surfaceLight_ (* iIn *), initialPrey_ (* W[0] *), opts:OptionsPattern[]] :=

(* Set up a module to hold some temporary parameter values *)
Module[{
(* This variable will hold the solution to the differential equations describing the ecological dynamics *)
ecoDynamics,

(* Set up a variable that will become true if the prey population cannot persist *)
extinction = False},

(* Now onto the equation solving!*)
ecoDynamics = NDSolve[
	{w'[t] == wDot[w[t], hW, lW, a, pW, kappa[kW, kH, kP, kH2, kP2, w[t], h[t], p[t], 0, 0], iIn, gammaA, eta, phi, 0, 0, sigmaA, p[t], h[t], a2, 0, 0]},
	
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


(* ::Subsection::Closed:: *)
(*Function *)


averageEcoState[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), initialPrey_ (* W[0] *), initialHeterotroph_ (* H[0] *), initialPhototroph_ (* P[0] *), opts:OptionsPattern[]] :=
Module[{ecoDynamics, extinction, averagePrey, averageHeterotroph, averagePhototroph},

(* Get the ecological dynamics *)
{ecoDynamics, extinction} = residentEcoDynamics[retention, replication, surfaceLight, initialPrey, initialHeterotroph, initialPhototroph, opts];


(* If the mixotroph went extinct, the mixotroph population sizes are 0. Maybe want to calculate the average prey population size with no mixotrophs, starting from initialPrey? Don't do this now *)
If[extinction,
	averagePrey = "Unknown";
	averageHeterotroph = averagePhototroph = 0,

(* If the mixotroph did not go extinct, calculate the average ecological state. Average the value over the last 1000 simulated time steps *)

{averagePrey, averageHeterotroph, averagePhototroph} = 
Table[NIntegrate[stateVar/. ecoDynamics, {t, 3000, 4000}, Evaluate[FilterRules[{opts}, Options[NIntegrate]]]]/(4000 - 3000),
	{stateVar, {w[t], h[t], p[t]}}]];
	
(* Return the average population sizes and whether mixotroph went extinct *)
Return[{extinction, averagePrey, averageHeterotroph, averagePhototroph}]]


invasionDerivatives[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), 
initialPrey_ (* W[0] *), initialHeterotroph_ (* H[0] *), initialPhototroph_ (* P[0] *), opts:OptionsPattern[]] :=

Module[{dMutantGrowthdRetention, dMutantGrowthdReplication},

(* Derivative w.r.t. retention *)
Module[{mutants, mutantDynamics, popSizes, derivatives},

mutants = Select[retention + {-1, 1} * 1/10000, 0 <= # <= 1&];

mutantDynamics = Map[{#, mutantAndResidentEcoDynamics[retention, replication, #, replication, surfaceLight, initialPrey, 0.99 * initialHeterotroph, 0.99 * initialPhototroph,  0.01 * initialHeterotroph,  0.01 * initialPhototroph, opts]}&,
	mutants];

(* Did introduction of any of the mutants cause an extinction? If so, the derivative of invasion fitness is maybe unclear. Set it to NA *)
If[AnyTrue[mutantDynamics[[All, 2, 2]], #&], dMutantGrowthdRetention = "NA",

(* Find the mutant's growth rate *)
popSizes = Map[{#[[1]], Table[NIntegrate[stateVar/. #[[2, 1]], {t, 3000, 4000}(*, Evaluate[FilterRules[{opts}, Options[NIntegrate]]]*)]/(4000 - 3000), {stateVar, {h[t], p[t], h2[t], p2[t]}}]}&, mutantDynamics];
derivatives = Map[(Total[#[[2, 3;;4]]] / Total[#[[2]]] - 0.01) (* Change in fraction of population compared to resident trait = mutant trait *) / (#[[1]] - retention) (* change in trait *)&, popSizes];

dMutantGrowthdRetention = Mean[derivatives] ] ];

(* Derivative w.r.t. replication *)
Module[{mutants, mutantDynamics, popSizes, derivatives},

mutants = Select[replication + {-1, 1} * 0.0001, 0 <= # <= 1&];

mutantDynamics = Map[{#, mutantAndResidentEcoDynamics[retention, replication, retention, #, surfaceLight, initialPrey, 0.99 * initialHeterotroph, 0.99 * initialPhototroph,  0.01 * initialHeterotroph,  0.01 * initialPhototroph, opts]}&,
	mutants];
	
(* Did introduction of any of the mutants cause an extinction? If so, the derivative of invasion fitness is maybe unclear. Set it to NA *)
If[AnyTrue[mutantDynamics[[All, 2, 2]], #&], dMutantGrowthdReplication = "NA",

(* Find the mutant's growth rate *)
popSizes = Map[{#[[1]], Table[NIntegrate[stateVar/. #[[2, 1]], {t, 3000, 4000}(*, Evaluate[FilterRules[{opts}, Options[NIntegrate]]]*)]/(4000 - 3000), {stateVar, {h[t], p[t], h2[t], p2[t]}}]}&, mutantDynamics];
derivatives = Map[(Total[#[[2, 3;;4]]] / Total[#[[2]]] - 0.01) (* Change in fraction of population compared to resident trait = mutant trait *) / (#[[1]] - replication) (* change in trait *)&, popSizes];

dMutantGrowthdReplication = Mean[derivatives] ] ];

Return[{dMutantGrowthdRetention, dMutantGrowthdReplication}]]


(* ::Section:: *)
(*Evolution*)


lightLevels = {25, 35, 50};


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEAR = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEAR = Flatten[ecoDynamicsCEAR, {{1}, {2,3}, {4}}];


(* ::Text:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


extinctionPointsCEAR = Map[Select[#[[4]] == True&], ecoDynamicsCEAR];
(* Also find the points where the mixotroph population does not go extinct *)
noExtinctionPointsCEAR = Map[Select[#[[4]] == False&], ecoDynamicsCEAR];


(* ::Text:: *)
(*For the points where the mixotroph goes extinct, check whether the prey also go extinct*)


Map[{#[[1]], #[[2]], #[[3]], residentEcoDynamics[#[[2]], #[[3]], #[[1]], 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]}&, extinctionPointsCEAR, {2}];
Map[w[#[[4, 1, 1, 2, 1, 1, 2]]] /. #[[4, 1]]&, %, {2}]
Map[Plot[w[t] /. #[[4, 1]], {t, 0, #[[4, 1, 1, 2, 1, 1, 2]]}, PlotRange -> Full, AxesLabel -> {"Time", "Prey"}]&, %%, {2}]


(* ::Text:: *)
(*Plot the points where the population goes extinct (green dots, overlaid on colors showing mixotroph population sizes)*)


Table[Show[{
(* Plot of the mixotroph population sizes *)
ListDensityPlot[Map[{#[[2]], #[[3]], #[[6]] + #[[7]]}&, ecoDynamicsCEAR[[i]]], 
FrameLabel -> {"Retention", "Replication"}, PlotRangePadding -> 0.02,
PlotLabel -> "(" <> ToUpperCase[Alphabet[]][[i]] <> ") Mixotroph population (cells/mL)\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", 
ImageSize -> 400, PlotLegends -> Automatic, PlotRange -> All, ClippingStyle -> Black, ColorFunction -> "ValentineTones"],

(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsCEAR[[i, All, {2, 3}]], PlotStyle -> {Green, PointSize[Large]}]}],
{i, Length[ecoDynamicsCEAR]}]


(* ::Text:: *)
(*Find the derivative of invasion fitness at the points where the resident is not extinct*)


derivsForPlotsCEAR = 
Map[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEAR, {2}];


(* ::Text:: *)
(*Quick vector plot to get an idea of the magnitude/direction of derivative of invasion fitness. Colors indicate magnitude; arrows indicate direction*)


Table[Show[{
(* Plot invasion fitness derivatives *)
ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsCEAR[[i]]], 
	VectorPoints -> noExtinctionPointsCEAR[[i, All, 2;;3]],
	ColorFunction -> "Aquamarine", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsCEAR[[i, All, {2, 3}]], PlotStyle -> {Pink, PointSize[Large]}]}],
{i, Length[derivsForPlotsCEAR]}]


(* ::Text:: *)
(*Yeah, that's really hard to see arrow direction. Let's make all the arrows the same size. Colors still indicate magnitude of vectors. White squares now indicate points where the mixotroph goes extinct.*)
(*Also I did a pastel color scheme this time because I was having trouble telling the other one apart.*)


Table[
Show[ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsCEAR[[i]]], 
	VectorPoints -> noExtinctionPointsCEAR[[i, All, 2;;3]],
	ColorFunction -> "Pastel", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium, VectorScale -> {Automatic, Automatic, None}], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCEAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCEAR]}]


(* ::Text:: *)
(*Alright, now we'll make some stream plots.*)


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Pastel"(*"Aquamarine"*)][Sqrt[(p/maxPrey)](**7/10 + 3/10*)]


Table[
Show[ListStreamDensityPlot[Map[{#[[1]], {#[[2]], #[[3,1]]}}&, derivsForPlotsCEAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full},
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False, 
	StreamStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCEAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCEAR]}]


extinctionPointsCEAR


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCE = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCE = Flatten[ecoDynamicsCE, {{1}, {2,3}, {4}}];


(* ::Text:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


extinctionPointsCE = Map[Select[#[[4]] == True&], ecoDynamicsCE];
(* Also find the points where the mixotroph population does not go extinct *)
noExtinctionPointsCE = Map[Select[#[[4]] == False&], ecoDynamicsCE];


(* ::Text:: *)
(*Plot the points where the population goes extinct (green dots, overlaid on colors showing mixotroph population sizes)*)


Table[Show[{
(* Plot of the mixotroph population sizes *)
ListDensityPlot[Map[{#[[2]], #[[3]], #[[6]] + #[[7]]}&, ecoDynamicsCE[[i]]], 
FrameLabel -> {"Retention", "Replication"}, PlotRangePadding -> 0.02,
PlotLabel -> "(" <> ToUpperCase[Alphabet[]][[i]] <> ") Mixotroph population (cells/mL)\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", 
ImageSize -> 400, PlotLegends -> Automatic, PlotRange -> All, ClippingStyle -> Black, ColorFunction -> "ValentineTones"],

(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsCE[[i, All, {2, 3}]], PlotStyle -> {Green, PointSize[Large]}]}],
{i, Length[ecoDynamicsCE]}]


(* ::Text:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCE = 
Map[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCE, {2}];


(* ::Text:: *)
(*Quick vector plot to get an idea of the magnitude/direction of derivative of invasion fitness. Colors indicate magnitude; arrows indicate direction*)


Table[Show[{
(* Plot invasion fitness derivatives *)
ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsCE[[i]]], 
	VectorPoints -> noExtinctionPointsCE[[i, All, 2;;3]],
	ColorFunction -> "Aquamarine", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsCE[[i, All, {2, 3}]], PlotStyle -> {Pink, PointSize[Large]}]}],
{i, Length[derivsForPlotsCE]}]


(* ::Text:: *)
(*Yeah, that's really hard to see arrow direction. Let's make all the arrows the same size. Colors still indicate magnitude of vectors. White squares now indicate points where the mixotroph goes extinct.*)
(*Also I did a pastel color scheme this time because I was having trouble telling the other one apart.*)


Table[
Show[ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsCE[[i]]], 
	VectorPoints -> noExtinctionPointsCE[[i, All, 2;;3]],
	ColorFunction -> "Pastel", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium, VectorScale -> {Automatic, Automatic, None}], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCE]}]


(* ::Text:: *)
(*Alright, now we'll make some stream plots.*)


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Aquamarine"][Sqrt[p/maxPrey]]


Table[
Show[ListStreamDensityPlot[Map[{#[[1]], {#[[2]], #[[3,1]]}}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 0,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False, 
	StreamStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCE]}]


(* ::Subsection:: *)
(*Put the plots together*)


(* ::Text:: *)
(*Color scheme (prey density)*)


Histogram[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]]]


Histogram[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]]


maxPrey = Max[{Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]],
	Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]}];
preyDensityColors[p_] := ColorData["BeachColors"][Log[1 + p^2]/Log[1 + maxPrey^2]]


preyDensityColors[p_] := If[p <= 16.6, ColorData["BeachColors"][(9/10) * Log[1 + p]/Log[1 + 16.6]], ColorData["BeachColors"][9/10 + p/(maxPrey * 10)]]


BarLegend[{preyDensityColors[#]&, {0, maxPrey}}]


(* ::Text:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}];
(*SetOptions[SwatchLegend, LabelStyle \[Rule] {FontSize -> 12, FontFamily -> "Arial"}];
SetOptions[PointLegend, LabelStyle \[Rule] {FontSize -> 12, FontFamily -> "Arial"}];*)


streamPlotsCEAR = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlotsCEAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Attack tradeoff\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCEAR[[i]]], StreamStyle -> Black(*, StreamPoints \[Rule] noExtinctionPointsCEAR[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsCEAR[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsCEAR]}]

streamPlotsCE = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1,2]], #[[3, 1]]}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i + Length[lightLevels]]]] <> ") No tradeoff\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCE[[i]]], StreamStyle -> Black(*, StreamPoints \[Rule] noExtinctionPointsCEAR[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsCE]}]


Grid[{streamPlotsCEAR, streamPlotsCE, 
{Row[{PointLegend[{Black}, {"Direction of selection"}, LegendMarkers -> Graphics[{Arrowheads[0.5], Arrow[{{0, 0}, {1, 0}}]}], LegendMarkerSize -> 75],
BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Row", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {450, 25}],
SwatchLegend[{Black}, {"Mixotroph extinction"}, LegendMarkerSize -> {20, 20}]}],
SpanFromLeft, SpanFromLeft}}, Spacings -> {2, 3}]
Export["stream_plots.pdf", %]


(* ::Subsection:: *)
(*Conversion trade-off only (old)*)


ecoDynamics = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 1, 1, 1, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6]],
{light, {25, 50}}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamics = Flatten[ecoDynamics, {{1}, {2,3}, {4}}];


Table[ListDensityPlot[ecoDynamics[[i, All, {2, 3, 5}]], PlotLegends -> Automatic, FrameLabel -> {"Retention", "Replication"}], {i, Length[ecoDynamics]}]


extinctionPoints = Map[Select[#[[4]] == True&], ecoDynamics]


derivsForPlots = 
Map[{{#[[2]], #[[3]]}, {invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True],
#[[5]]}}&, noExtinctionPoints, {2}];


Table[ListStreamPlot[Map[{#[[1]], #[[2, 1]]}&, derivsForPlots[[i]]]], {i, Length[derivsForPlots]}]


Table[
Show[ListDensityPlot[ecoDynamics[[i, All, {2, 3, 5}]], PlotLabel -> "Light = " <> ToString[{15, 20, 25}[[i]]] <> " umol photons/m^2/s",
FrameLabel -> {"Retention", "Replication"}],
ListPlot[extinctionPoints[[i, All, {2, 3}]], PlotStyle -> Green]],
{i, 3}]
(*Export["mmexport1.pdf", %]*)


ecoDynamics[[1, 12]]


noExtinctionPoints = Map[Select[#[[4]] == False&], ecoDynamics];


invasionDerivatives[1/10, 0, 15, 4.565001044002136`,0.49142161562573927`,0.04807157546521535, conversionTradeOff -> True, PrecisionGoal -> 10]


derivsForPlots = 
Map[{{#[[2]], #[[3]]}, {invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6],
#[[5]]}}&, ecoDynamics, {2}];


Min[derivsForPlots[[All, All, 2, 2]]]


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamics][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Aquamarine"][(p/maxPrey)*7/10]


maxPrey


Table[ListStreamDensityPlot[derivsForPlots[[i]],
FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"}, ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
StreamStyle -> Black, PlotRange -> {{0, 1}, {0, 1}},
PlotLabel -> {"(A)", "(B)", "(B)"}[[i]] <> "No trade-off \n surface light = " <> ToString[{25, 35, 50}[[i]]] <> " \[Mu]mol photons/m^2/s", InterpolationOrder -> 1, ImageSize -> 300],
{i, 2}]


Table[ListStreamDensityPlot[acDerivs[[i]],
FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"}, ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
StreamStyle -> Black, PlotRange -> {{0, 1}, {0, 1}},
PlotLabel -> {"(C)", "(C)", "(D)"}[[i]] <> "Attack trade-off \n surface light = " <> ToString[{25, 35, 50}[[i]]] <> " \[Mu]mol photons/m^2/s", InterpolationOrder -> 1, ImageSize -> 300],
{i, {1,3}}]


BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLabel -> "Prey cells/mL", LegendLayout -> "Row"]


?Rasterize


?GraphicsGrid


Rasterize[Grid[{%122, %121, {%116, SpanFromLeft}}, Alignment -> Center], ImageResolution -> 300]
Export["stream_plots.pdf", %]


?BarLegend


(* ::Subsection:: *)
(*Putting the plots together*)


(* ::Subsection:: *)
(*Attack trade-off & CE trade-off*)


ecoDynamicsAT = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 1/10, 1, 2, conversionTradeOff -> True, attackTradeOff -> True]],
{light, {15, 20, 25}}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsAT = Flatten[ecoDynamicsAT, {{1}, {2,3}, {4}}];


extinctionPointsAT = Map[Select[#[[4]] == True&], ecoDynamicsAT]


averageEcoState[1, 1, 20, 0.0000001, 4.6, 1.5, conversionTradeOff -> True, attackTradeOff -> True]


GraphicsColumn[Table[ListDensityPlot[ecoDynamicsAT[[2, 1;;120, {2, 3, i}]], PlotLegends -> Automatic, PlotRange -> All, PlotLabel -> {"Prey", "Heterotroph", "Phototroph"}[[i - 4]]], {i, {5, 6, 7}}], ImageSize -> 300]
Export["Attacktradeoff_ecoEquil.pdf", %]


Table[ecoDynamicsAT[[2, 11*i + j + 1, {2, 3, 5, 6, 7}]], {i, 0, 10}, {j, 0, 10}] // MatrixForm


residentEcoDynamics[1, 1, 15, 0, 0, 3, conversionTradeOff -> True, attackTradeOff -> True]


Plot[p[t]/.%48[[1]], {t, 0, 4000}]


Plot[p[t]/.%114[[1]], {t, 0, 4000}, PlotRange -> All]


noExtinctionPoints = Map[Select[(#[[4]] == False) && (Chop[#[[6]] + #[[7]]] > 0)&], ecoDynamicsAT];


noExtinctionPoints[[2]]


derivsForPlotsAT = 
Map[{{#[[2]], #[[3]]}, {invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], conversionTradeOff -> True, attackTradeOff -> True, PrecisionGoal -> 10, AccuracyGoal -> 10],
#[[5]]}}&, noExtinctionPoints, {2}];


Min[derivsForPlotsAT[[All, All, 2, 2]]]


derivsForPlotsAT[[2]]


preyDensityColors[p_] := ColorData["Aquamarine"][p/(Max[derivsForPlotsAT[[All, All, 2, 2]]]*10/7)]


Table[ListStreamDensityPlot[derivsForPlotsAT[[i]],
FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"}, ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
StreamStyle -> Black, PlotRange -> {{0, 1}, {0, 1}},
PlotLabel -> "AT & CE:\nSelection stream plot, plotted over prey density \n surface light = " <> ToString[{15, 20, 25}[[i]]] <> " \[Mu]mol photons/m^2/s"],
{i, 3}]


Table[ListStreamDensityPlot[derivsForPlotsAT[[i]],
FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"}, ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
StreamStyle -> Black, PlotRange -> {{0, 1}, {0, 1}},
PlotLabel -> "AT & CE:\nSelection stream plot, plotted over prey density \n surface light = " <> ToString[{15, 20, 25}[[i]]] <> " \[Mu]mol photons/m^2/s"],
{i, 3}]


BarLegend[{preyDensityColors[#]&, {0, Max[derivsForPlotsAT[[All, All, 2, 2]]]}}, LegendLayout -> "Row"]


derivsForPlots2 = 
Map[{{#[[2]], #[[3]]}, {invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], conversionTradeOff -> True],
#[[5]]}}&, Map[Select[#[[4]] == True&], ecoDynamicsAT], {2}];
