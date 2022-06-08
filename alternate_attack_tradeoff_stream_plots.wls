#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
(*Parameters*)


(* ::Text:: *)
(*Note on the parameter c_HP in the paper: In the code we have replaced c_HP with its reciprocal, gamma_C. This doesn't come up anywhere for the user, but just be aware that setting gamma_C = 0.7 in the functions below is equivalent to setting c_HP to 10/7. *)


(* ::Text:: *)
(*These parameters produce evolution to different ESSs at different light levels (25 vs. 1000) with conversion trade +/- attack trade off (same ESSs with either)*)
(*Removing conversion trade-off means we lose this, partly I think because of changes in ability of the mixotroph to survive/maintain decent population sizes at high retention rates*)


pW = 3;
gammaP = 2/3;
sigmaP = 0;
kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;
hW = 50;
hP = hP2 = 50;
lW = 0.5;
lH = lH2 = 0.25;
lP = lP2 = 0.25;
dMax = 0.3;
sigmaD = 0; (* No more decay trade-off *)
a = a2 = 0.15;
gammaA = 0; (* We want this to be 0 consistently, because we've decided the phototroph doesn't eat *)
c = 0.1;
sigmaC = 0; (* We don't use this, the equations just require it to be defined *)


(* ::Section::Closed:: *)
(*Model equations (alternate attack trade-off)*)


(* ::Text:: *)
(*Here we model the attack rate tradeoff as a = a_max * ( 1 - (phi + eta)/2 ).*)


(* ::Text:: *)
(*Prey*)


wDot[w_, hW_, lW_, a_, pW_, kappa_, iIn_, gammaA_, eta_, phi_, eta2_, phi2_, sigmaA_, p_, h_, a2_, h2_, p2_]:= 
(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-(a(1-sigmaA (phi + eta)/2)) w h-(gammaA a(1-sigmaA (phi + eta)/2)) w p-(a2(1-sigmaA (phi2 + eta2)/2)) w h2 - (gammaA a2(1-sigmaA (phi2 + eta2)/2)) w p2 


(* ::Text:: *)
(*Resident heterotroph*)


hDot[eta_, c_, gammaC_, kappa_, gammaP_, pW_, sigmaP_, p_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, lH_, dMax_, sigmaD_]:=
(1-eta) (c/(gammaC c))((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi (a(1-sigmaA (phi + eta)/2)) w h + (1-phi)c (a(1-sigmaA (phi + eta)/2)) w h+ (1-eta) c (gammaA a(1-sigmaA (phi + eta)/2)) w p -lH h + (dMax (1-sigmaD eta))p


(* ::Text:: *)
(*Resident phototroph*)


pDot[eta_, gammaP_, sigmaP_, pW_, p_, kappa_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, gammaC_, sigmaC_, c_, dMax_, sigmaD_, lP_] :=
eta ((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi (a(1-sigmaA (phi + eta)/2)) w h +eta (gammaA a(1-sigmaA (phi + eta)/2)) (gammaC c) p w - (dMax (1-sigmaD eta)) p - lP p


(* ::Text:: *)
(*Mutant heterotroph*)


h2Dot[eta2_, c_, gammaC_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, lH2_, dMax_, sigmaD_] :=
(1-eta2) (c/(gammaC c))((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 (a2(1-sigmaA (phi2 + eta2)/2)) w h2 + (1-phi2) c (a2(1-sigmaA (phi2 + eta2)/2)) w h2 +(1-eta2) c (gammaA a2(1-sigmaA (phi2 + eta2)/2)) w p2 - lH2 h2+ (dMax(1-sigmaD eta2)) p2


(* ::Text:: *)
(*Mutant phototroph*)


p2Dot[eta2_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, c_, gammaC_, sigmaC_, dMax_, sigmaD_, lP2_] :=
eta2 ((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 (a2(1-sigmaA (phi2 + eta2)/2)) w h2+eta2 (gammaA a2(1-sigmaA (phi2 + eta2)/2)) (gammaC c) p2 w-(dMax(1-sigmaD eta2)) p2-lP2 p2


(* ::Text:: *)
(*Light absorbance*)


kappa[kW_, kH_, kP_, kH2_, kP2_, w_, h_, p_, h2_, p2_]:=
kW * w + kH * h + kP * p + kH2 * h2 + kP2 * p2 


(* ::Section::Closed:: *)
(*Functions for finding eco-evolutionary dynamics*)


(* ::Subsection::Closed:: *)
(*Set default options for analyses*)


(* ::Text:: *)
(*Set the default options. Changes these when calling the functions below, not here!*)


(* ::Text:: *)
(*Default option values:*)


(* ::Item:: *)
(*Attack tradeoff (option attackTradeOff) : off . Retention and replication traits do not affect attack tradeoff .*)


(* ::Item:: *)
(*Conversion efficiency tradeoff (option conversionTradeOff) : off . The conversion between phototrophic and heterotrophic biomass is 1. (We don' t use this, we always use a conversion efficiency of c_HP instead .)*)


(* ::Item:: *)
(*NDSolve default options*)


Options[residentEcoDynamics] = Options[averageEcoState] =
Options[minimumMixotroph] =
Options[mutantAndResidentEcoDynamics] =
Options[invasionDerivatives] =
(* The options *)
Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[NDSolve]];


(* ::Subsection::Closed:: *)
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
	WhenEvent[(h[t] <= 10^-6) && (p[t] <= 10^-6), extinction = True; "StopIntegration"],
	
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
 extinction = ((h[4000] <= 10^-6) && (p[4000] <= 10^-6)) /. ecoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{ecoDynamics, extinction}]]


(* ::Subsection::Closed:: *)
(*Function to find the average ecological state, with no mutants present*)


averageEcoState[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), initialPrey_ (* W[0] *), initialHeterotroph_ (* H[0] *), initialPhototroph_ (* P[0] *), opts:OptionsPattern[]] :=
Module[{ecoDynamics, extinction, averagePrey, averageHeterotroph, averagePhototroph, averageMixotroph},

(* Get the ecological dynamics *)
{ecoDynamics, extinction} = residentEcoDynamics[retention, replication, surfaceLight, initialPrey, initialHeterotroph, initialPhototroph, opts];


(* If the mixotroph went extinct, the mixotroph population sizes are 0. Maybe want to calculate the average prey population size with no mixotrophs, starting from initialPrey? Don't do this now *)
If[extinction,
	averagePrey = "Unknown";
	averageHeterotroph = averagePhototroph = averageMixotroph = 0,

(* If the mixotroph did not go extinct, calculate the average ecological state. Average the value over the last 1000 simulated time steps *)

{averagePrey, averageHeterotroph, averagePhototroph, averageMixotroph} = 
Table[NIntegrate[stateVar/. ecoDynamics, {t, 3000, 4000}, Evaluate[FilterRules[{opts}, Options[NIntegrate]]]]/(4000 - 3000),
	{stateVar, {w[t], h[t], p[t], h[t] + p[t]}}]];
	
(* Return the average population sizes and whether mixotroph went extinct *)
Return[{extinction, averagePrey, averageHeterotroph, averagePhototroph, averageMixotroph}]]


(* ::Subsection::Closed:: *)
(*Function to find the minimum mixotroph population size, with no mutant present*)


minimumMixotroph[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), initialPrey_ (* W[0] *), initialHeterotroph_ (* H[0] *), initialPhototroph_ (* P[0] *), opts:OptionsPattern[]] :=
Module[{ecoDynamics, extinction, minMixotroph},

(* Get the ecological dynamics *)
{ecoDynamics, extinction} = residentEcoDynamics[retention, replication, surfaceLight, initialPrey, initialHeterotroph, initialPhototroph, opts];


(* If the mixotroph went extinct, the mixotroph population sizes are 0. Maybe want to calculate the average prey population size with no mixotrophs, starting from initialPrey? Don't do this now *)
If[extinction,
	minMixotroph = 0,

(* If the mixotroph did not go extinct, find the minimum mixotroph population size reached over the last 1000 simulated time steps *)

	minMixotroph = 
		NMinimize[{(h[t] + p[t])/. ecoDynamics, 3000 <= t <= 4000}, t, Evaluate[FilterRules[{opts}, Options[NMinimize]]]][[1]]]; (* [[1]] gives us the minimum [[2]] would give argMin]*)
		
(* Return the average population sizes and whether mixotroph went extinct *)
Return[minMixotroph]]


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
	WhenEvent[((h[t] <= 10^-6) && (p[t] <= 10^-6)) && ((h2[t] <= 10^-6) && (p2[t] <= 10^-6)), extinction = True; "StopIntegration"],
	
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
 extinction = ((h[4000] <= 10^-6) && (p[4000] <= 10^-6) && (h2[4000] <= 10^-6) && (p2[4000] <= 10^-6)) /. ecoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{ecoDynamics, extinction}]]


(* ::Subsection::Closed:: *)
(*Function to find the derivative of invasion fitness*)


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
popSizes = Map[{#[[1]], Table[NIntegrate[stateVar/. #[[2, 1]], {t, 3000, 4000}]/(4000 - 3000), {stateVar, {h[t] + p[t], h2[t] + p2[t]}}]}&, mutantDynamics];
derivatives = Map[(#[[2, 2]] / Total[#[[2]]] - 0.01) (* Change in fraction of population compared to resident trait = mutant trait *) / (#[[1]] - retention) (* change in trait *)&, popSizes];

dMutantGrowthdRetention = Mean[derivatives] ] ];

(* Derivative w.r.t. replication *)
Module[{mutants, mutantDynamics, popSizes, derivatives},

mutants = Select[replication + {-1, 1} * 0.0001, 0 <= # <= 1&];

mutantDynamics = Map[{#, mutantAndResidentEcoDynamics[retention, replication, retention, #, surfaceLight, initialPrey, 0.99 * initialHeterotroph, 0.99 * initialPhototroph,  0.01 * initialHeterotroph,  0.01 * initialPhototroph, opts]}&,
	mutants];
	
(* Did introduction of any of the mutants cause an extinction? If so, the derivative of invasion fitness is maybe unclear. Set it to NA *)
If[AnyTrue[mutantDynamics[[All, 2, 2]], #&], dMutantGrowthdReplication = "NA",

(* Find the mutant's growth rate *)
popSizes = Map[{#[[1]], Table[NIntegrate[stateVar/. #[[2, 1]], {t, 3000, 4000}(*, Evaluate[FilterRules[{opts}, Options[NIntegrate]]]*)]/(4000 - 3000), 
	{stateVar, {h[t] + p[t], h2[t] + p2[t]}}]}&, mutantDynamics];
derivatives = Map[(#[[2, 2]] / Total[#[[2]]] - 0.01) (* Change in fraction of population compared to resident trait = mutant trait *) / (#[[1]] - replication) (* change in trait *)&, popSizes];

dMutantGrowthdReplication = Mean[derivatives] ] ];

Return[{dMutantGrowthdRetention, dMutantGrowthdReplication}]]


(* ::Section::Closed:: *)
(*Find ecological equilibrium & derivative of invasion fitness (alternate attack trade-off)*)


(* ::Text:: *)
(*Retention and replication additively decrease attack rate*)


lightLevels = {25, 35, 50};


(* ::Text:: *)
(*We'll only do this for the attack tradeoff present case, because obviously if the attack tradeoff isn't present the results don't depend on how we model the attack tradeoff.*)


(* ::Subsection::Closed:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsAdditiveAR = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsAdditiveAR = Flatten[ecoDynamicsAdditiveAR, {{1}, {2,3}, {4}}];


(* ::Subsection::Closed:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


minPopSizesAdditiveAR = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesAdditiveAR = Flatten[minPopSizesAdditiveAR, {{1}, {2, 3}}];


extinctionPointsAdditiveAR = Map[Select[#[[4]] <= 10^-6&], minPopSizesAdditiveAR];


(* ::Text:: *)
(*Where does the population not go extinct?*)


noExtinctionPointsAdditiveAR = Pick[ecoDynamicsAdditiveAR, Map[#[[4]] > 10^-6&, minPopSizesAdditiveAR, {2}]];


(* ::Subsection::Closed:: *)
(*Plot the points where the population goes extinct (green dots, overlaid on colors showing mixotroph population sizes)*)


Table[Show[{
(* Plot of the mixotroph population sizes *)
ListDensityPlot[Map[{#[[2]], #[[3]], #[[6]] + #[[7]]}&, ecoDynamicsAdditiveAR[[i]]], 
FrameLabel -> {"Retention", "Replication"}, PlotRangePadding -> 0.02,
PlotLabel -> "(" <> ToUpperCase[Alphabet[]][[i]] <> ") Mixotroph population (cells/mL)\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", 
ImageSize -> 400, PlotLegends -> Automatic, PlotRange -> All, ClippingStyle -> Black, ColorFunction -> "ValentineTones"],

(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsAdditiveAR[[i, All, {2, 3}]], PlotStyle -> {Green, PointSize[Large]}]}],
{i, Length[ecoDynamicsAdditiveAR]}]


(* ::Subsection::Closed:: *)
(*Find the derivative of invasion fitness at the points where the resident is not extinct*)


derivsForPlotsAdditiveAR = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsAdditiveAR, {2}];


(* ::Text:: *)
(*Quick vector plot to get an idea of the magnitude/direction of derivative of invasion fitness. Colors indicate magnitude; arrows indicate direction*)


Table[
Show[ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsAdditiveAR[[i]]], 
	VectorPoints -> noExtinctionPointsAdditiveAR[[i, All, 2;;3]],
	ColorFunction -> "Pastel", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium, VectorScale -> {Automatic, Automatic, None}], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsAdditiveAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsAdditiveAR]}]


(* ::Subsection::Closed:: *)
(*Some stream plots*)


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsAdditiveAR][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Pastel"(*"Aquamarine"*)][Sqrt[(p/maxPrey)](**7/10 + 3/10*)]


Table[
Show[ListStreamDensityPlot[Map[{#[[1]], {#[[2]], #[[3,1]]}}&, derivsForPlotsAdditiveAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full},
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False, 
	StreamStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsAdditiveAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsAdditiveAR]}]


(* ::Section::Closed:: *)
(*Plot the original trade-off for comparison*)


(* ::Subsection::Closed:: *)
(*Load the model for the original attack trade-off (retention and replication multiplicatively decrease attack rate)*)


Get[NotebookDirectory[] <> "model_and_functions.wls"]


(* ::Subsection::Closed:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsMultiplicativeAR = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsMultiplicativeAR = Flatten[ecoDynamicsMultiplicativeAR, {{1}, {2,3}, {4}}];


(* ::Subsection::Closed:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


minPopSizesMultiplicativeAR = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesMultiplicativeAR = Flatten[minPopSizesMultiplicativeAR, {{1}, {2, 3}}];


extinctionPointsMultiplicativeAR = Map[Select[#[[4]] <= 10^-6&], minPopSizesMultiplicativeAR];


(* ::Text:: *)
(*Where does the population not go extinct?*)


noExtinctionPointsMultiplicativeAR = Pick[ecoDynamicsMultiplicativeAR, Map[#[[4]] > 10^-6&, minPopSizesMultiplicativeAR, {2}]];


(* ::Subsection::Closed:: *)
(*Plot the points where the population goes extinct (green dots, overlaid on colors showing mixotroph population sizes)*)


Table[Show[{
(* Plot of the mixotroph population sizes *)
ListDensityPlot[Map[{#[[2]], #[[3]], #[[6]] + #[[7]]}&, ecoDynamicsMultiplicativeAR[[i]]], 
FrameLabel -> {"Retention", "Replication"}, PlotRangePadding -> 0.02,
PlotLabel -> "(" <> ToUpperCase[Alphabet[]][[i]] <> ") Mixotroph population (cells/mL)\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", 
ImageSize -> 400, PlotLegends -> Automatic, PlotRange -> All, ClippingStyle -> Black, ColorFunction -> "ValentineTones"],

(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsMultiplicativeAR[[i, All, {2, 3}]], PlotStyle -> {Green, PointSize[Large]}]}],
{i, Length[ecoDynamicsMultiplicativeAR]}]


(* ::Subsection::Closed:: *)
(*Find the derivative of invasion fitness at the points where the resident is not extinct*)


derivsForPlotsMultiplicativeAR = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsMultiplicativeAR, {2}];


(* ::Text:: *)
(*Quick vector plot to get an idea of the magnitude/direction of derivative of invasion fitness. Colors indicate magnitude; arrows indicate direction*)


Table[
Show[ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsMultiplicativeAR[[i]]], 
	VectorPoints -> noExtinctionPointsMultiplicativeAR[[i, All, 2;;3]],
	ColorFunction -> "Pastel", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium, VectorScale -> {Automatic, Automatic, None}], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsMultiplicativeAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsMultiplicativeAR]}]


(* ::Subsection:: *)
(*Make some stream plots*)


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsMultiplicativeAR][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Pastel"(*"Aquamarine"*)][Sqrt[(p/maxPrey)](**7/10 + 3/10*)]


Table[
Show[ListStreamDensityPlot[Map[{#[[1]], {#[[2]], #[[3,1]]}}&, derivsForPlotsMultiplicativeAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full},
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False, 
	StreamStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsMultiplicativeAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsMultiplicativeAR]}]


(* ::Section:: *)
(*Plot the no attack trade-off case for comparison*)


(* ::Subsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCE = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCE = Flatten[ecoDynamicsCE, {{1}, {2,3}, {4}}];


(* ::Subsection:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


minPopSizesCE = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCE = Flatten[minPopSizesCE, {{1}, {2, 3}}];


extinctionPointsCE = Map[Select[#[[4]] <= 10^-6&], minPopSizesCE];


(* ::Text:: *)
(*Where does the population not go extinct?*)


noExtinctionPointsCE = Pick[ecoDynamicsCE, Map[#[[4]] > 10^-6&, minPopSizesCE, {2}]];


(* ::Subsection:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCE = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCE, {2}];


(* ::Text:: *)
(*Quick vector plot to get an idea of the magnitude/direction of derivative of invasion fitness. Colors indicate magnitude; arrows indicate direction*)


Table[
Show[ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsCE[[i]]], 
	VectorPoints -> noExtinctionPointsCE[[i, All, 2;;3]],
	ColorFunction -> "Pastel", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium, VectorScale -> {Automatic, Automatic, None}], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCE]}]


(* ::Subsection:: *)
(*Make some stream plots*)


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Pastel"][Sqrt[p/maxPrey]]


Table[
Show[ListStreamDensityPlot[Map[{#[[1]], {#[[2]], #[[3,1]]}}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 0,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False, 
	StreamStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCE]}]


(* ::Section:: *)
(*Final plots*)


(* ::Subsection:: *)
(*Set up plot properties*)


(* ::Subsubsection:: *)
(*Color scheme (prey density)*)


ClearAll[preyDensityColors]


maxPrey = Max[{noExtinctionPointsAdditiveAR[[All, All, 5]], noExtinctionPointsMultiplicativeAR[[All, All, 5]], noExtinctionPointsCE[[All, All, 5]]}];


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


(* ::Subsubsection:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black}];
(*SetOptions[SwatchLegend, LabelStyle \[Rule] {FontSize -> 12, FontFamily -> "Arial"}];
SetOptions[PointLegend, LabelStyle \[Rule] {FontSize -> 12, FontFamily -> "Arial"}];*)


(* ::Subsection::Closed:: *)
(*Make individual plots*)


streamPlotsAdditiveAR = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1,2]], #[[3, 1]]}&, derivsForPlotsAdditiveAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> Pane[Style["(" <> Capitalize[Alphabet[][[i + 2 * Length[lightLevels]]]] <> ")", {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black}], ImageSize -> 300], ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsAdditiveAR[[i]]], StreamStyle -> Black(*, StreamPoints \[Rule] noExtinctionPointsAdditiveAR[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsAdditiveAR[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsAdditiveAR]}]


streamPlotsMultiplicativeAR = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1,2]], #[[3, 1]]}&, derivsForPlotsMultiplicativeAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> Pane[Style["(" <> Capitalize[Alphabet[][[i + Length[lightLevels]]]] <> ")", {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black}], ImageSize -> 300], ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsMultiplicativeAR[[i]]], StreamStyle -> Black(*, StreamPoints \[Rule] noExtinctionPointsAdditiveAR[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsMultiplicativeAR[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsMultiplicativeAR]}]


streamPlotsCE = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel ->Pane[Style["(" <> Capitalize[Alphabet[][[i]]] <> ")", {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black}], ImageSize -> 300], ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCE[[i]]], StreamStyle -> Black(*, StreamPoints \[Rule] noExtinctionPointsCE[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsCE]}]


(* ::Subsection:: *)
(*Combine plots with legend and export*)


(* ::Text:: *)
(*In order to get the plot labels to not have extra line breaks, first go into print settings and change the paper size to legal and the orientation to landscape*)


Column[{Row[{Grid[

	(* Row labels *)
	Join[Map[Rotate[Style[#, {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black, TextAlignment -> Center}], 90Degree]&, 
		{{"", ""}, {"No tradeoff", "Replication probability (\[Eta])"}, {"Multiplicative Attack tradeoff", "Replication probability (\[Eta])"}, {"Additive Attack tradeoff", "Replication probability (\[Eta])"}, {"", ""}}, {2}],
		
	(* Column labels *)
	Join[Map[Style[#, {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black, TextAlignment -> Center}]&, 
		{Table[ToString[light] <> " \[Mu]mol photons/" <> ToString[Superscript["m", "2"], StandardForm] <> "/s", {light, lightLevels}]},
		{2}],
		
	(* Some extremely fancy rasterization of the plots without rasterizing their labels (removes the labels, raterizes, and re-adds the labels *)
	Map[
		Labeled[Rasterize[Show[#, FrameLabel -> None, PlotLabel -> None](*, RasterSize \[Rule] 500*), ImageSize -> Medium], 
		(*Append[(FrameLabel /. Options[#])[[All, 1]], *)Style[(PlotLabel /. Options[#]), TextAlignment -> Center](*]*), Top(*{Left, Bottom, Top}*),
		LabelStyle -> Append[(LabelStyle /. Options[#]), (TextAlignment -> Center)], RotateLabel -> True]&,
	{streamPlotsCE, streamPlotsMultiplicativeAR, streamPlotsAdditiveAR}, {2}],

	(* Finish up adding the row and column labels *)	
	{Table[Style["Retention probability (\[Phi])", {FontSize -> 28, FontFamily -> "Arial", FontColor -> Black, TextAlignment -> Center}], {i, 3}]}], 2],
	
	Spacings -> {1, 1}],
 
	(* Legend *)
		Column[{
		BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Column", LegendLabel -> "Prey density\n(cells/mL)", LegendMarkerSize -> {25, 450}],
		PointLegend[{Black, Black}, {"Direction of\nselection", "Consumer\nextinction"}, 
			LegendMarkers -> {Graphics[{Arrowheads[0.5], Arrow[{{0, 0}, {1, 0}}]}], Graphics[Rectangle[{0, 0}]]},
			LegendMarkerSize -> {{75, 75}, {20, 20}}]}, 
	Frame -> True, Alignment -> {Center, Center}]}, 
Spacer[50]],
Graphics[{Thickness[0.0025], Arrowheads[0.025], Arrow[{{0, 0}, {100, 0}}], Inset[Style["Increasing light", 28], Background -> White]}, ImageSize -> {1200, 50}]}]

Export["Alternate_attack_tradeoff.pdf", %]


(* ::Section::Closed:: *)
(*Plot the attack rate at each point for both trade-offs*)


Row[{ContourPlot[1 - phi * eta, {phi, 0, 1}, {eta, 0, 1}, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"}, PlotLabel -> "(A) Multiplicative trade-off", 
	ColorFunction -> "GrayTones", ColorFunctionScaling -> None, LabelStyle -> {Black, "Arial", 18}, ImageSize -> 400, Contours -> Range[0, 1, 1/10]],
	
ContourPlot[1 - (phi + eta)/2, {phi, 0, 1}, {eta, 0, 1}, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"}, PlotLabel -> "(B) Additive trade-off", 
	ColorFunction -> "GrayTones", ColorFunctionScaling -> None, LabelStyle -> {Black, "Arial", 18}, ImageSize -> 400, Contours -> Range[0, 1, 1/10]],
	
BarLegend["GrayTones", Range[0, 1, 1/10], LegendLabel -> "Fraction of maximum\nattack rate", LabelStyle -> {Black, "Arial", 18}, LegendMarkerSize -> 300]}, Spacer[20]]

Export["Attack_rate_vs_tradeoff.pdf", %]


(* ::Section:: *)
(*Grab the min/maximum invasion fitness for something I need later*)


Sort[Map[1 - Exp[-2 * Norm[#]]&, Cases[Flatten[derivsForPlotsMultiplicativeAR[[All, All, 2]], 1], Except[{_String, _String}]]  * 1/10000 (* Undo the derivative *) ]]
Histogram[%, PlotRange -> {Full, Full}]
{Min[%%], Max[%%]}


Sort[Map[1 - Exp[-2 * Norm[#]]&, Cases[Flatten[derivsForPlotsCE[[All, All, 2]], 1], Except[{_String, _String}]]  * 1/10000 (* Undo the derivative *) ]]
Histogram[%, PlotRange -> {Full, Full}]
{Min[%%], Max[%%]}



