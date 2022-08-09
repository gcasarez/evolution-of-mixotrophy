#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
(*Model equations*)


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
(*Set default options for analyses*)


(* ::Text:: *)
(*Set the default options. Changes these when calling the functions below, not here!*)


(* ::Text:: *)
(*Default option values:*)


(* ::Item:: *)
(*Attack tradeoff (option attackTradeOff) : false . Retention and replication traits do not affect attack rate . If true, the traits decrease the attack rate.*)


(* ::Item:: *)
(*Conversion efficiency tradeoff (option conversionTradeOff) : false . The conversion between phototrophic and heterotrophic biomass is 1. If true, the conversion is  c_HP instead. We always set this option to true when calling the function, but preserve the default as false for backwards-compatibility with earlier versions of the model.*)


(* ::Item:: *)
(*NDSolve default options*)


Options[residentEcoDynamics] = 
Options[averageEcoState] =
Options[minimumMixotroph] =
Options[mutantAndResidentEcoDynamics] =
Options[invasionDerivatives] =
Options[rareMutantEcoDynamics] =
Options[simulateEvolution] =
(* The options *)
Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[NDSolve]];


(* ::Item:: *)
(*The function to find the minimum mixotroph also takes the options for NMinimize*)


Options[minimumMixotroph] = Join[Options[minimumMixotroph], Options[NMinimize]];


(* ::Text:: *)
(*For the functions that directly simulate ecological dynamics (residentEcoDynamics and mutantAndResidentEcoDynamics), there is an extra option:*)


(* ::Item:: *)
(*Maximum time to simulate ecological dynamics (option maxTime): 4000 days.*)


(* ::Subitem:: *)
(*Note that all but one function that calls either of residentEcoDynamics or mutantAndResidentEcoDynamics uses a maxTime of 4000. *)


(* ::Subitem:: *)
(*The function simulateEvolution calls residentEcoDynamics with a max time of 8000 days, due to needing to allow the resident to equilibrate before introducing the mutant.*)


Options[residentEcoDynamics] = Join[Options[residentEcoDynamics], {maxTime -> 4000}];
Options[mutantAndResidentEcoDynamics] = Join[Options[mutantAndResidentEcoDynamics], {maxTime -> 4000}];


(* ::Text:: *)
(*Options for function to simulate evolution*)


Options[simulateEvolution] = Join[{attackTradeOff -> False, conversionTradeOff -> False, maxMutations -> 500, initialRetention -> 0, initialReplication -> 0}, 
	Options[NDSolve], Options[NMinimize], Options[NMaximize], Options[FindRoot]];


(* ::Section:: *)
(*Ecological dynamics*)


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

(* Get the time to simulate ecological dynamics for *)
tmax = OptionValue[maxTime],

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
	
	(* Find the population sizes over tmax days *)
	{w, h, p}, {t, 0, tmax},
	
	(* Optional input to adjust the differential equation solver *)
	FilterRules[{opts}, Options[NDSolve]]]; 
	
(* Clean up the formatting of ecoDynamics *)
ecoDynamics = Flatten[ecoDynamics, 1];

(* Check that extinction did not occur on the last time step *)
If[extinction == False,
 extinction = ((h[tmax] <= 10^-6) && (p[tmax] <= 10^-6)) /. ecoDynamics];

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

(* Get the time to simulate ecological dynamics for *)
tmax = OptionValue[maxTime],

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
	
	(* Find the population sizes over tmax days *)
	{w, h, p, h2, p2}, {t, 0, tmax},
	
	(* Optional input to adjust the differential equation solver *)
	FilterRules[{opts}, Options[NDSolve]]]; 
	
(* Clean up the formatting of ecoDynamics *)
ecoDynamics = Flatten[ecoDynamics, 1];

(* Check that extinction did not occur on the last time step *)
If[extinction == False,
 extinction = ((h[tmax] <= 10^-6) && (p[tmax] <= 10^-6) && (h2[tmax] <= 10^-6) && (p2[tmax] <= 10^-6)) /. ecoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{ecoDynamics, extinction}]]


(* ::Section::Closed:: *)
(*Invasion fitness*)


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


(* ::Section:: *)
(*Simulate evolution*)


(* ::Subsection::Closed:: *)
(*Function to find rare mutant dynamics*)


(* ::Text:: *)
(*I recommend setting initialMutantHeterotroph and initialMutantPhototroph to the resident's values for best computational accuracy. (This function should work as long these are in the same proportions, or even if they aren't and there is enough burn-in. In practice, using the exact same numbers helps with the precision of calculations.)*)


rareMutantEcoDynamics[preyEcology_, heterotrophEcology_, phototrophEcology_, mutantRetention_ (* phi2 *), mutantReplication_ (* eta2 *), surfaceLight_ (* iIn *), 
initialMutantHeterotroph_, initialMutantPhototroph_, tmin_, tmax_ (* time to end simulation*), opts:OptionsPattern[]] :=

(* Set up a module to hold some temporary parameter values *)
Module[{
(* This variable will hold the solution to the differential equations describing the mutant ecological dynamics *)
mutantEcoDynamics,

(* Match input with the parameters used in equations *)
phi2 = mutantRetention, eta2 = mutantReplication,
iIn = surfaceLight, 

(* Implement trade-offs *)
sigmaA = If[OptionValue[attackTradeOff], 1, 0], (* 1 for attack rate trade-off, 0 for no attack rate trade-off *)
gammaC = If[OptionValue[conversionTradeOff], 0.7, 1], (* 0.7 for conversion efficiency trade-off, 1 for no conversion efficiency trade-off *)

(* Set up a variable that will become true if both the resident and mutant mixotroph populations cannot persist *)
extinction = False},

(* Quick check to make sure we're not starting with the resident extinct *)
If[(heterotrophEcology[0] == 0) && (phototrophEcology[0] == 0),
Print["Resident appears to be extinct?"];
extinction = True];

(* Now onto the equation solving!*)
mutantEcoDynamics = NDSolve[
	{h2'[t] == h2Dot[eta2, c, gammaC, gammaP, sigmaP, pW, p2[t], 
		kappa[kW, kH, kP, 0, 0, preyEcology[t], heterotrophEcology[t], phototrophEcology[t], 0, 0], 
		hP2, iIn, phi2, a2, preyEcology[t], h2[t], gammaA, sigmaA, lH2, dMax, sigmaD],
	p2'[t] == p2Dot[eta2, gammaP, sigmaP, pW, p2[t], 
		kappa[kW, kH, kP, 0, 0, preyEcology[t], heterotrophEcology[t], phototrophEcology[t], 0, 0], 
		hP2, iIn, phi2, a2, preyEcology[t], h2[t], gammaA, sigmaA, c, gammaC, sigmaC, dMax, sigmaD, lP2],
	
	(* Initial conditions for the prey, heterotroph, and phototroph population sizes *)
	h2[tmin] == initialMutantHeterotroph, p2[tmin] == initialMutantPhototroph,
	
	(* If the entire mixotroph population goes to 0, stop solving the differential equations and set the variable "extinction" to true *)
	WhenEvent[((h2[t] <= 10^-6) && (p2[t] <= 10^-6)), extinction = True; "StopIntegration"],
	
	(* If any single population goes below 0, set it to 0. This is just to be safe; technically the equations shouldn't allow this, but sometimes small dips below 0 happen when solving numerically *)
	WhenEvent[h2[t] < 0, h2[t] -> 0], WhenEvent[p2[t] < 0, p2[t] -> 0]},
	
	(* Find the mutant population sizes over time (relative to the initial mutant population) *)
	{h2, p2}, {t, tmin, tmax},
	
	(* Optional input to adjust the differential equation solver *)
	FilterRules[{opts}, Options[NDSolve]]]; 
	
(* Clean up the formatting of ecoDynamics *)
mutantEcoDynamics = Flatten[mutantEcoDynamics, 1];

(* Check that extinction did not occur on the last time step *)
If[extinction == False,
 extinction = ((h2[tmax] <= 0) && (p2[tmax] <= 0)) /. mutantEcoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{mutantEcoDynamics, extinction}]]


(* ::Subsection::Closed:: *)
(*Helper function that estimate the period of an interpolating function*)


(* ::Text:: *)
(*This function will be used later to calculate the period of resident population cycling.*)
(*We need the period to more accurately calculate resident growth rate, so that we can measure mutant population changes over a whole number of cycles.*)


Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];


(* ::Text:: *)
(*Error messages for the function. It requires an interpolating function with*)
(* - a 1-dimensional domain*)
(* - at least two peaks, so the period can be estimated from the distance between the peaks*)


interpolatingFunctionPeriod::DomainError = "Domain is not 1-dimensional.";
interpolatingFunctionPeriod::NotPeriodic = "Peaks detected: `1`. Check that function is periodic and include a larger domain if so.";
interpolatingFunctionPeriod::FewPeaks = "Cannot discard first 3 estimates of period, because only `1` estimates can be obtained.";


(* ::Text:: *)
(*This function takes an interpolating function as input and estimates the period of its cycles based on the median distance between successive peaks.*)
(*If the optional input ignoreEarlyPeaks is true (default value is true), then the the distances between the first 4 peaks are ignored when calculating the median period.*)


interpolatingFunctionPeriod[ifun_, ignoreEarlyPeaks_:True] :=
	Module[{peakLocations, peakSeparation}, 
	
	(* Check that the interpolating function domain is 1-dimensional *)
	If[Dimensions[InterpolatingFunctionDomain[ifun]] != {1, 2},
		Message[interpolatingFunctionPeriod::DomainError];
		Return[]];
		
	(* Find the locations of the peaks of the interpolating function *)
	peakLocations = Pick[Flatten[InterpolatingFunctionGrid[ifun]], PeakDetect[InterpolatingFunctionValuesOnGrid[ifun]], 1];
	
	(* Check that there are at least 2 peaks (so we know we've captured a full period) *)
	If[Length[peakLocations] <= 1,
		Message[interpolatingFunctionPeriod::NotPeriodic, Length[peakLocations]];
		Return[0]];
	
	(* Find the distance between each peak and the next *)
	peakSeparation = Differences[peakLocations];
	
	(* If ignoreEarlyPeaks is true, and we have enough peaks, we'll throw out the first 3 peaks *)
	If[ignoreEarlyPeaks,
		If[Length[peakSeparation] >= 4,
			peakSeparation = peakSeparation[[4;;All]],
			Message[interpolatingFunctionPeriod::FewPeaks, Length[peakSeparation]]
		]
	];
	
	(* Return estimate of period *)
	Return[Median[peakSeparation]]]


(* ::Subsection::Closed:: *)
(*Function to simulate evolution *)


(* ::Text:: *)
(*Takes the surface light as input and simulates evolution from heterotrophy (default) or from another set of initial retention and replication values, if the options initialRetention and initialReplication are given.*)
(*Use attackTradeOff -> True and conversionTradeOff -> True to turn on the attack tradeoff and to set a conversion factor of cHP between heterotrophic and phototrophic biomass. (cHP set elsewhere with the other parameters.)*)


simulateEvolution[surfaceLight_, opts:OptionsPattern[]] :=
Module[{
(* Variables we'll use inside this function *)
residentRetention, residentReplication, (* the resident's traits *)
residentEcology, (* the ecological dynamics with just the resident present *)
extinction, (* will be True if the entire resident population goes extinct *)
minResident, maxResident, (* minimum and maximum resident population size *)
period, (* period of population cycling *)
mutants, (* list of possible mutants that might arise *)
fixationProbabilities, (* list of fixation probabilities of possible mutants *)
retentionVsTime, replicationVsTime,  (* the resident's traits over time (as evolution occurs and the resident changes) *)
timeVsTime, timeStepStart, (* troubleshooting variables, for measuring how long each step of the simulation takes *)
cyclingVsTime, (* troubleshooting variable, stores whether cycling occurred *)
maxMutationSteps = OptionValue[maxMutations]}, (* the maximum number of times we'll simulate a mutant fixing in the population *)

(* Initial lists that will hold the resident's traits over time *)
retentionVsTime = Table["NA", maxMutationSteps + 1];
replicationVsTime = Table["NA", maxMutationSteps + 1];

(* Initialize list of runtimes and list of whether cycling occurred (troubleshooting *)
timeVsTime = Table["NA", maxMutationSteps + 1];
cyclingVsTime = Table["NA", maxMutationSteps + 1];

(* Initialize the simulation with the initial resident's retention and replication rates (default 0, purely heterotrophic) *)
residentRetention = OptionValue[initialRetention];
residentReplication = OptionValue[initialReplication];
retentionVsTime[[1]] = residentRetention;
replicationVsTime[[1]] = residentReplication;
timeVsTime[[1]] = 0;

(* Introduce mutants into the population and see whether they invade and become the new resident. *)
(* Do is a for loop. We'll be looping over i = 2:1001, which comes at the very end of all this code, because Mathematica is weird about for loop syntax *)
Do[
	(* **** Find resident population properties before mutant invasion **** *)
	timeStepStart = SessionTime[];
	
	(* Find the resident's ecological dynamics *)
	{residentEcology, extinction} = residentEcoDynamics[residentRetention, residentReplication, surfaceLight, 1, 1, 1, 
		Evaluate[FilterRules[{opts, maxTime -> 8000}, Options[residentEcoDynamics]]]];
		
	(* If extinction comes back false, check that the resident population didn't drop to <= 10^-6 (we also consider that extinction) *)
	If[!(extinction),
		minResident = NMinimize[{(h[t] + p[t])/. residentEcology, 3000 <= t <= 8000}, t, Evaluate[FilterRules[{opts}, Options[NMinimize]]]][[1]];
		extinction = (minResident <= 10^-6)];
		
	(* If the resident cannot survive, stop! *)
	If[extinction,
	Print["Resident with (retention = " <> ToString[residentRetention] <> ", replication = " <> ToString[residentReplication] <> 
		") went extinct. No further mutations introduced."];
	Break[]];

	(* Find the maximum resident population size. We will use this to see if the population is cycling *)
	maxResident = NMaximize[{(h[t] + p[t])/. residentEcology, 3000 <= t <= 8000}, t, Evaluate[FilterRules[{opts}, Options[NMaximize]]]][[1]];
	
	(* Check for cycling. Find the period of the cycle if it is present *)
	If[minResident <= 0.99 * maxResident, 
		(* The difference between minimum and maximum population sizes suggests cycling. *)
			cyclingVsTime[[i - 1]] = True;
			(* Estimate the period of cycling using the heterotroph population *)
			period = interpolatingFunctionPeriod[(h /. residentEcology)],
		(* No cycling in the resident population *)
			cyclingVsTime[[i - 1]] = False;
			(* Set "period" to an arbitrary number (50 days). *)
			(* We need "period" to be > 0 to calculate mutant growth, but the exact number doesn't matter if there is no cycling. *)
			(* We chose 50 days becausie it's similar to the periods for cycling populations *)
			period = 50
		];
	
	(* **** Simulation mutant invasion **** *)
	
	(* Create a list of possible mutants *)
	mutants = Map[{residentRetention, residentReplication} + #&, {{-1/100, 0}, {1/100, 0}, {0, -1/100}, {0, 1/100}}];

	(* Remove mutants that cannot arise (retention or replication is below 0 or greater than 1) *)
	mutants = Select[mutants, (0 <= #[[1]] <= 1) && (0 <= #[[2]] <= 1)&];
	
	(* Calculate the fixation probability of each mutant. *)
	fixationProbabilities = Map[ (* Apply the following function to all mutants *)
		Module[{mutantDynamics, mutantExtinction, minMutant, selectionCoefficient, prFix}, (* Variables we'll use to calculate fixation prob. *)
		
			(* Find the mutant's ecological dynamics when rare *)
			{mutantDynamics, mutantExtinction} = rareMutantEcoDynamics[(w /. residentEcology),  (h /. residentEcology), (p /. residentEcology), 
				#[[1]] (* Mutant retention *), #[[2]] (* Mutant replication *), surfaceLight, (h[4000] /. residentEcology), (p[4000] /. residentEcology), 
				4000 (* start time for finding eco. dynamics *), 8000 (* end time *), Evaluate[FilterRules[{opts}, Options[rareMutantEcoDynamics]]]];
				
			(* If extinction comes back false, check that mutant population is positive and non-tiny (greater than 10^-10)  *)
			If[!(mutantExtinction),
				minMutant = NMinimize[{(h2[t] + p2[t])/. mutantDynamics, 4000 <= t <= 8000}, t, Evaluate[FilterRules[{opts}, Options[NMinimize]]]][[1]];
				mutantExtinction = (minMutant <= 10^-10)];

			If[mutantExtinction,
				(* If the mutant went extinct, the probability that it fixes is 0 *)
					prFix = 0,
				(* If the mutant didn't go extinct... *)
					(*  Use mutant's growth rate to approximate the selection coefficient *)
					selectionCoefficient = (( (h2[8000] + p2[8000]) / (h2[8000 - period] + p2[8000 - period]) /. mutantDynamics) - 1) / period;
					(* Use the selection coefficient to approximate fixation probability *)
					prFix = 1 - Exp[-2 * selectionCoefficient];
					(* If the fixation probability is negative, set it to 0 *)
					prFix = Max[prFix, 0]
				];
				
			(* Return the fixation probability *)
			prFix]& (* end of module and end of function *)
			
		(* Map the above function over all the mutants *)
		, mutants];
		
		
	(* If all the fixation probabilities are 0, stop the simulation, the resident is locally uninvadable *)
	If[Total[fixationProbabilities] == 0,
		Print["Resident (retention = " <> ToString[N[residentRetention]] <> ", replication = " <> ToString[N[residentReplication]] <> ") is locally uninvadable"];
	
		(* Fill the last entries of the simulation with the uninvadable resident *)
		retentionVsTime[[i;;maxMutationSteps + 1]] = residentRetention;
		replicationVsTime[[i;;maxMutationSteps + 1]] = residentReplication;
	
		(* End the simulation (get out of the for loop *)
		Break[]];

	(* Randomly choose a mutation to fix, based on the fixation probabilities *)
	{residentRetention, residentReplication} = RandomChoice[fixationProbabilities/Total[fixationProbabilities] -> mutants];
	
	(* Record the new resident's trait values *)
	retentionVsTime[[i]] = residentRetention;
	replicationVsTime[[i]] = residentReplication;
	
	(* Record the time this step of the simulation took *)
	timeVsTime[[i]] = SessionTime[] - timeStepStart,
		
	(* End of the for-loop. We write the variable to iterate over (i) and the values to iterate over, 2 to maxMutationSteps + 1 *)
	{i, 2, maxMutationSteps + 1}];

(* Return retention and replication values over time *)
Return[{retentionVsTime, replicationVsTime, timeVsTime, cyclingVsTime}]]
