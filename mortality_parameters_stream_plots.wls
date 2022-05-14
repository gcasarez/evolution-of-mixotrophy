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
gammaP = 2/3;(*2/3;*)(*2.9/3;*)
sigmaP = 0;
kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;
hW = 50;(*50;*)
hP = hP2 = 50;(*10;*)(*50;*)
lW = 0.5;
(*lH = lH2 = 0.25;
lP = lP2 = 0.25;*)
dMax = 0.3;(*0.3;*)
sigmaD = 0; (* No more decay trade-off *)
a = a2 = 0.15;
gammaA = 0; (* We want this to be 0 consistently, because we've decided the phototroph doesn't eat *)
(*sigmaA = 0;*) (* 0 for no attack rate trade-off, 1 for attack rate trade-off *)
c = 0.1;
(*gammaC = 1;*) (* 1 for no conversion efficiency trade-off, 0.7 for conversion efficiency trade-off *)
sigmaC = 0; (* Added sigmaC because it looks like it's not defined here *)
(*iIn = 25; *)(* Incident light *)


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
(*Phototroph mortality rate is higher*)


lH = lH2 = 0.15;
lP = lP2 = 0.35;


light = 50;


(* ::Subsection::Closed:: *)
(*Attack & conversion trade-offs*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARhighlP = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARhighlP = Flatten[ecoDynamicsCEARhighlP, {{1,2}, {3}}];


(* ::Text:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


extinctionPointsCEARhighlP = Select[ecoDynamicsCEARhighlP, (* Extinction *) (#[[4]] == True) || (* Population < 10^-6 *) ((#[[6]] < 10^-6) && (#[[7]] < 10^-6))&];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEARhighlP = Select[ecoDynamicsCEARhighlP, (#[[4]] == False) && ((#[[6]] >= 10^-6) || ( #[[7]] >= 10^-6))&];


(* ::Text:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEARhighlP = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEARhighlP];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and down*)


Select[derivsForPlotsCEARhighlP, #[[1]] == {1, 0}&][[All, 2]]


Export["higherPhototrophMortalityCEAR.wdx", derivsForPlotsCEARhighlP]


(* ::Subsection::Closed:: *)
(*Conversion trade-off only*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEhighlP = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEhighlP = Flatten[ecoDynamicsCEhighlP, {{1,2}, {3}}];


(* ::Text:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


extinctionPointsCEhighlP = Select[ecoDynamicsCEhighlP, (* Extinction *) (#[[4]] == True) || (* Population < 10^-6 *) ((#[[6]] < 10^-6) && (#[[7]] < 10^-6))&];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEhighlP = Select[ecoDynamicsCEhighlP, (#[[4]] == False) && ((#[[6]] >= 10^-6) || ( #[[7]] >= 10^-6))&];


(* ::Text:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEhighlP = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEhighlP];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and down*)


Select[derivsForPlotsCEhighlP, #[[1]] == {1, 0}&][[All, 2]]


Export["higherPhototrophMortalityCE.wdx", derivsForPlotsCEhighlP]


(* ::Section:: *)
(*Heterotroph mortality rate is higher*)


lH = lH2 = 0.35;
lP = lP2 = 0.15;


(*light = 50;*)


(* ::Subsection::Closed:: *)
(*Attack & conversion trade-offs*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARhighlH = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARhighlH = Flatten[ecoDynamicsCEARhighlH, {{1,2}, {3}}];


(* ::Text:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


extinctionPointsCEARhighlH = Select[ecoDynamicsCEARhighlH, (* Extinction *) (#[[4]] == True) || (* Population < 10^-6 *) ((#[[6]] < 10^-6) && (#[[7]] < 10^-6))&];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEARhighlH = Select[ecoDynamicsCEARhighlH, (#[[4]] == False) && ((#[[6]] >= 10^-6) || ( #[[7]] >= 10^-6))&];


(* ::Text:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEARhighlH = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEARhighlH];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and up.*)


Select[derivsForPlotsCEARhighlH, #[[1]] == {1, 0}&][[All, 2]]


Export["higherHeterotrophMortalityCEAR.wdx", derivsForPlotsCEARhighlH]


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEhighlH = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEhighlH = Flatten[ecoDynamicsCEhighlH, {{1,2}, {3}}];


(* ::Text:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


extinctionPointsCEhighlH = Select[ecoDynamicsCEhighlH, (* Extinction *) (#[[4]] == True) || (* Population < 10^-6 *) ((#[[6]] < 10^-6) && (#[[7]] < 10^-6))&];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEhighlH = Select[ecoDynamicsCEhighlH, (#[[4]] == False) && ((#[[6]] >= 10^-6) || ( #[[7]] >= 10^-6))&];


(* ::Text:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEhighlH = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEhighlH];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and up*)


Select[derivsForPlotsCEhighlH, #[[1]] == {1, 0}&][[All, 2]]


Export["higherHeterotrophMortalityCE.wdx", derivsForPlotsCEhighlH]


(* ::Section:: *)
(*Default mortality rates (lH = lP)*)


lH = lH2 = 0.25;
lP = lP2 = 0.25;


(* ::Subsection::Closed:: *)
(*Attack & conversion trade-offs*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARdefault = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARdefault = Flatten[ecoDynamicsCEARdefault, {{1,2}, {3}}];


(* ::Text:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


extinctionPointsCEARdefault = Select[ecoDynamicsCEARdefault, (* Extinction *) (#[[4]] == True) || (* Population < 10^-6 *) ((#[[6]] < 10^-6) && (#[[7]] < 10^-6))&];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEARdefault = Select[ecoDynamicsCEARdefault, (#[[4]] == False) && ((#[[6]] >= 10^-6) || ( #[[7]] >= 10^-6))&];


(* ::Text:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEARdefault = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEARdefault];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and down.*)


Select[derivsForPlotsCEARdefault, #[[1]] == {1, 0}&][[All, 2]]


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Text:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEdefault = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEdefault = Flatten[ecoDynamicsCEdefault, {{1,2}, {3}}];


(* ::Text:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


extinctionPointsCEdefault = Select[ecoDynamicsCEdefault, (* Extinction *) (#[[4]] == True) || (* Population < 10^-6 *) ((#[[6]] < 10^-6) && (#[[7]] < 10^-6))&];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEdefault = Select[ecoDynamicsCEdefault, (#[[4]] == False) && ((#[[6]] >= 10^-6) || ( #[[7]] >= 10^-6))&];


(* ::Text:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEdefault = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEdefault];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and up*)


Select[derivsForPlotsCEdefault, #[[1]] == {1, 0}&][[All, 2]]


(* ::Section:: *)
(*Plots*)


(* ::Text:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions -> {"OperatorSubstitution"->False}}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions->{"OperatorSubstitution"->False}}];


(* ::Text:: *)
(*Color scheme (prey density) and legend*)


maxPrey = Max[{noExtinctionPointsCEhighlP[[All, 5]], noExtinctionPointsCEdefault[[All, 5]], noExtinctionPointsCEhighlH[[All, 5]],
noExtinctionPointsCEARhighlP[[All, 5]], noExtinctionPointsCEARdefault[[All, 5]], noExtinctionPointsCEARhighlH[[All, 5]]}];


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


(*BarLegend[{preyDensityColors[#]&, {0, maxPrey}}]*)
preyDensityLegend = BarLegend[{preyDensityColors[Exp[#] - 1]&, {Log[1], Log[1 + maxPrey]}}, 
Ticks -> Table[{Log[i + 1], ToString[i]}, {i, {0, 1, 2.5, 5, 10, 15, 25, 40}}],
LegendLayout -> "Column", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {25, 450}];


(* ::Text:: *)
(*Set up plot labels*)


plotLabels = (*{"(A) High phototroph mortality\nLow heterotroph mortality\nNo tradeoff", 
"(B) Equal heterotroph and\nphototroph mortality\nNo tradeoff",
"(C) Low phototroph mortality\nHigh heterotroph mortality\nNo tradeoff",
"(D) High phototroph mortality\nLow heterotroph mortality\nAttack tradeoff", 
"(E) Equal heterotroph and\nphototroph mortality\nAttack tradeoff",
"(F) Low phototroph mortality\nHigh heterotroph mortality\nAttack tradeoff"};*)
Table[Pane[Style[letter, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}], ImageSize -> 300],
{letter, {"(A)", "(B)", "(C)", "(D)", "(E)", "(F)"}}];


(* ::Text:: *)
(*Mark points where the population went extinct with black boxes*)


extinctionGraphics = 
Table[Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPoints[[All, 2;;3]]]], 
	{extinctionPoints, 
	{extinctionPointsCEhighlP, extinctionPointsCEdefault, extinctionPointsCEhighlH, 
	extinctionPointsCEARhighlP, extinctionPointsCEARdefault, extinctionPointsCEARhighlH}}];


(* ::Text:: *)
(*Collect derivatives in a list*)


derivsForPlots = {derivsForPlotsCEhighlP, derivsForPlotsCEdefault, derivsForPlotsCEhighlH,
	derivsForPlotsCEARhighlP, derivsForPlotsCEARdefault, derivsForPlotsCEARhighlH};


(* ::Text:: *)
(*Make plots*)


plots = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlots[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> plotLabels[[i]], ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlots[[i]]], StreamStyle -> Black, StreamColorFunction -> None],
	
	extinctionGraphics[[i]]}],
{i, Length[derivsForPlots]}]


(* ::Text:: *)
(*Combine plots and information for export*)


Row[{Grid[

	(* Row labels *)
	Join[Map[Rotate[Style[#, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}], 90Degree]&, 
		{{""}, {"No tradeoff"}, {"Attack tradeoff"}}, {2}],
		
	(* Column labels *)
	Join[Map[Style[#, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}]&, 
		{{"High phototroph mortality, low heterotroph mortality", "Equal phototroph and heterotroph mortality", "Low phototroph mortality, high heterotroph mortality"}},
		{2}],
	 
	(* Some extremely fancy rasterization of the plots without rasterizing their labels (removes the labels, raterizes, and re-adds the labels *)
	Map[
		Labeled[Rasterize[Show[#, FrameLabel -> None, PlotLabel -> None, ImagePadding -> 1](*, RasterSize \[Rule] 500*), ImageSize -> Medium], 
		Append[(FrameLabel /. Options[#])[[All, 1]], Style[(PlotLabel /. Options[#]), 14, TextAlignment -> Center]], {Left, Bottom, Top},
		LabelStyle -> Append[(LabelStyle /. Options[#]), (TextAlignment -> Center)], RotateLabel -> True]&,
	{plots[[1;;3]], plots[[4;;6]]}, {2}]
		(* close out the two joinings of row and column labels *)], 2], 
Spacings -> {2, 3}], 

	(* Legend *)
	Column[{
		BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Column", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {25, 450}],
		PointLegend[{Black, Black}, {"Direction of selection", "Consumer extinction"}, 
			LegendMarkers -> {Graphics[{Arrowheads[0.5], Arrow[{{0, 0}, {1, 0}}]}], Graphics[Rectangle[{0, 0}]]},
			LegendMarkerSize -> {{75, 75}, {20, 20}}]},
	Frame -> True, Alignment -> {Center, Center}]}, 
Spacer[50]]

Export["stream_plots_varying_mortality.pdf", %]






