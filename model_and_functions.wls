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


Options[simulateEvolution] = Join[{attackTradeOff -> False, conversionTradeOff -> False, maxMutations -> 1000, initialRetention -> 0, initialReplication -> 0}, 
	Options[NDSolve], Options[NMinimize], Options[FindRoot]];
Options[minimumMixotroph] = Join[Options[minimumMixotroph], Options[NMinimize]];


Options[fixationFromMeanOffspring] = Join[{uHguess -> 0, uPguess -> 0}, Options[FindRoot]];
Options[getFixationProbability] = Join[Options[meanOffspring], Options[fixationFromMeanOffspring]];


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


(* ::Section::Closed:: *)
(*Simulate evolution*)


(* ::Subsection:: *)
(*Function to constrain mutations to biologically possible values*)


(* ::Text:: *)
(*Simple function to take a list of mutant retention and replication values and makes those < 0 equal to 0, and those > 1 equal to 1*)


constrainMutation = Map[Max[Min[#, 1], 0]&];


(* ::Subsection:: *)
(*Function to simulation population evolution*)


(* ::Text:: *)
(*A function to simulate  mutations arising in a population of mixotrophs, and fixing if the mutants appear to be outcompeting the resident*)


simulateEvolution[surfaceLight_, opts:OptionsPattern[]] :=
Module[{
(* Variables we'll use inside this function *)
residentRetention, residentReplication, mutantRetention, mutantReplication, (* the resident and mutant's traits *)
residentEcology, mutantAndResidentEcology, (* the ecological dynamics of mutants and residents *)
extinction, (* will be True if the entire mixotroph population went extinct *)
averagePrey, averageHeterotroph, averagePhototroph, (* average prey, heterotrophs, and phototrophs when the resident is alone *)
averageResidentMixotroph, averageMutantMixotroph, (* average mixotroph population sizes when the mutant is present *)
initialPrey, initialResidentHeterotroph, initialResidentPhototroph, initialMutantHeterotroph, initialMutantPhototroph, (* initial conditions for ecological simulations *)
mutation, (* size & direction of a mutation *)
mutantProportion, (* proportion of the mixotroph population that are mutants *)
retentionVsTime = Table["NA", 1001], replicationVsTime = Table["NA", 1001]},  (* the resident's traits over time (as evolution occurs and the resident changes) *)

(* Initialize the simulation with a resident that is purely heterotrophic (rentention = replication = 0) *)
residentRetention = residentReplication = 0;
retentionVsTime[[1]] = residentRetention;
replicationVsTime[[1]] = residentReplication;

(* Find the initial resident's ecological dynamics *)
{residentEcology, extinction} = residentEcoDynamics[residentRetention, residentReplication, surfaceLight, 1, 1, 1, opts];

(* If the resident cannot survive, stop! *)
If[extinction,
Print["Mixotroph with (phi = 0, eta = 0) does not survive when surface light is " <> ToString[surfaceLight] <> " \[Mu]mol photons/m^2/s"];
Return[{retentionVsTime, replicationVsTime}]];

(* Find the "average" ecological state with just the resident (integrate each of the state variables over the last 1000 days and divide by 1000 to get the average *)
{averagePrey, averageHeterotroph, averagePhototroph} = Table[
	NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000), 
	{stateVar, {w[t], h[t], p[t]}}];
	
(* Set the initial values of the resident and mutant populations, so that the mutant is 1% of the resident population *)
initialPrey = averagePrey;
{initialResidentHeterotroph, initialResidentPhototroph} = 0.99 * {averageHeterotroph, averagePhototroph};
{initialMutantHeterotroph, initialMutantPhototroph} = 0.01 * {averageHeterotroph, averagePhototroph};



(* Introduce mutations and see whether the mutants can invade. (Do is a for loop. We'll be looping over i = 2:1001, which comes at the very end of all this code, because Mathematica is weird about for loop syntax) *)
Do[
	(* Obtain a mutant! Start by initializing a vector of mutations to retention and replication *)
	mutation = {0, 0};

	(* Randomly choose one of the two traits (retention & replication) to mutate. Draw that mutation from a normal distribution with mean 0, standard deviation 0.005 *)
	mutation[[RandomChoice[{1, 2}]]] = Random[NormalDistribution[0, 0.01]]; 

	(* Mutate retention/replication, making sure that they stay between 0 and 1 *)
	{mutantRetention, mutantReplication} = constrainMutation[{residentRetention, residentReplication} + mutation];
	
	(* Get the ecological dynamics of the mutant and resident *)
	{mutantAndResidentEcology, extinction} = mutantAndResidentEcoDynamics[residentRetention, residentReplication, mutantRetention, mutantReplication, surfaceLight, 
	initialPrey, initialResidentHeterotroph, initialResidentPhototroph, initialMutantHeterotroph, initialMutantPhototroph, opts];
	
	(* If the resident and mutant go extinct when they are together, stop the analysis *)
	If[extinction,
		Print["Resident and mutant go extinct when together. Resident (retention, replication) = (" <> ToString[residentRetention] <> ", " <> ToString[residentReplication] <> 
		"), mutant (retention, replication) = (" <> ToString[mutantRetention] <> ", " <> ToString[mutantReplication] <> ")"];
		Break[]];
	
	(* Find the "average" ecological state with the mutant and resident *)
	(* This is slightly different from how we previously did it *)
	(* Previously, we looked at the last time step, and if the mutant phototroph was > 1.02 * initial mutant phototroph AND mutant heterotroph was > 1.02 * mutant heterotroph, we said the mutant invaded *)
	{averageResidentMixotroph, averageMutantMixotroph} = Table[
		NIntegrate[stateVar/. mutantAndResidentEcology, {t, 3000, 4000}]/(4000 - 3000), 
		{stateVar, {h[t] + p[t], h2[t] + p2[t]}}];
		
		
	(* Calculate the proportion of the population that are mutants *)
	mutantProportion = (averageMutantMixotroph) / (averageMutantMixotroph + averageResidentMixotroph);

	(* The mutant invades (becomes the resident) if it has increased beyond its initial 1% of the population *)
	If[mutantProportion > 0.01,
		{residentRetention, residentReplication} = {mutantRetention, mutantReplication};
	
		(* Find the new resident's ecological state *)
		{residentEcology, extinction} = residentEcoDynamics[residentRetention, residentReplication, surfaceLight, averagePrey, averageHeterotroph, averagePhototroph, opts];
		
		(* If the mutant goes extinct on its own, record the trait values where this happened and stop the analysis *)
		If[extinction, 
			Print["New resident went extinct; retention = " <> ToString[residentRetention] <> ", replication = " <> ToString[residentReplication]]; 
			retentionVsTime[[i]] = residentRetention;
			replicationVsTime[[i]] = residentReplication;
			Break[]];
			
		(* Find the averages of the state variables *)
		{averagePrey, averageHeterotroph, averagePhototroph} = Table[
			NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000), 
			{stateVar, {w[t], h[t], p[t]}}];
	
		(* Set the initial values of the resident and mutant populations, so that the new mutant will be 1% of the new resident population *)
		initialPrey = averagePrey;
		{initialResidentHeterotroph, initialResidentPhototroph} = 0.99 * {averageHeterotroph, averagePhototroph};
		{initialMutantHeterotroph, initialMutantPhototroph} = 0.01 * {averageHeterotroph, averagePhototroph}];
		
	(* Record the current resident's trait values *)
	retentionVsTime[[i]] = residentRetention;
	replicationVsTime[[i]] = residentReplication, (* comma here to signal the end of the for loop *)
		
(* End of the for-loop. We write the variable to iterate over (i) and the values to iterate over, 2 to 501 (we introduce 500 mutations) *)
{i, 2, 1001}];

(* Return retention and replication values over time *)
Return[{retentionVsTime, replicationVsTime}]]


(* ::Section:: *)
(*Simulate evolution accounting for fixation probability*)


(* ::Text:: *)
(*Assume stage structure & no cycling*)


(* ::Subsection:: *)
(*Helper functions*)


(* ::Subsubsection::Closed:: *)
(*Mean offspring*)


(* ::Text:: *)
(*Function to find the mean offspring in each state (heterotrophic, phototrophic) produced by a mixotroph in each state.*)
(*This function works for residents or mutants (use it for residents to normalize mutant reproduction).*)


meanOffspring[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), 
averagePrey_ (* W[0] *), averageHeterotroph_ (* H[0] *), averagePhototroph_ (* P[0] *), lightAbsorption_ (* kappa *), opts:OptionsPattern[]] :=
Module[{sigmaA, gammaC, attackRate, d, pP, lambdaHH, lambdaPH, lambdaHP, lambdaPP},

(* Implement trade-offs *)
sigmaA = If[OptionValue[attackTradeOff], 1, 0]; (* 1 for attack rate trade-off, 0 for no attack rate trade-off *)
gammaC = If[OptionValue[conversionTradeOff], 0.7, 1]; (* 0.7 for conversion efficiency trade-off, 1 for no conversion efficiency trade-off *)

(* Attack rate *)
attackRate = a * (1 - sigmaA * (retention * replication));

(* Plastid decay rate *)
d = dMax * (1 - sigmaD * replication);

(* Phototroph's photosynthetic rate *)
pP = gammaP * pW;

(* Mean number each type of "offspring" (including self) produced *)
lambdaHH = 1 + (1 - retention) * c * attackRate * averagePrey - lH; (* heterotroph producing heterotrophic offspring *)
lambdaPH = (1 - replication) * (1/gammaC) * pP / lightAbsorption * Log[(hP + surfaceLight) / (hP + surfaceLight * Exp[-1 * lightAbsorption])] + d; (* phototroph producing heterotrophic offspring *)
lambdaHP = retention * attackRate * averagePrey; (* heterotroph producing phototrophic offspring *)
lambdaPP = 1 + replication * pP / lightAbsorption * Log[(hP + surfaceLight) / (hP + surfaceLight * Exp[-1 * lightAbsorption])] - d - lP; (* phototroph producing phototrophic offspring *)

(* Return the mean number of offspring produced, in matrix form, such that the new number of {H, P} = matrix * {H, P} *)
(* THE ORDER IS NOT A MISTAKE!! lambdaPH belongs in the first row, lambdaHP in the second *)
Return[{{lambdaHH, lambdaPH}, {lambdaHP, lambdaPP}}]
]


(* ::Subsubsection::Closed:: *)
(*Calculate fixation probability given mean offspring*)


(* ::Text:: *)
(*Function to  calculate the probability of fixation given the mean offspring produced by each state and the average resident population in each state.*)
(*This function finds the extinction probabilities and then passes them off to fixationFromExtinctionProbability (below), which checks them and returns fixation probabilities.*)


fixationFromMeanOffspring[{{lambdaHH_, lambdaPH_}, {lambdaHP_, lambdaPP_}}, {averageHeterotroph_, averagePhototroph_},  opts:OptionsPattern[]] :=
Module[{uHeqn, uPeqn},
uHeqn = (uH == Exp[-(lambdaHH * (1 - uH) + lambdaHP * (1 - uP))]); (* Pr(a mutant heterotroph & all its descendents go extinct) *)
uPeqn = (uP == Exp[-(lambdaPH * (1 - uH) + lambdaPP * (1 - uP))]); (* Pr(a mutant phototroph & all descendents go extinct) *)

(* Solve for extinction probabilities *)
extinctionProbs = FindRoot[{uHeqn, uPeqn}, {{uH, OptionValue[uHguess]}, {uP, OptionValue[uPguess]}}, Evaluate[FilterRules[{opts}, Options[FindRoot]]]];

(* Calculate expected fixation probability *)
Return[fixationFromExtinctionProbability[{uH, uP} /. extinctionProbs, {averageHeterotroph, averagePhototroph}]]
]


(* ::Subsubsection::Closed:: *)
(*Calculate fixation probability given extinction probabilities*)


(* ::Text:: *)
(*This function takes in the extinction probabilities of a single heterotroph and phototroph, as well as the likelihood they initially arise in each state, and calculates the expected probability of fixation.*)


fixationFromExtinctionProbability::imaginaryNumbers = "One or both of the extinction probabilities is a complex number. Fixation probability cannot be calculated.";
fixationFromExtinctionProbability::notProbability = "Both extinction probabilities are real numbers, but one or both is outside of the range [0, 1]. Fixation probability cannot be calculated.";
fixationFromExtinctionProbability::notPopSizes = "Total population is <= 0 or at least one subpopulation is negative. Fixation probability cannot be calculated.";


fixationFromExtinctionProbability[{uH_, uP_}, {h_, p_}] := 
Module[{},
(* Check whether uH and uP are actual extinction probabilities (real, between 0 and 1) *)
If[Im[{uH, uP}] != {0, 0},
	Message[fixationFromExtinctionProbability::imaginaryNumbers];
	Return["Unknown"]];
	
If[(uH < 0) || (Chop[uH - 1] > 0) || (uP < 0) || (Chop[uP - 1] > 0),
	Message[fixationFromExtinctionProbability::notProbability];
	Print[{uH, uP}];
	Return["Unknown"]];
	
(* Check whether h and p are real population sizes (nonnegative and sum to a positive number) *)
If[(h < 0)|| (p < 0) || (h + p <= 0),
	Message[fixationFromExtinctionProbability::notPopSizes];
	Return["Unknown"]];
	
(* If the extinction probabilities are real probabilities, find the expected fixation probability *)
(* We take the expectation across the probability that a mutant arises in the heterotrophic or phototrophic state *)
Return[h / (h + p) * Chop[(1 - uH)] + p / (h + p) * Chop[(1 - uP)]]]


(* ::Subsection::Closed:: *)
(*Function to calculate the fixation probability of a mutant, given resident ecology*)


(* ::Text:: *)
(*Function to calculate the fixation probability of a mutant from start to finish. Requires environmental variables (including resident population state) and mutant traits.*)


getFixationProbability[retention_ (* phi *), replication_ (* eta *), surfaceLight_ (* iIn *), 
averagePrey_ (* W[0] *), averageHeterotroph_ (* H[0] *), averagePhototroph_ (* P[0] *), lightAbsorption_ (* kappa *), residentGrowthRate_, opts:OptionsPattern[]] :=
Module[{sigmaA, gammaC, attackRate, (*lightAbsorption,*) pP, offspring},

(* Implement trade-offs *)
sigmaA = If[OptionValue[attackTradeOff], 1, 0]; (* 1 for attack rate trade-off, 0 for no attack rate trade-off *)
gammaC = If[OptionValue[conversionTradeOff], 0.7, 1]; (* 0.7 for conversion efficiency trade-off, 1 for no conversion efficiency trade-off *)

(* Attack rate *)
attackRate = a * (1 - sigmaA * (retention * replication));

(* Plastid decay rate *)
d = dMax * (1 - sigmaD * replication);

(* Phototroph's photosynthetic rate *)
pP = gammaP * pW;

(* Mean number each type of "offspring" (including self) produced *)
offspring = meanOffspring[retention, replication, surfaceLight, averagePrey, averageHeterotroph, averagePhototroph, lightAbsorption, 
	Evaluate[FilterRules[{opts}, Options[meanOffspring]]]];

(* Normalize offspring by resident growth rate *)
offspring = offspring / residentGrowthRate;

(* Calculate expected fixation probability *)
Return[fixationFromMeanOffspring[offspring, {averageHeterotroph, averagePhototroph}, Evaluate[FilterRules[{opts}, Options[fixationFromMeanOffspring]]]]]
]


(* ::Subsection::Closed:: *)
(*Function that uses the fixation probability to model evolution*)


simulateEvolution[surfaceLight_, opts:OptionsPattern[]] :=
Module[{
(* Variables we'll use inside this function *)
residentRetention, residentReplication, (* the resident's traits *)
residentEcology, (* the ecological dynamics with just the resident present *)
extinction, (* will be True if the entire resident population goes extinct *)
minResident, (* minimum resident population size *)
averagePrey, averageHeterotroph, averagePhototroph, (* average prey, heterotrophs, phototrophs, and light absorption when mutant is absent (or rare) *)
averageAbsorption, (* average light absorption when mutant is absent (or rare) *)
residentProjectionMatrix, (* matrix giving the average heterotrophs & phototrophs produced by each type of resident *)
residentGrowthRate, (* resident population growth rate *)
mutants, (* list of possible mutants that might arise *)
prFix, (* list of fixation probabilities of possible mutants *)
retentionVsTime, replicationVsTime,  (* the resident's traits over time (as evolution occurs and the resident changes) *)
timePrev, timeRes, timeEco, timeAvgPop, timeMinPop, timeFix, timeVsTime, (* troubleshooting variables, store the time at different points in the simulation *)
cyclingVsTime, (* "troubleshooting" variable, stores whether cycling occurred *)
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
	timeRes = timePrev = SessionTime[];
	
	(* Find the initial resident's ecological dynamics *)
	{residentEcology, extinction} = residentEcoDynamics[residentRetention, residentReplication, surfaceLight, 1, 1, 1, 
		Evaluate[FilterRules[{opts}, Options[residentEcoDynamics]]]];
		
	(*Print["Ecological dynamics: " <> ToString[SessionTime[]]];*)
	timeEco = SessionTime[] - timePrev;
	timePrev = timePrev + timeEco;
	
	(* If extinction comes back false, check that the resident population didn't drop to <= 10^-6 (we also consider that extinction) *)
	If[!(extinction),
		minResident = NMinimize[{(h[t] + p[t])/. residentEcology, 3000 <= t <= 4000}, t, Evaluate[FilterRules[{opts}, Options[NMinimize]]]][[1]];
		extinction = (minResident <= 10^-6)];
		
	(*Print["Min pop. size: " <> ToString[SessionTime[]]];*)
	timeMinPop = SessionTime[] - timePrev;
	timePrev = timeMinPop + timePrev;

	(* If the resident cannot survive, stop! *)
	If[extinction,
	Print["Resident with (retention = " <> ToString[residentRetention] <> ", replication = " <> ToString[residentReplication] <> 
		") went extinct. No further mutations introduced."];
	Break[]];

	(* Find the "average" ecological state with just the resident (integrate each of the state variables over the last 1000 days and divide by 1000 to get the average *)
	{averagePrey, averageHeterotroph, averagePhototroph, averageAbsorption} = Table[
		NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000), 
		{stateVar, {w[t], h[t], p[t], kappa[kW, kH, kP, kH, kP, w[t], h[t], p[t], 0, 0]}}];
	
	(* Just for now -- debugging *)
	(* averageAbsorption = kappa[kW, kH, kP, kH, kP, averagePrey, averageHeterotroph, averagePhototroph, 0, 0]; *)
	
	(*Print["Avg pop size: " <> ToString[SessionTime[]]];*)
	timeAvgPop = SessionTime[] - timePrev;
	timePrev = timePrev + timeAvgPop;
		
	(* Turn tiny negative numbers into 0s *)
	{averagePrey, averageHeterotroph, averagePhototroph, averageAbsorption} = Chop[{averagePrey, averageHeterotroph, averagePhototroph, averageAbsorption}];
	
	(* Check for cycling (debugging) *)
	If[minResident < 0.99 * (averageHeterotroph + averagePhototroph), 
		cyclingVsTime[[i - 1]] = True
		(*Print["Cycling at (retention = " <> ToString[residentRetention] <> ", replication = " <> ToString[residentReplication] <> ")"]*),
		cyclingVsTime[[i - 1]] = False];
	
	(* Get the population growth rate to normalize the mutant rates by (approximation, assumes constant pop. size and poisson-distributed number of offspring) *)
	residentProjectionMatrix = meanOffspring[residentRetention, residentReplication, surfaceLight, 
		averagePrey, averageHeterotroph, averagePhototroph, averageAbsorption, Evaluate[FilterRules[{opts}, Options[meanOffspring]]]];
	residentGrowthRate = Eigenvalues[residentProjectionMatrix, 1][[1]];
	
	(* **** Simulation mutant invasion **** *)
	
	(* Create a list of possible mutants *)
	mutants = Map[{residentRetention, residentReplication} + #&, {{-0.01, 0}, {0.01, 0}, {0, -0.01}, {0, 0.01}}];

	(* Remove mutants that cannot arise (retention or replication is below 0 or greater than 1) *)
	mutants = Select[mutants, (0 <= #[[1]] <= 1) && (0 <= #[[2]] <= 1)&];

	(* Find the fixation probabilities of each mutation *)
	prFix = Map[getFixationProbability[#[[1]], #[[2]], surfaceLight, averagePrey, averageHeterotroph, averagePhototroph, averageAbsorption, residentGrowthRate,
		Evaluate[FilterRules[{opts}, Options[getFixationProbability]]]]&, mutants];
		
	(*Print["Fixation probability: " <> ToString[SessionTime[]]];*)
	timeFix = SessionTime[] - timePrev;
	timePrev = timePrev + timeFix;
		
	(* Check that prFix = Unknown is not getting into the output *)
	If[AnyTrue[prFix, # == "Unknown"&], Print["??? "]; Print[prFix]; Break[]];

	(* If all the fixation probabilities are 0, stop the simulation, the resident is locally uninvadable *)
	If[Total[prFix] == 0,
		Print[prFix];
		Print["Resident is locally uninvadable"];
		Print[residentRetention];
		Print[residentReplication];
	
		(* Fill the last entries of the simulation with the uninvadable resident *)
		retentionVsTime[[i;;maxMutationSteps + 1]] = residentRetention;
		replicationVsTime[[i;;maxMutationSteps + 1]] = residentReplication;
	
		(* End the simulation (get out of the for loop *)
		Break[]];

	(* Randomly choose a mutation to fix, based on the fixation probabilities *)
	{residentRetention, residentReplication} = RandomChoice[prFix/Total[prFix] -> mutants];
	
	(* Record the new resident's trait values *)
	retentionVsTime[[i]] = residentRetention;
	replicationVsTime[[i]] = residentReplication;
	
	timeVsTime[[i]] = SessionTime[] - timeRes;
	
	(* Print timing info *)
	(*Print["Resident " <> ToString[i-1] <> ", (retention = " <> ToString[residentRetention] <> 
		", replication = " <> ToString[residentReplication] <> "): " <> ToString[N[timeVsTime[[i]]]]];
		
	Print["    Time eco: " <> ToString[N[timeEco]] <> ", time min pop: " <> ToString[N[timeMinPop]] <> 
		", time avg pop: " <> ToString[N[timeAvgPop]] <> ", time fix: " <> ToString[N[timeFix]]]*)
	, (* comma here to signal the end of the for-loop *)
		
	(* End of the for-loop. We write the variable to iterate over (i) and the values to iterate over, 2 to 501 (we introduce 500 mutations) *)
	{i, 2, maxMutationSteps + 1}];

(* Return retention and replication values over time *)
Return[{retentionVsTime, replicationVsTime, timeVsTime, cyclingVsTime}]]


(* ::Subsection::Closed:: *)
(*Test the function that simulates evolution*)


testRun = Monitor[simulateEvolution[50, attackTradeOff -> False, conversionTradeOff -> True, maxMutations -> 500, initialRetention -> 0, initialReplication -> 0.1], i];


varyReplication = Monitor[Table[
	{eta, simulateEvolution[50, attackTradeOff -> False, conversionTradeOff -> True, maxMutations -> 500, initialRetention -> 0, initialReplication -> eta]}, 
	{eta, 0, 1, 1/10}], eta];
	


testRun[[1]]


backMuts = Pick[MapThread[{#1, #2}&, {Transpose[testRun[[1;;2, 1;;499]]], Transpose[testRun[[1;;2, 2;;500]]]}],
 Table[testRun[[1;;2, i]] == testRun[[1;;2, i + 2]], {i, 1, 499}]];


backMuts[[1]]


testSim = mutantAndResidentEcoDynamics[0.03, 0.01, 0.04, 0.01, 25, 
1, 1, 1,  1/100, 1/100, attackTradeOff -> True, conversionTradeOff -> True];
testSim2 = mutantAndResidentEcoDynamics[0.04, 0.01, 0.03, 0.01, 25, 
1, 1, 1,  1/100, 1/100, attackTradeOff -> True, conversionTradeOff -> True];


Column[{
Plot[{(h[t] + p[t]) /. testSim[[1]], (h2[t] + p2[t]) /. testSim[[1]]}, {t, 0, 4000}, 
	FrameLabel -> {"Time (days)", "Pop. size"}, Frame -> {True, True, False, False}, PlotStyle -> {Purple, Blue},
	LabelStyle -> {16, Black}, ImageSize -> 300,
	PlotLegends -> {"Resident (\[Phi] = 0.03, \[Eta] = 0.01)", "Mutant (\[Phi] = 0.04, \[Eta] = 0.01)"}],
	
Plot[{(h[t] + p[t]) /. testSim2[[1]], (h2[t] + p2[t]) /. testSim2[[1]]}, {t, 0, 4000}, 
	FrameLabel -> {"Time (days)", "Pop. size"}, Frame -> {True, True, False, False}, PlotStyle -> {Blue, Red},
	LabelStyle -> {16, Black}, ImageSize -> 300,
	PlotLegends -> {"Mutant/new resident (\[Phi] = 0.04, \[Eta] = 0.01)", "Back mutant (\[Phi] = 0.03, \[Eta] = 0.01)"}]
}]

(*Export["Mutant_and_backmutant_sample.png", %]*)


Plot[Log[(h2[t] + p2[t])] /. testSim2[[1]], {t, 0, 4000}]


Plot[Log[(h[t] + p[t])] /. testSim[[1]], {t, 0, 4000}]


Row[{Show[
ListLinePlot[Transpose[testRun][[All, 1;;2]], Frame -> True, FrameLabel -> {"Retention", "Replication"}, AspectRatio -> 1, PlotRange -> {{0, 1}, {0, 1}},
	LabelStyle -> {16, Black}, PlotStyle -> Directive[Thin, Gray], ImageSize -> 500],

ListPlot[
	MapThread[Style[#1, #2, PointSize[Medium]]&, {Transpose[testRun][[All, 1;;2]],
		Table[ColorData["CandyColors"][i/(Dimensions[testRun][[2]] - 1)], {i, 0, (Dimensions[testRun][[2]] - 1)}]}]]
],

PointLegend[{ColorData["CandyColors"][0], ColorData["CandyColors"][0.25], 
ColorData["CandyColors"][0.5], ColorData["CandyColors"][0.75],
ColorData["CandyColors"][1]}, {"Initial population", ".", ".", ".", "Final population"}, LabelStyle -> {16, Black}]}]

(*Export["sample_trajectory_fix_prob.png", %];*)


testRun[[4]]


Column[{
ListLinePlot[{testRun[[1]], testRun[[2]]}, PlotLegends -> {"Retention", "Replication"}, Frame -> {True, True, False, False}, 
	PlotLabel -> "Traits vs. time",
	FrameLabel -> {"Mutation", "Trait"}, ImageSize -> 300, LabelStyle -> {16, Black}],
	
ListLinePlot[testRun[[3]], PlotRange -> Full, PlotStyle -> Directive[Thin, Black],
	PlotLabel -> "Computational time",
	Frame -> {True, True, False, False}, FrameLabel -> {"Mutation", "Time (sec)"}, ImageSize -> 300, LabelStyle -> {16, Black}]
}]

(*Export["sample_time_trajectory_fix_prob.png", %]*)


Histogram[testRun[[3]], {"Log", 50}, "Probability",
PlotRange -> Full, Frame -> {True, True, False, False},
PlotLabel -> "Time to simulate 1 mutant's invasion",
FrameLabel -> {"Computational time (sec)", "Frequency"}, LabelStyle -> {16, Black}]

Export["sample_time_hist_fix_prob.png", %]





Table[ListLinePlot[{varyReplication[[i, 2, 1]], varyReplication[[i, 2, 2]]}, PlotLegends -> {"Retention", "Replication"}, Frame -> {True, True, False, False}, 
	PlotLabel -> "Initial replication = " <> ToString[N[varyReplication[[i, 1]]]],
	FrameLabel -> {"Mutation", "Trait"}, ImageSize -> 300, LabelStyle -> {16, Black}], {i, Length[varyReplication]}]


Table[ListLinePlot[Transpose[varyReplication[[i, 2, 1;;2]]], Frame -> {True, True, False, False}, 
	PlotLabel -> "Initial replication = " <> ToString[N[varyReplication[[i, 1]]]], PlotRange -> {{0, 1}, {0, 1}},
	FrameLabel -> {"Retention", "Replication"}, ImageSize -> 300, LabelStyle -> {16, Black}], {i, Length[varyReplication]}]


(* ::Subsection::Closed:: *)
(*Simulation ecological dynamics to check mutant fixation*)


(* ::Text:: *)
(*No attack trade-off, retention = 1, replication = 0.1, light = 50*)
(*Introduce mutants with slightly higher and lower replication probabilities*)


(* ::Subsubsection:: *)
(*Quick function to find ecological dynamics *)


Options[longEcoDynamics] =
Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[NDSolve]];


longEcoDynamics[residentRetention_ (* phi *), residentReplication_ (* eta *), mutantRetention_ (* phi2 *), mutantReplication_ (* eta2 *), surfaceLight_ (* iIn *), 
initialPrey_ (* W[0] *), initialResidentHeterotroph_ (* H[0] *), initialResidentPhototroph_ (* P[0] *),  
initialMutantHeterotroph_ (* H2[0] *), initialMutantPhototroph_ (* P2[0] *), 
time_ (* days to simulate dynamics over *), opts:OptionsPattern[]] :=

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
	
	(* Find the population sizes over time *)
	{w, h, p, h2, p2}, {t, 0, time},
	
	(* Optional input to adjust the differential equation solver *)
	FilterRules[{opts}, Options[NDSolve]]]; 
	
(* Clean up the formatting of ecoDynamics *)
ecoDynamics = Flatten[ecoDynamics, 1];

(* Check that extinction did not occur on the last time step *)
If[extinction == False,
 extinction = ((h[time] <= 10^-6) && (p[time] <= 10^-6) && (h2[time] <= 10^-6) && (p2[time] <= 10^-6)) /. ecoDynamics];

(* Return the ecological dynamics and whether extinction occured *)
Return[{ecoDynamics, extinction}]]


(* ::Subsubsection::Closed:: *)
(*Mutant with slightly higher replication probability*)


ecoEx1 = longEcoDynamics[1, 1/10, 1, 1/10 + 1/100, 50, (* traits and light *)
4, 99/100, 99/100, 1/100, 1/100, (* initial population sizes *)
50000 (* days to simulate dynamics over *), conversionTradeOff -> True, attackTradeOff -> False];


(* ::Subsubsection::Closed:: *)
(*Mutant with slightly lower replication probability*)


ecoEx2 = longEcoDynamics[1, 1/10, 1, 1/10 - 1/100, 50, (* traits and light *)
4, 99/100, 99/100, 1/100, 1/100, (* initial population sizes *)
50000 (* days to simulate dynamics over *), conversionTradeOff -> True, attackTradeOff -> False];


(* ::Subsubsection::Closed:: *)
(*Plot ecological dynamics*)


plotEx1 = Plot[{(h2[t] + p2[t]) / (h[t] + p[t] + h2[t] + p2[t]) /. ecoEx1[[1]]},
{t, 0, 50000},
PlotRange -> {{0, 50000}, {0, 1}},
Frame -> {True, True, False, False},
FrameLabel -> {"Time (days)", "Mutant fraction of population"},
PlotLabel -> "Mutant with higher replication (0.11)",
LabelStyle -> {"Arial", Black, 14}];


plotEx2 = Plot[{(h2[t] + p2[t]) / (h[t] + p[t] + h2[t] + p2[t]) /. ecoEx2[[1]]},
{t, 0, 50000},
PlotRange -> {{0, 50000}, {0, 0.05}},
Frame -> {True, True, False, False},
FrameLabel -> {"Time (days)", "Mutant fraction of population"},
PlotLabel -> "Mutant with lower replication (0.09)",
LabelStyle -> {"Arial", Black, 14}, PlotStyle -> Orange];


Labeled[GraphicsRow[{plotEx2, plotEx1}, ImageSize -> 800],
{Style["Mutant dynamics with resident retention = 1, replication = 0.1", FontFamily -> "Arial", FontColor -> Black, FontSize -> 18],
Style["Matches streamplot predictions that higher replication is more fit.\nNo attack tradeoff, light = 50 \[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s.", 
	Italic, FontFamily -> "Arial", FontColor -> Black, FontSize -> 16]},
{Top, Bottom}]
Export["mutant_fixation_retention_1_replication_0.1.pdf", %]

