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
Options[invasionDerivatives] = Options[simulateEvolution] = 
(* The options *)
Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[NDSolve]];


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


(* ::Subsection:: *)
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


(* ::Section:: *)
(*Invasion fitness*)


(* ::Subsection:: *)
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
(*Function to constrain mutations to biologically possible values*)


(* ::Text:: *)
(*Simple function to take a list of mutant retention and replication values and makes those < 0 equal to 0, and those > 1 equal to 1*)


constrainMutation = Map[Max[Min[#, 1], 0]&];


(* ::Subsection::Closed:: *)
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