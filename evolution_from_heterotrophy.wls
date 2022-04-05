#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Parameters*)


(* ::Text:: *)
(*Note on the parameter c_HP in the paper: In the code we have replaced c_HP with its reciprocal, gamma_C. This doesn't come up anywhere for the user, but just be aware that setting gamma_C = 0.7 in the functions below is equivalent to setting c_HP to 10/7. *)


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


(* ::Section:: *)
(*Resident and mutant dynamics*)


(* ::Text:: *)
(*Set the default options : both trade - offs turned off. Turn these on or off when calling residentEcoDynamics, not here!*)
(*The options for NDSolve are also included, in case the differential equation solving run by residentEcoDynamics has to be adjusted*)


Options[residentEcoDynamics] = Options[mutantAndResidentEcoDynamics] = Options[simulateEvolution] = Join[{attackTradeOff -> False, conversionTradeOff -> False}, Options[NDSolve]];


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


(* ::Subsection::Closed:: *)
(*Function to simulate evolution*)


(* ::Text:: *)
(*Simple function to take a list of mutant retention and replication values and makes those < 0 equal to 0, and those > 1 equal to 1*)


constrainMutation = Map[Max[Min[#, 1], 0]&];


(* ::Text:: *)
(*A function to simulate  mutations arising in a population of mixotrophs, and fixing if the mutants appear to be outcompeting the resident*)


simulateEvolution[surfaceLight_, opts:OptionsPattern[]] :=
Module[{
(* Variables we'll use inside this function *)
residentRetention, residentReplication, mutantRetention, mutantReplication, (* the resident and mutant's traits *)
residentEcology, mutantAndResidentEcology, (* the ecological dynamics of mutants and residents *)
extinction, (* will be True if the entire mixotroph population went extinct *)
averagePrey, averageHeterotroph, averagePhototroph, (* average prey, heterotrophs, and phototrophs when the resident is alone *)
averageResidentHeterotroph, averageResidentPhototroph, averageMutantHeterotroph, averageMutantPhototroph, (* average mixotroph population sizes when the mutant is present *)
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
Print["Mixotroph with (phi = 0, eta = 0) does not survive when surface light is " <> ToString[surfaceLight] <> " \[Mu]mol photons/m^2/s"]
Return[{retentionVsTime, replicationVsTime}]];

(* Find the "average" ecological state with just the resident (integrate each of the state variables over the last 1000 days and divide by 1000 to get the average *)
{averagePrey, averageHeterotroph, averagePhototroph} = Table[
	NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000), 
	{stateVar, {w[t], h[t], p[t]}}];
	
(* Set the initial values of the resident and mutant populations, so that the mutant is 1% of the resident population *)
initialPrey = averagePrey;
{initialResidentHeterotroph, initialResidentPhototroph} = 0.99 * {averageHeterotroph, averagePhototroph};
{initialMutantHeterotroph, initialMutantPhototroph} = 0.01 * {averageHeterotroph, averagePhototroph};



(* Introduce mutations and see whether the mutants can invade. (Do is a for loop. We'll be looping over i = 2:501, which comes at the very end of all this code, because Mathematica is weird about for loop syntax) *)
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
	{averageResidentHeterotroph, averageResidentPhototroph, averageMutantHeterotroph, averageMutantPhototroph} = Table[
		NIntegrate[stateVar/. mutantAndResidentEcology, {t, 3000, 4000}]/(4000 - 3000), 
		{stateVar, {h[t], p[t], h2[t], p2[t]}}];
		
		
	(* Calculate the proportion of the population that are mutants *)
	mutantProportion = (averageMutantHeterotroph + averageMutantPhototroph) / (averageResidentHeterotroph + averageResidentPhototroph + averageMutantHeterotroph + averageMutantPhototroph);

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
(*Run simulations*)


(* ::Subsection::Closed:: *)
(*Function for classifying simulation outcomes*)


classifyOutcome[{retention_, replication_}] :=
Which[(replication == "NA" && retention == "NA"), "Extinction", (* NAs indicate extinction *)
	replication == 0, (* If replication evolves to 0, the mixotroph is either a heterotroph (does not retain plastid) or a strict kleptoplast (always retains plastid) *)
	Which[retention == 0, "Strict heterotroph",
	retention == 1, "Obligately kleptoplastidic mixotroph",
	0 < retention < 1, "No label"],
	replication == 1, (* If replication evolves to 1, the mixotroph is either a phototroph (does not retain plastid), a kleptoplastidic phototroph (always retains plastid) or an occasional kleptoplastidic phototroph (sometimes retains plastid *)
	Which[retention == 0, "Strict phototroph",
	retention > 0, "Facultatively kleptoplastidic phototroph"],
	0 < replication < 1, "No label"] (* We don't have a name for any state where replication is intermediate *) 


(* ::Text:: *)
(*A list of the outcomes we observe (we'll still be able to tell if things didn't get one of these labels, though)*)


(* outcomeLabels = {"Heterotroph", "Strict Kleptoplast", "Phototroph", "Kleptoplastidic Phototroph", "Occasional Kleptoplastidic Phototroph", "Extinction"}; *)
outcomeLabels = {(*"Strict kleptoplast", *)"Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"};


(* ::Text:: *)
(*Colors for plotting*)


(* myColors = "CandyColors"; *)
myColors = Append[(ColorData["CandyColors"] /@ Subdivide[3])[[{(*1,*) 2, 4}]], White];


(* ::Subsection:: *)
(*Light levels to investigate*)


(* ::Text:: *)
(*Light levels to simulate evolution from heterotrophy*)


simLightLevels = {22, 23, 24, 25, 26, 27, 28, 29, 30};
higherLightLevels = {50, 100, 150, 200, 250, 300};


(* ::Text:: *)
(*We also wanted to investigate slightly lower light levels. It turns out the heterotroph can't survive at these light levels, so we're putting them in a separate list.*)
(*(To avoid wasting simulation time, we'll just confirm that the heterotroph can't survive at these light levels, and thus there can be no evolution from strict heterotrophy.)*)


lowerLightLevels = {20, 21};


(* ::Text:: *)
(*Finally, we'll put all the light levels together for plotting.*)


lightLevels = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 50, 100, 150, 200, 250, 300};


(* ::Text:: *)
(*A simple function to label outcomes of evolution according to the classification scheme in the paper*)


(* ::Subsection::Closed:: *)
(*No trade-offs*)


evolutionSimsNT = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 10], 10]&, lightLevels];


Export["evolutionSimsNT.wdx", evolutionSimsNT]


evolutionSimsNT = Import["evolutionSimsNT.wdx"];


(* ::Subsubsection::Closed:: *)
(*Plot the trait values vs. time (just to see what's going on)*)


GraphicsRow[{
ListLinePlot[evolutionSimsNT[[3, All, 1]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}],
ListLinePlot[evolutionSimsNT[[3, All, 2]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"}]
}, ImageSize -> 800]


ListPlot[Transpose[evolutionSimsNT[[3, 8]]], Frame -> True, FrameLabel -> {"Retention probability (phi)", "Replication probability (eta)"}, PlotRange -> {{0, 1}, {0, 1}},
PlotLabel -> "Example simulation"]


(* ::Subsubsection::Closed:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


outcomesNT = Map[classifyOutcome[#[[All, -1]]]&, evolutionSimsNT, {2}];
outcomeCountsNT = Outer[Count, outcomesNT, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsNT = outcomeCountsNT / Dimensions[evolutionSimsNT][[2]];


(* ::Subsubsection::Closed:: *)
(*Plot a bar chart showing how light affects the endpoint of evolution ("outcome")*)


noTradeoffsChart = BarChart[outcomeCountsNT, ChartLayout -> "Stacked", ChartLegends -> outcomeLabels, ChartStyle -> myColors, 
Frame -> {True, True, False, False},
(* Put the light levels on the x-axis *)
FrameTicks -> {Table[{i, lightLevels[[i]]}, {i, Length[lightLevels]}], Automatic},
FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "No trade-offs"]


(* ::Subsection:: *)
(*Conversion efficiency trade-off*)


(* ::Subsubsection:: *)
(*Simulate evolution (light > 21 micromol photons/m^2/s)*)


evolutionSimsCE = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True], 100]&, simLightLevels];


Export["evolutionSimsCE_10e-6_light22-30.wdx", evolutionSimsCE];


evolutionSimsCEmost = Import["evolutionSimsCE_10^-6_1-100.wdx"];
evolutionSimsCE29 = Import["evolutionSimsCE_10^-6_1-100_light29.wdx"];


evolutionSimsCE = Join[evolutionSimsCEmost[[1;;7]], evolutionSimsCE29, evolutionSimsCEmost[[{8}]]];
Export["evolutionSimsCE_10e-6_light22-30.wdx", evolutionSimsCE];


evolutionSimsCE = Import["evolutionSimsCE_10e-6_light22-30.wdx"];


(* ::Subsubsection:: *)
(*Simulation evolution (light > 30 micromol photons/m^2/s)*)


evolutionSimsCEhighLight = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True], 10]&, higherLightLevels];


Export["evolutionSimsCE_10e-6_light50_300.wdx", evolutionSimsCEhighLight];


evolutionSimsCEhighLight = Import["evolutionSimsCE_10e-6_light50_300.wdx"];


GraphicsGrid[Table[
{ListLinePlot[evolutionSimsCEhighLight[[i, All, 1]], AxesLabel -> {"Time (mutations)", "Retention"}, PlotLabel -> "Light = " <> ToString[higherLightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium],
ListLinePlot[evolutionSimsCEhighLight[[i, All, 2]], AxesLabel -> {"Time (mutations)", "Replication"}, PlotLabel -> "Light = " <> ToString[higherLightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium]},
{i, Length[higherLightLevels]}]]
Export["evolutionSimsCE_10e-6_light50_300.pdf", %]


(* ::Subsubsection::Closed:: *)
(*Plot the trait values vs. time (just to see what's going on)*)


GraphicsGrid[Table[{
ListLinePlot[evolutionSimsCE[[i, All, 1]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}, PlotLabel -> "Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day", PlotRange -> {All, {0, 1}}],
ListLinePlot[evolutionSimsCE[[i, All, 2]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"}, PlotLabel -> "Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day", PlotRange -> {All, {0, 1}}]}, 
{i, Length[simLightLevels]}], 
ImageSize -> 800]


(* ::Subsubsection:: *)
(*Simulate strict heterotroph at low light (light 20-21 micromol photons/m^2/s)*)


(* ::Text:: *)
(*Strict heterotroph never survives (extinction = True in all cases)*)


AllTrue[Table[residentEcoDynamics[0, 0, 20, 1, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False], 100][[All, 2]], #&]


AllTrue[Table[residentEcoDynamics[0, 0, 21, 1, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False], 100][[All, 2]], #&]


(* ::Text:: *)
(*Add these cases to the higher light simulated outcomes*)


evolutionSimsCE = Prepend[evolutionSimsCE, Table[{Table["NA", 1001], Table["NA", 1001]}, 100]];
evolutionSimsCE = Prepend[evolutionSimsCE, Table[{Table["NA", 1001], Table["NA", 1001]}, 100]];


(* ::Subsubsection:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


outcomesCE = Map[classifyOutcome[#[[All, -1]]]&, evolutionSimsCE, {2}];
outcomeCountsCE = Outer[Count, outcomesCE, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsCE = outcomeCountsCE / Dimensions[evolutionSimsCE][[2]];


(* ::Subsubsection:: *)
(*Plot a bar chart showing how light affects the endpoint of evolution ("outcome")*)


conversionTradeOffChart = BarChart[outcomeCountsCE, ChartLayout -> "Stacked", ChartLegends -> outcomeLabels, ChartStyle -> myColors, 
Frame -> {True, True, False, False},
(* Put the light levels on the x-axis *)
FrameTicks -> {Table[{i, lightLevels[[i]]}, {i, Length[lightLevels]}], Automatic},
FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "No trade-off"] (* This is the title because we're assuming everything has the CE trade-off, so we're just calling it normal *)


(* ::Subsubsection:: *)
(*Plot a line chart showing endpoint of evolution vs light*)


outcomesAndLightCE = MapThread[{#1, #2}&, {{lightLevels, lightLevels, lightLevels}, Transpose[outcomeCountsCE]}, 2];


linePlotCE = ListLinePlot[outcomesAndLightCE, PlotStyle -> {myColors[[1]], myColors[[2]], Black},
PlotMarkers -> Automatic,
PlotLegends -> Placed[LineLegend[{"Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"}], {0.65, 0.5}],
Frame -> True, FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "Simulated evolution from strict heterotrophy\n No tradeoff", 
ImageSize -> 450, LabelStyle -> {Black, 12, "Arial"}]


(* ::Subsection::Closed:: *)
(*Attack rate trade-off*)


evolutionSimsAR = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 10, attackTradeOff -> True], 10]&, lightLevels];


Export["evolutionSimsAR.wdx", evolutionSimsAR]


evolutionSimsAR = Import["evolutionSimsAR.wdx"];


(* ::Subsubsection::Closed:: *)
(*Plot the trait values vs. time (just to see what's going on)*)


GraphicsGrid[Table[{
ListLinePlot[evolutionSimsAR[[i, All, 1]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}, PlotLabel -> "Light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/day"],
ListLinePlot[evolutionSimsAR[[i, All, 2]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"}, PlotLabel -> "Light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/day"]}, 
{i, 3}], 
ImageSize -> 800]


(* ::Subsubsection::Closed:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


outcomesAR = Map[classifyOutcome[#[[All, -1]]]&, evolutionSimsAR, {2}];
outcomeCountsAR = Outer[Count, outcomesAR, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsAR = outcomeCountsAR / Dimensions[evolutionSimsAR][[2]];


(* ::Subsubsection::Closed:: *)
(*Plot a bar chart showing how light affects the endpoint of evolution ("outcome")*)


attackTradeOffChart = BarChart[outcomeCountsAR, ChartLayout -> "Stacked", ChartLegends -> outcomeLabels, ChartStyle -> myColors, 
Frame -> {True, True, False, False},
(* Put the light levels on the x-axis *)
FrameTicks -> {Table[{i, lightLevels[[i]]}, {i, Length[lightLevels]}], Automatic},
FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "Attack rate trade-off"]


(* ::Subsection:: *)
(*Conversion efficiency and attack rate trade-offs*)


(* ::Subsubsection:: *)
(*Simulate evolution (light > 21 micromol photons/m^2/s)*)


evolutionSimsCEAR = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True], 100]&, simLightLevels];


Export["evolutionSimsCEAR_10e-6_light22-30.wdx", evolutionSimsCEAR];


evolutionSimsCEARmost = Import["evolutionSimsCEAR_10e-6_sims1-100.wdx"];
evolutionSimsCEAR29 = Import["evolutionSimsCEAR_10e-6_sims1-100_light29.wdx"];


evolutionSimsCEAR = Join[evolutionSimsCEARmost[[1;;7]], evolutionSimsCEAR29, evolutionSimsCEARmost[[{8}]]];
Export["evolutionSimsCEAR_10e-6_light22-30.wdx", evolutionSimsCEAR];


evolutionSimsCEAR = Import["evolutionSimsCEAR_10e-6_light22-30.wdx"];


(* ::Subsubsection:: *)
(*Plot the trait values vs. time (just to see what's going on)*)


GraphicsGrid[Table[{
ListLinePlot[evolutionSimsCEAR[[i, All, 1]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}, PlotLabel -> "Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day", PlotRange -> {Automatic, {0, 1}}],
ListLinePlot[evolutionSimsCEAR[[i, All, 2]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"}, PlotLabel -> "Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day"]}, 
{i, Length[simLightLevels]}], 
ImageSize -> 800]


(* ::Subsubsection:: *)
(*Simulate strict heterotroph at low light (light 20-21 micromol photons/m^2/s)*)


(* ::Text:: *)
(*Strict heterotroph never survives (extinction = True in all cases)*)


AllTrue[Table[residentEcoDynamics[0, 0, 20, 1, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True], 100][[All, 2]], #&]


AllTrue[Table[residentEcoDynamics[0, 0, 21, 1, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True], 100][[All, 2]], #&]


(* ::Text:: *)
(*Add these cases to the higher light simulated outcomes*)


evolutionSimsCEAR = Prepend[evolutionSimsCEAR, Table[{Table["NA", 1001], Table["NA", 1001]}, 100]];
evolutionSimsCEAR = Prepend[evolutionSimsCEAR, Table[{Table["NA", 1001], Table["NA", 1001]}, 100]];


(* ::Subsubsection:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


outcomesCEAR = Map[classifyOutcome[#[[All, -1]]]&, evolutionSimsCEAR, {2}];
outcomeCountsCEAR = Outer[Count, outcomesCEAR, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsCEAR = outcomeCountsCEAR / Dimensions[evolutionSimsCEAR][[2]];


(* ::Subsubsection:: *)
(*Plot a bar chart showing how light affects the endpoint of evolution ("outcome")*)


conversionAndAttackTradeOffChart = BarChart[outcomeCountsCEAR, ChartLayout -> "Stacked", ChartLegends -> outcomeLabels, ChartStyle -> myColors, 
Frame -> {True, True, False, False},
(* Put the light levels on the x-axis *)
FrameTicks -> {Table[{i, lightLevels[[i]]}, {i, Length[lightLevels]}], Automatic},
FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "Attack rate trade-off"] (* This is the title because we're assuming everything has the CE trade-off, so we're just calling it normal *)


(* ::Subsubsection:: *)
(*Plot a line chart showing endpoint of evolution vs light*)


outcomesAndLightCEAR = MapThread[{#1, #2}&, {{lightLevels, lightLevels, lightLevels}, Transpose[outcomeCountsCEAR]}, 2];


linePlotCEAR = ListLinePlot[outcomesAndLightCEAR, PlotStyle -> {myColors[[1]], myColors[[2]], Black},
PlotMarkers -> Automatic,
PlotLegends -> {"Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"},
Frame -> True, FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "Simulated evolution from strict heterotrophy\n Attack rate tradeoff",
LabelStyle -> {Black, 12, "Arial"}, ImageSize -> 450]


(* ::Section:: *)
(*Figures: evolutionary simulations*)


(* ::Subsection:: *)
(*Plot all trade-off variations together*)


?Row


(* Add letter labels to plots *)
linePlots = MapThread[Show[#1, PlotLabel ->  #2 <> (PlotLabel /. Options[#1[[1]]])]&, {{linePlotCE, linePlotCEAR}, {"(A) ", "(B) "}}];

(* Remove legend from the second plot *)
linePlots = {linePlots[[1]], linePlots[[2, 1]]};

(* Put the plots in a row with a single, shared legend *)
Row[linePlots, Spacer[30]]

(* Export as pdf *)
Export["evo_from_heterotrophy.pdf", %];


(* ::Subsection:: *)
(*Plot all evolutionary trajectories together*)


(* ::Text:: *)
(*Conversion trade-off*)


trajectoryPlotsCE = Table[{
ListLinePlot[evolutionSimsCE[[i + Length[lowerLightLevels], All, 1]], Frame -> True, PlotRange -> {Automatic, {0, 1}}, PlotStyle -> Directive[Thin, Gray, Opacity[0.5]],
FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}, PlotRangePadding -> 0.025, ImageSize -> 350,
PlotLabel -> "(" <> ToUpperCase[Alphabet[][[1 + 2*(i - 1)]]] <> ") Retention: Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day"],

ListLinePlot[evolutionSimsCE[[i + Length[lowerLightLevels], All, 2]], Frame -> True, PlotRange -> {Automatic, {0, 1}}, PlotStyle -> Directive[Thin, Gray, Opacity[0.5]],
FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"},  PlotRangePadding -> 0.025, ImageSize -> 350,
PlotLabel -> "(" <> ToUpperCase[Alphabet[][[2*i]]] <> ") Replication: Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day"]}, 
{i, Length[simLightLevels]}];


trajectoryPlotsCE = GraphicsGrid[Transpose[Flatten[Partition[trajectoryPlotsCE, 2], {2, 3}]], Frame -> {None, All, {{1, 4}, {1, 2}} -> True}];


Export["trajectories_ce.pdf", trajectoryPlotsCE]


(* ::Text:: *)
(*Attack & conversion trade-offs*)


trajectoryPlotsCEAR = Table[{
ListLinePlot[evolutionSimsCEAR[[i + Length[lowerLightLevels], All, 1]], Frame -> True, PlotRange -> {Automatic, {0, 1}}, PlotStyle -> Directive[Thin, Gray, Opacity[0.5]],
FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}, PlotRangePadding -> 0.025, ImageSize -> 350,
PlotLabel -> "(" <> ToUpperCase[Alphabet[][[1 + 2*(i - 1)]]] <> ") Retention: Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day"],

ListLinePlot[evolutionSimsCEAR[[i + Length[lowerLightLevels], All, 2]], Frame -> True, PlotRange -> {Automatic, {0, 1}}, PlotStyle -> Directive[Thin, Gray, Opacity[0.5]],
FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"},  PlotRangePadding -> 0.025, ImageSize -> 350,
PlotLabel -> "(" <> ToUpperCase[Alphabet[][[2*i]]] <> ") Replication: Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day"]}, 
{i, Length[simLightLevels]}];


trajectoryPlotsCEAR = GraphicsGrid[Transpose[Flatten[Partition[trajectoryPlotsCEAR, 2], {2, 3}]], Frame -> {None, All, {{1, 4}, {1, 2}} -> True}];


Export["trajectories_cear.pdf", trajectoryPlotsCEAR]


(* ::Text:: *)
(*Return the label style options to their defaults*)


SetOptions[Labeled, LabelStyle -> {}];


(* ::Section:: *)
(*Figures: example ecology and evolution plots*)


(* ::Subsection:: *)
(*Colors and fonts*)


(* ::Text:: *)
(*Colors for plotting mixotroph population dynamics*)


(*colorH = ColorData["HTML"]["Pink"];
colorH2 = ColorData["HTML"]["DeepPink"];
colorP = Darker[ColorData["HTML"]["LightGreen"]];
colorP2 = Darker[ColorData["HTML"]["Green"]];*)

colorW = RGBColor["#77aa77"];
colorH = colorH2 = RGBColor["#fdc552"];
colorP = colorP2 = RGBColor["#4573ba"];
colorRetention = RGBColor["#a02b2b"];
colorReplication = RGBColor["#2c89a0"];


(* ::Text:: *)
(*Label font options (this saves so much trouble to specify here!)*)


SetOptions[Labeled, LabelStyle -> {FontColor -> Black, FontSize -> 12, FontFamily -> "Arial"}];


(* ::Subsection:: *)
(*Ecology plots*)


(* ::Subsubsection:: *)
(*Run the simulation*)


(* ::Text:: *)
(*First simulate a heterotroph alone (resident) to get the initialization for a simulation with mutants*)


{residentEcology, extinction} = residentEcoDynamics[0, 0, 25, 4, 1, 1, conversionTradeOff -> True, attackTradeOff -> True]


Plot[{w[t]/.residentEcology, h[t]/.residentEcology, p[t]/.residentEcology}, {t, 0, 300}, PlotRange -> All]


{avgPrey, avgHeterotroph, avgPhototroph} = Table[NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000), {stateVar, {w[t], h[t], p[t]}}];


(* ::Text:: *)
(*Simulate a resident heterotrophic population being invaded by a mutant that can sometimes replicate the plastid*)


{withMutantEcology, extinction} = mutantAndResidentEcoDynamics[0, 0, 1/100, 0, 25, avgPrey, (99/100) * avgHeterotroph, (99/100) * avgPhototroph,  (1/100) * avgHeterotroph, (1/100) * avgPhototroph, conversionTradeOff -> True, attackTradeOff -> True]


?GridLines


(* ::Subsubsection:: *)
(*Plot prey ecological dynamics*)


preyPlot = Plot[w[t-100] /. withMutantEcology, {t, 100, 1000}, PlotRange -> {{0, 1000}, {0, 20}}, PlotStyle -> colorW, Frame -> True, 
GridLines -> {{{100, Directive[Thick, Black, Dashed]}}, None},
FrameLabel -> {{"Prey (cells/mL)", None}, {"Time (days)", None}}, ImageSize -> 350, ImagePadding -> {{50, 50}, {50, 10}}, PlotRangePadding -> 0];
preyPlot = Show[preyPlot, Plot[w[t+3900] /. residentEcology, {t, 0, 100}, PlotStyle -> colorW]]


(* ::Subsubsection:: *)
(*Plot mixotroph ecological dynamics*)


(* ::Text:: *)
(*Plot legend*)


mixotrophPlotLegend = LineLegend[{colorH, Directive[colorH2, Dashed], colorP, Directive[colorP2, Dashed]}, {"Resident heterotroph", "Mutant heterotroph", "Resident phototroph", "Mutant phototroph"}];


(* ::Text:: *)
(*This code (commented out) plots the resident and mutant with different y-axes. Not using it because we decided to show the whole time series of invasion, so the mutant population grows to approximately the resident's original size*)


(*residentPlot = Plot[{h[t] /. withMutantEcology, p[t] /. withMutantEcology}, {t, 0, 1000}, PlotRange \[Rule] {{0, 1000}, {0, 3}}, PlotStyle \[Rule] {Directive[Thick, colorH], Directive[Thick, colorP]}];
mutantPlot = Plot[{h2[t] /. withMutantEcology, p2[t] /. withMutantEcology}, {t, 0, 1000}, PlotRange \[Rule] {{0, 1000}, {0, 3}}, PlotStyle \[Rule] {Directive[Thin, Dashed, colorH2], Directive[Thin, Dashed, colorP2]}];

residentRange = (PlotRange /. AbsoluteOptions[residentPlot, PlotRange])[[2]];
mutantRange = (PlotRange /. AbsoluteOptions[mutantPlot, PlotRange])[[2]];

residentTicks = N[FindDivisions[residentRange, 5]];
mutantTicks = N[FindDivisions[mutantRange, 5]];
mutantTicks = {Rescale[#, mutantRange, residentRange], ToString[#]}& /@ mutantTicks;

Show[residentPlot, Graphics[GeometricTransformation[mutantPlot[[1]], RescalingTransform[{{0, 1}, mutantRange}, {{0, 1}, residentRange}]]],
Frame \[Rule] True, FrameTicks \[Rule] {{residentTicks, mutantTicks}, {Automatic, Automatic}},
FrameLabel \[Rule] {{"Resident biomass (cells/mL)", "Mutant biomass (cells/mL)"}, {None, None}}]*)


(* ::Text:: *)
(*This code plots heterotrophs (mutant and resident) and phototrophs together, with different y-axes. Useful because there aren't that many phototrophs for these particular trait values.*)


heterotrophPlot = Plot[{h[t - 100] /. withMutantEcology // Chop, h2[t - 100] /. withMutantEcology // Chop}, {t, 100, 1000}, 
GridLines -> {{{100, Directive[Thick, Black, Dashed]}}, None},
PlotRange -> {{0, 1000}, {0, 0.5}}, PlotStyle -> {Directive[Thick, colorH], Directive[Thick, Dashed, colorH2]}];
phototrophPlot = Plot[{p[t - 100] /. withMutantEcology // Chop, p2[t - 100] /. withMutantEcology // Chop}, {t, 100, 1000}, PlotRange -> {{0, 1000}, {0, 0.05}}, PlotStyle -> {Directive[Thick, colorP], Directive[Thick, Dashed, colorP2]}];
heterotrophPlot = Show[heterotrophPlot, Plot[h[t + 3900] /. residentEcology // Chop, {t, 0, 100}, PlotRange -> {{0, 1000}, {0, 0.5}}, PlotStyle -> Directive[Thick, colorH]]];
phototrophPlot = Show[phototrophPlot, Plot[p[t + 3900] /. residentEcology // Chop, {t, 0, 100}, PlotRange -> {{0, 1000}, {0, 0.05}}, PlotStyle -> Directive[Thick, colorP]]];

heterotrophRange = (PlotRange /. AbsoluteOptions[heterotrophPlot, PlotRange])[[2]];
phototrophRange = (PlotRange /. AbsoluteOptions[phototrophPlot, PlotRange])[[2]];

heterotrophTicks = N[FindDivisions[heterotrophRange, 5]];
phototrophTicks = N[FindDivisions[phototrophRange, 5]];
phototrophTicks = {Rescale[#, phototrophRange, heterotrophRange], ToString[#]}& /@ phototrophTicks;

mixotrophPlot = Show[heterotrophPlot, Graphics[GeometricTransformation[phototrophPlot[[1]], RescalingTransform[{{0, 1}, phototrophRange}, {{0, 1}, heterotrophRange}]]],
Frame -> True, FrameTicks -> {{heterotrophTicks, phototrophTicks}, {Automatic, Automatic}},
FrameLabel -> {{"Heterotroph biomass (cells/mL)", "Phototroph biomass (cells/mL)"}, {"Time (days)", None}}, ImageSize -> 350, ImagePadding -> {{50, 50}, {50, 10}}, PlotRangePadding -> 0];


(* ::Subsubsection:: *)
(*Plot the proportion of the population that is the mutant*)


mutantProportionPlot = Plot[((h2[t - 100] + p2[t - 100])/ (h[t - 100] + p[t - 100] + h2[t - 100] + p2[t - 100]))/. withMutantEcology, {t, 100, 1000}, PlotRange -> {{0, 1000}, {0, 1}}, PlotStyle -> Black, Frame -> True,
GridLines -> {{{100, Directive[Thick, Black, Dashed]}}, None},
FrameLabel -> {"Time (days)", "Mutant proportion of population"}, ImageSize -> 350, ImagePadding -> {{50, 50}, {50, 10}}, PlotRangePadding -> 0];
mutantProportionPlot = Show[mutantProportionPlot, Plot[0, {t, 0, 100}, PlotStyle -> Black]];


(* ::Subsubsection:: *)
(*Combine plots into a column*)


ecologyPlots = Column[{Labeled[preyPlot, "(A) Prey ecology", Top, Alignment -> {Left, Baseline}],
Labeled[mixotrophPlot, "(B) Mixotroph ecology", Top, Alignment -> Left],
Labeled[mutantProportionPlot, "(C) Mutant spread", Top, Alignment -> Left]}]


(* ::Subsection:: *)
(*Evolution plots*)


(* ::Text:: *)
(*For these plots, we'll use the evolution simulations we already did earlier in this file (load them in here in case we did them last time and don't want to re-run)*)


evolutionSimsCEAR = Import["evolutionSimsCEAR_10e-6_light22-30.wdx"];


simLightLevels[[4]]


(* ::Text:: *)
(*Plot the plastid retention probability over evolutionary time*)


retentionPlot = ListLinePlot[evolutionSimsCEAR[[4, 1, 1]], Frame -> True, FrameLabel -> {"Time (mutations that have arisen)", "Retention"}, ImageSize -> 350, 
ImagePadding -> {{50, 50}, {50, 10}}, PlotStyle -> colorRetention, PlotRange -> {All, {0, 1}}, PlotRangePadding -> 0];


(* ::Text:: *)
(*Plot the plastid replication probability over evolutionary time*)


replicationPlot = ListLinePlot[evolutionSimsCEAR[[4, 1, 2]], Frame -> True, FrameLabel -> {"Time (mutations that have arisen)", "Replication"}, ImageSize -> 350, 
ImagePadding -> {{50, 50}, {50, 10}}, PlotStyle -> colorReplication, PlotRange -> {All, {0, 1}}, PlotRangePadding -> 0];


(* ::Text:: *)
(*Plot the evolutionary trajectory in trait space so we can see what all the different (retention, replication) combinations were*)


traitSpacePlot = ListLinePlot[Transpose[evolutionSimsCEAR[[4, 1]]], Frame -> True, FrameLabel -> {"Retention", "Replication"}, ImageSize -> 350, 
ImagePadding -> {{50, 50}, {50, 10}}, PlotStyle -> Black, PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0];


(* ::Subsubsection:: *)
(*Combine plots into a column*)


evolutionPlots = Column[{Labeled[retentionPlot, "(D) Retention evolution", Top, Alignment -> {Left, Baseline}],
Labeled[replicationPlot, "(E) Replication evolution", Top, Alignment -> Left],
Labeled[traitSpacePlot, "(F) Evolutionary trajectory in trait space", Top, Alignment -> Left]}]


(* ::Subsection:: *)
(*Combine ecology and evolution plots into the figure*)


Row[{mixotrophPlotLegend, ecologyPlots, evolutionPlots}]
Export["example_timeseries.pdf", %]



