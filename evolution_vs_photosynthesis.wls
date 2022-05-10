#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
(*Parameters*)


(* ::Text:: *)
(*Note on the parameter c_HP in the paper: In the code we have replaced c_HP with its reciprocal, gamma_C. This doesn't come up anywhere for the user, but just be aware that setting gamma_C = 0.7 in the functions below is equivalent to setting c_HP to 10/7. *)


pW = 3;
(*gammaP = 2/3;*)(*2/3;*)(*2.9/3;*)
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
(*iIn = 50;*) (* Incident light *)


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
(*Note: for this script only, we added an additional option to specify the photosynthetic rate individually for each simulation.*)
(*The option is photoEfficiency and represents gammaW, or the maximum photosynthetic rate of the phototroph relative to the maximum photosynthetic rate of the prey.*)
(*The default value is photoEfficiency = pP/pW = 2/3.*)


(* ::Subsection::Closed:: *)
(*Set default options*)


(* ::Text:: *)
(*Set the default options : both trade - offs turned off. Turn these on or off when calling residentEcoDynamics, not here!*)
(*The options for NDSolve are also included, in case the differential equation solving run by residentEcoDynamics has to be adjusted*)
(*We'll also add an extra option for the photosynthetic rate of the phototroph relative to the photosynthetic rate of the prey (photoEfficiency = gammaP, defaults to 2/3)*)


Options[residentEcoDynamics] = Options[mutantAndResidentEcoDynamics] = Options[simulateEvolution] = 
	Join[{attackTradeOff -> False, conversionTradeOff -> False, photoEffiency -> 2/3}, Options[NDSolve]];


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
gammaP = OptionValue[photoEffiency], (* The photosynthetic rate of phototroph = gammaP * pW, where pW is prey's photosynthetic rate *)

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
gammaP = OptionValue[photoEffiency], (* The photosynthetic rate of phototroph = gammaP * pW, where pW is prey's photosynthetic rate *)

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


(* ::Text:: *)
(*A simple function to label outcomes of evolution according to the classification scheme in the paper*)


classifyOutcome[{retention_, replication_}] :=
Which[(replication == "NA" && retention == "NA"), "Extinction", (* NAs indicate extinction *)
	replication == 0, (* If replication evolves to 0, the mixotroph is either a heterotroph (does not retain plastid) or a strict kleptoplast (always retains plastid) *)
		Which[retention == 0, "Strict heterotroph",
		retention == 1, "Obligately kleptoplastidic mixotroph",
		0 < retention < 1, "No label"],
	replication == 1, (* If replication evolves to 1, the mixotroph is either a phototroph (does not retain plastid), a kleptoplastidic phototroph (always retains plastid) or an occasional kleptoplastidic phototroph (sometimes retains plastid *)
		Which[retention == 0, "Strict phototroph",
		retention > 0, "Facultatively kleptoplastidic phototroph"],
	0 < replication < 1, (* If replication is intermediate, the mixotroph is either a heterotroph with some replication ability (does not retain plastid) or unlabeled *)
		Which[retention == 0, "Heterotroph with latent replication capacity",
		retention > 0, "No label"]] (* We don't have a name for any state where replication is intermediate and retention is > 0 *) 


(* ::Text:: *)
(*A list of the outcomes we observe (we'll still be able to tell if things didn't get one of these labels, though)*)


(* outcomeLabels = {"Heterotroph", "Strict Kleptoplast", "Phototroph", "Kleptoplastidic Phototroph", "Occasional Kleptoplastidic Phototroph", "Extinction"}; *)
outcomeLabels = {"Strict heterotroph", "Heterotroph with latent replication capacity", "Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"};


(* ::Text:: *)
(*Colors for plotting*)


(* myColors = "CandyColors"; *)
myColors = Append[(ColorData["CandyColors"] /@ Subdivide[3])[[{1, 3, 2, 4}]], White];


(* ::Subsection:: *)
(*Photosynthesis levels to investigate*)


(* ::Text:: *)
(*Maximum photosynthetic rates to use in simulations*)


pPlevels = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2}
gammaPlevels = pPlevels/pW (* Gamma P is the photosynthetic rate relative to the prey's photosynthetic rate *)


(* ::Subsection:: *)
(*Conversion efficiency trade-off*)


(* ::Subsubsection::Closed:: *)
(*Simulate evolution*)


evolutionSimsCE = Map[ParallelTable[simulateEvolution[25, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, photoEffiency -> #], 10]&, gammaPlevels];


Export["evolutionSimsCE_pP_light25.wdx", evolutionSimsCE];


(* ::Subsubsection:: *)
(*Plot the trait values vs. time (just to see what's going on)*)


GraphicsGrid[Table[{
ListLinePlot[evolutionSimsCE[[i, All, 1]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}, 
	PlotLabel -> "pP = " <> ToString[pPlevels[[i]]] <> " \[Mu]mol photons/m^2/day", PlotRange -> {{0, 1000}, {0, 1}}, PlotRangePadding -> 0.05],
ListLinePlot[evolutionSimsCE[[i, All, 2]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"}, 
PlotLabel -> "pP = " <> ToString[pPlevels[[i]]] <> " \[Mu]mol photons/m^2/day", PlotRange -> {{0, 1000}, {0, 1}}, PlotRangePadding -> 0.05]},
{i, Length[gammaPlevels]}], 
ImageSize -> 800]


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
(* Put the phototroph's maximum photosynthetic rate on the x-axis *)
FrameTicks -> {Table[{i, pPlevels[[i]]}, {i, Length[pPlevels]}], Automatic},
FrameLabel -> {"pP", "Fraction of simulations"},
PlotLabel -> "No trade-off"] (* This is the title because we're assuming everything has the CE trade-off, so we're just calling it normal *)


(* ::Subsubsection:: *)
(*Plot a line chart showing endpoint of evolution vs light*)


linePlotCE = ListLinePlot[outcomesAndpPCE[[{3, 4, 1, 2}]], PlotStyle -> myColors[[{3, 4, 1, 2}]], (* This ordering is just to match the main text order *)
PlotMarkers -> Graphics`PlotMarkers[][[{1, 2, 4, 5}]], (* Skip a symbol that represented extinction in the main text plot *)
PlotLegends -> LineLegend[outcomeLabels[[{3, 4, 1, 2}]]], (* Order the names to match the order the outcomes were plotted in *)
Frame -> True, 
FrameLabel -> {"Maximum photosynthetic rate of phototroph (" <> ToString[Superscript["day", "-1"], StandardForm] <> ")", "Fraction of simulations"},
PlotLabel -> "Endpoint of evolution from heterotrophy", 
ImageSize -> 450, LabelStyle -> {Black, 12, "Arial"}]

Export["evolution_from_heteroptrophy_vs_pP_light25_ce.pdf", %]



