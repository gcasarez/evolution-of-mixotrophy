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


(* ::Section::Closed:: *)
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


Options[residentEcoDynamics] = Options[rareMutantEcoDynamics] =  
	Join[{attackTradeOff -> False, conversionTradeOff -> False, photoEffiency -> 2/3}, Options[NDSolve]];


Options[residentEcoDynamics] = Append[Options[residentEcoDynamics], maxTime -> 4000];


Options[simulateEvolution] = Join[{attackTradeOff -> False, conversionTradeOff -> False,  photoEffiency -> 2/3,
	maxMutations -> 500, initialRetention -> 0, initialReplication -> 0}, 
	Options[NDSolve], Options[NMinimize], Options[NMaximize], Options[FindRoot]];


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
(*Function to simulate evolution*)


(* ::Text:: *)
(*Simple function to take a list of mutant retention and replication values and makes those < 0 equal to 0, and those > 1 equal to 1*)


(* ::Subsubsection::Closed:: *)
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
gammaP = OptionValue[photoEffiency], (* The photosynthetic rate of phototroph = gammaP * pW, where pW is prey's photosynthetic rate *)

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


(* ::Subsubsection::Closed:: *)
(*Helper function that estimate the period of an interpolating function*)


(* ::Text:: *)
(*This function will be used later to calculate the period of resident population cycling.*)
(*We need the period to more accurately calculate resident growth rate, so that we can measure mutant population changes over a whole number of cycles.*)


Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];
ParallelNeeds["DifferentialEquations`InterpolatingFunctionAnatomy`"];


(* ::Text:: *)
(*Error messages for the function. It requires an interpolating function with*)
(* - a 1-dimensional domain*)
(* - at least two peaks, so the period can be estimated from the distance between the peaks*)


interpolatingFunctionPeriod::DomainError = "Domain is not 1-dimensional. Dimensions of domain: `1`";
interpolatingFunctionPeriod::NotPeriodic = "Peaks detected: `1`. Check that function is periodic and include a larger domain if so.";
interpolatingFunctionPeriod::FewPeaks = "Cannot discard first 3 estimates of period, because only `1` estimates can be obtained.";


(* ::Text:: *)
(*This function takes an interpolating function as input and estimates the period of its cycles based on the median distance between successive peaks.*)
(*If the optional input ignoreEarlyPeaks is true (default value is true), then the the distances between the first 4 peaks are ignored when calculating the median period.*)


interpolatingFunctionPeriod[ifun_, ignoreEarlyPeaks_:True] :=
	Module[{peakLocations, peakSeparation}, 
	
	(* Check that the interpolating function domain is 1-dimensional *)
	If[Dimensions[InterpolatingFunctionDomain[ifun]] != {1, 2},
		Message[interpolatingFunctionPeriod::DomainError, Dimensions[InterpolatingFunctionDomain[ifun]]];
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


(* ::Text:: *)
(*This function returns a list of an interpolating function's peaks*)


interpolatingFunctionPeaks::DomainError = "Domain is not 1-dimensional. Dimensions of domain: `1`. If you expected a 1-dimensional domain, check that the Differential Equations package InterpolatingFunctionAnatomy is loaded.";
interpolatingFunctionPeaks::LacksPeaks = "Peaks detected: `1`. Check that function is periodic and include a larger domain if so.";
interpolatingFunctionPeaks::PeakSpacing = "Variable peak spacing: std/mean = `1`. Check that function is periodic.";


interpolatingFunctionPeaks[ifun_] :=
	Module[{peakLocations, peakSeparation}, 
	
	(* Check that the interpolating function domain is 1-dimensional *)
	If[Dimensions[InterpolatingFunctionDomain[ifun]] != {1, 2},
		Message[interpolatingFunctionPeaks::DomainError, Dimensions[InterpolatingFunctionDomain[ifun]]];
		Return[]];
		
	(* Find the locations of the peaks of the interpolating function *)
	peakLocations = Pick[Flatten[InterpolatingFunctionGrid[ifun]], PeakDetect[InterpolatingFunctionValuesOnGrid[ifun]], 1];
	
	(* Check that the function appears to be periodic: has multiple peaks *)
	If[Length[peakLocations] <= 1,
		Message[interpolatingFunctionPeaks::LacksPeaks, Length[peakLocations]]];
	
	(* Check that the function appears to be periodic: peaks are evenly spaced *)
	peakSeparation = Differences[peakLocations];
	If[StandardDeviation[peakSeparation] / Mean[peakSeparation] > 2/100,
		(* Just give a warning if this check is not met *)
		Message[interpolatingFunctionPeaks::PeakSpacing, StandardDeviation[peakSeparation] / Mean[peakSeparation]]];
		
	(* Return peak locations *)
	Return[peakLocations]]


(* ::Subsubsection::Closed:: *)
(*Function to simulate evolution *)


(* ::Text:: *)
(*Takes the surface light as input and simulates evolution from heterotrophy (default) or from another set of initial retention and replication values, if the options initialRetention and initialReplication are given.*)
(*Use attackTradeOff -> True and conversionTradeOff -> True to turn on the attack tradeoff and to set a conversion factor of cHP between heterotrophic and phototrophic biomass. (cHP set elsewhere with the other parameters.)*)


simulateEvolution::TooFewCycles = "Population cycles after day 4000: `1`. Cannot estimate mutant growth after 20 population cycles. Using constant population approximation instead.";


simulateEvolution[surfaceLight_, opts:OptionsPattern[]] :=
Module[{
(* Variables we'll use inside this function *)
residentRetention, residentReplication, (* the resident's traits *)
residentEcology, (* the ecological dynamics with just the resident present *)
extinction, (* will be True if the entire resident population goes extinct *)
minResident, maxResident, (* minimum and maximum resident population size *)
totalResident, (* total resident population size over time (heterotroph + phototroph *)
peaks, (* locations of peaks of resident population cycling *)
growthEstimateDays, (* days to use to estimate change in mutant population over time *)
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
		Evaluate[FilterRules[{opts, maxTime -> 6000}, Options[residentEcoDynamics]]]];
		
	(* If extinction comes back false, check that the resident population didn't drop to <= 10^-6 (we also consider that extinction) *)
	If[!(extinction),
		minResident = NMinimize[{(h[t] + p[t])/. residentEcology, 3000 <= t <= 6000}, t, Evaluate[FilterRules[{opts}, Options[NMinimize]]]][[1]];
		extinction = (minResident <= 10^-6)];
		
	(* If the resident cannot survive, stop! *)
	If[extinction,
	Print["Resident with (retention = " <> ToString[residentRetention] <> ", replication = " <> ToString[residentReplication] <> 
		") went extinct. No further mutations introduced."];
	Break[]];

	(* Find the maximum resident population size. We will use this to see if the population is cycling *)
	maxResident = NMaximize[{(h[t] + p[t])/. residentEcology, 3000 <= t <= 6000}, t, Evaluate[FilterRules[{opts}, Options[NMaximize]]]][[1]];
	
	(* Check for cycling. If it is present, set growthEstimateDays so that we calculate population sizes at the same point in the cycle *)
	If[minResident <= 0.99 * maxResident, 
		(* The difference between minimum and maximum population sizes suggests cycling. *)
			cyclingVsTime[[i - 1]] = True;
			
			(* Get the mixotroph population *)
			totalResident =  Interpolation[Transpose[{InterpolatingFunctionGrid[ (h /. residentEcology)] (* time points *), 
				InterpolatingFunctionValuesOnGrid[(h /. residentEcology)] + InterpolatingFunctionValuesOnGrid[(p /. residentEcology)] (* pop. sizes *)}],
				InterpolationOrder -> 1];
				
			(* Get the location of resident population peaks that occur after mutant invasion *)
			peaks = interpolatingFunctionPeaks[totalResident];
			peaks = Select[peaks, # >= 4000&];
			
			(* If there are enough peaks, store the 16th peak and the 20th peak, so growth rate can be estimated over 5 population cycles *)
			If[Length[peaks] >= 20,
				growthEstimateDays = {peaks[[16]], peaks[[20]]},
				
				(* If there are very few peaks, probably something is wrong, so print a warning and estimate growth rate over days 4750 - 5000 *)
				Message[simulateEvolution::TooFewCycles, Length[peaks]];
				cyclingVsTime[[i - 1]] = "Ignored";
				growthEstimateDays = {4750, 5000}
			],
				
		(* No cycling in the resident population. Estimate growth rate over days 4750 - 5000. *)
			cyclingVsTime[[i - 1]] = False;
			growthEstimateDays = {4750, 5000}
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
				4000 (* start time for finding eco. dynamics *), 6000 (* end time *), Evaluate[FilterRules[{opts}, Options[rareMutantEcoDynamics]]]];
				
			(* If extinction comes back false, check that mutant population is positive and non-tiny (greater than 10^-10)  *)
			If[!(mutantExtinction),
				minMutant = NMinimize[{(h2[t] + p2[t])/. mutantDynamics, 4000 <= t <= 6000}, t, Evaluate[FilterRules[{opts}, Options[NMinimize]]]][[1]];
				mutantExtinction = (minMutant <= 10^-10)];

			If[mutantExtinction,
				(* If the mutant went extinct, the probability that it fixes is 0 *)
					prFix = 0,
				(* If the mutant didn't go extinct... *)
					(*  Use mutant's growth rate to approximate the selection coefficient *)
					selectionCoefficient = (( (h2[growthEstimateDays[[2]]] + p2[growthEstimateDays[[2]]]) / (h2[growthEstimateDays[[1]]] + p2[growthEstimateDays[[1]]]) /. mutantDynamics) - 1) / (growthEstimateDays[[2]] - growthEstimateDays[[1]]);
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


(* ::Section::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Photosynthesis levels to investigate*)


(* ::Text:: *)
(*Maximum photosynthetic rates to use in simulations*)


pPlevels = {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2}
gammaPlevels = pPlevels/pW (* Gamma P is the photosynthetic rate relative to the prey's photosynthetic rate *)


(* ::Subsection::Closed:: *)
(*Simulate evolution (Conversion efficiency trade-off)*)


evolutionSimsCE = ParallelMap[simulateEvolution[25, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False, photoEffiency -> #]&,
	Flatten[Table[gammaP, {gammaP, gammaPlevels}, 10]]];


evolutionSimsCE = MapThread[{#1, #2}&, {Flatten[Table[gammaP, {gammaP, gammaPlevels}, 10]], evolutionSimsCE}];


Export["evolutionSimsCE_pP.wdx", evolutionSimsCE];


(* ::Text:: *)
(*Sort simulations by photosynthesis level*)


gEvoSimsCE = GatherBy[SortBy[evolutionSimsCE, #[[1]]&], #[[1]]&];
simEndpointsCE = gEvoSimsCE[[All, All, 2, 1;;2, -1]];


(* ::Section::Closed:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


(* ::Subsection::Closed:: *)
(*Initial classification*)


outcomesCE = Map[classifyOutcome, simEndpointsCE, {2}];
outcomeCountsCE = Outer[Count, outcomesCE, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsCE = outcomeCountsCE / Dimensions[outcomesCE][[2]];


Map[Total, outcomeCountsCE]


(* ::Subsection::Closed:: *)
(*Cases of evolution to extinction*)


extinctionSimsCE = Pick[gEvoSimsCE, outcomesCE, "Extinction"];
extinctionSimsCE = Select[extinctionSimsCE, UnsameQ[#, {}] &];


(* ::Subsubsection::Closed:: *)
(*Find the resident that went extinct*)


extinctResidents = Map[FirstPosition["NA"], extinctionSimsCE[[All, All, 2, 1]], {2}] - 1;


(* ::Text:: *)
(*Resident retention and replication*)


MapThread[#1[[2, 1;;2, #2[[1]]]]&, {extinctionSimsCE, extinctResidents}, 2]


(* ::Text:: *)
(*All of these residents are very nearly facultatively kleptoplastidic phototrophs, but not quite*)


(* ::Subsubsection::Closed:: *)
(*Plot the trajectories*)


Map[ListLinePlot[Transpose[#[[2, 1;;2]]],
	PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.02,
	Frame -> True, FrameLabel -> {"Retention", "Replication"}]&, 
extinctionSimsCE, {2}]


(* ::Subsubsection::Closed:: *)
(*Let's call these simulations "extinction on trajectory toward FKP"*)


outcomesCE = Map[StringReplace["Extinction" -> "Extinction on trajectory toward FKP"], outcomesCE, {2}];


outcomeLabels = StringReplace[outcomeLabels, "Extinction" -> "Extinction on trajectory toward FKP"];


(* ::Subsection::Closed:: *)
(*Collapse purely heterotrophic phenotypes into a single category*)


(* ::Text:: *)
(*Rename heterotrophic phenotypes to "heterotroph"*)


outcomesCE = Map[StringReplace["Strict heterotroph" -> "Heterotroph"], outcomesCE, {2}];
outcomesCE = Map[StringReplace["Heterotroph with latent replication capacity" -> "Heterotroph"], outcomesCE, {2}];


outcomeLabels = StringReplace[outcomeLabels, "Strict heterotroph" -> "Heterotroph"];
outcomeLabels = StringReplace[outcomeLabels, "Heterotroph with latent replication capacity" -> "Heterotroph"];


(* ::Text:: *)
(*Find their positions in the list of outcomes*)


hPositions = Flatten[Position[outcomeLabels, "Heterotroph"]];


(* ::Text:: *)
(*Put heterotrophic outcome counts in a single category*)


outcomeCountsCE[[All, hPositions[[1]]]] = Map[Total, outcomeCountsCE[[All, hPositions]]];


outcomeCountsCE = Map[Delete[#, hPositions[[2]]]&, outcomeCountsCE];


(* ::Text:: *)
(*Collapse heterotrohic outcome labels into a single category*)


outcomeLabels = Delete[outcomeLabels, hPositions[[2]]];


(* ::Section::Closed:: *)
(*Plot a line chart showing endpoint of evolution vs light*)


(* ::Subsection::Closed:: *)
(*Make a set of lists showing how each outcome changes in frequency as the photosynthetic rate changes.*)


(* ::Text:: *)
(*E.g. the first list will be a set of {pP, strict heterotroph outcome} pairs*)


outcomesAndpPCE = MapThread[{#1, #2}&, {Table[pPlevels, Dimensions[outcomeCountsCE][[2]]], Transpose[outcomeCountsCE]}, 2];


(* ::Subsection::Closed:: *)
(*Plot markers and legend*)


myMarkers = Graphics`PlotMarkers[][[{1, 2, 7, 4}]]; (* Start at 4 to skip a symbol that represented extinction in the main text plot *)
myMarkers[[All, 2]] = {10, 15, 15, 12};


myLegend = PointLegend[myColors[[{3, 4, 4, 1}]], 
{outcomeLabels[[2]], outcomeLabels[[3]], outcomeLabels[[4]], outcomeLabels[[1]]},
LegendMarkers -> myMarkers,
LabelStyle -> {Black, "Arial", 12}]


(* ::Subsection::Closed:: *)
(*Plot everything but facultatively kleptoplastidic phototrophy and extinction*)


(* ::Text:: *)
(*We'll add that one later because the extinction at pP = 1.5 means I want to make the point look different there*)


linePlotCE = ListLinePlot[outcomesAndpPCE[[{2, 1}]], PlotStyle -> myColors[[{3, 1}]], (* This ordering is just to match the main text order *)
PlotMarkers -> myMarkers[[{1, 4}]], 
Frame -> True, 
FrameLabel -> {"Maximum photosynthetic rate of phototroph (" <> ToString[Superscript["day", "-1"], StandardForm] <> ")", "Fraction of simulations"},
PlotLabel -> "Endpoint of evolution from heterotrophy", 
ImageSize -> 5 * 72, LabelStyle -> {Black, 12, "Arial"}]


(* ::Subsection::Closed:: *)
(*Make another list plot with just facultatively kleptoplastidic phototrophy and extinction on the path to FKP*)


fkpOutcome = ListLinePlot[{
	outcomesAndpPCE[[3, 1;;6]],
	Join[outcomesAndpPCE[[4, 6;;7]], outcomesAndpPCE[[3, {8}]]],
	outcomesAndpPCE[[4, {7}]],
	outcomesAndpPCE[[3, 8;;9]]},
PlotStyle -> {myColors[[4]], Directive[{myColors[[4]], Dashed}], myColors[[4]], myColors[[4]]},
PlotMarkers -> {myMarkers[[2]], None, myMarkers[[3]], myMarkers[[2]]}]


(* ::Subsection::Closed:: *)
(*Combine plots and legend, export*)


Legended[Show[linePlotCE, fkpOutcome], myLegend]
Export["evolution_from_heteroptrophy_vs_pP_light25_ce.pdf", %]
