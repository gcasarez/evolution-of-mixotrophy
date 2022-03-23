#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*Plot settings*)


(* ::Text:: *)
(*Colors for plotting prey and mixotrophs*)


colorW = RGBColor["#77aa77"]; (* prey *)
colorM = RGBColor["#cd1aff"]; (* mixotroph *)


(* ::Text:: *)
(*Font for plots*)


Map[SetOptions[#, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black}]&,
	{Plot, ListPlot, ListLinePlot, LineLegend}];


(* ::Section::Closed:: *)
(*Model Equations*)


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
(*Functions to Find Equilibria*)


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


(* ::Subsection::Closed:: *)
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


(* ::Section::Closed:: *)
(*Find Stable Equilibria!*)


(* ::Text:: *)
(*We want to make a bifurcation diagram where light = 35 micromol photons/m^2/s, replication = 0.2, and there is no attack trade-off *)


(* ::Text:: *)
(*We'll show how the equilibrium population changes with the retention probability, by looking across a variety of retention values.*)


retentionValues = Table[retention, {retention, 0, 1, 1/20}];


(* ::Subsubsection::Closed:: *)
(*Step 1: Find places where the derivative = 0 (these could be stable or unstable equilibria). *)


(* ::Text:: *)
(*We look across a variety of initial population sizes (guessW, guessH, guessP), in case there are multiple equilibria.*)


candidateStableEquils = Table[Check[ (* Check is here to prevent answers that don't meet the specified accuracy/precision from being returned *)
equilibrium[35, retention, 2/10, {guessW, guessH, guessP}, attackTradeOff -> False, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6, MaxIterations -> 10000],
{}, FindRoot::lstol],
{retention, retentionValues}, {guessW, {1/1000, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40}}, {guessH, {1/1000, 1, 2, 3, 4, 5}}, {guessP, {0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40}}];


(* ::Text:: *)
(*Reformat the list of candidate stable equilibria so that all the equilibria for a particular retention probability are listed together.*)


candidateStableEquils = Flatten[candidateStableEquils, {{1}, {2, 3, 4}, {5}}];


(* ::Text:: *)
(*Remove any "candidate" equilibria that are just the empty set*)


candidateStableEquils = Map[Select[Length[#] > 0&], candidateStableEquils];


(* ::Text:: *)
(*Remove any duplicate equilibria. (Sometimes multiple initial population sizes will lead the solver to the same equilibrium.)*)


duplicatesRemoved = Table[DeleteDuplicates[candidateStableEquils[[i]], EuclideanDistance[({w, h, p} /. #1), ({w, h, p} /. #2)] < 1/1000&], {i, Length[candidateStableEquils]}];


(* ::Subsubsection::Closed:: *)
(*Step 2: Find the stability of the equilibria*)


(* ::Text:: *)
(*Calculate the Jacobian for each candidate equilibrium.*)


jacobiansList = Table[jacobian[35, retentionValues[[i]], 2/10, attackTradeOff -> False, conversionTradeOff -> True] /. duplicatesRemoved[[i]], {i, Length[candidateStableEquils]}];


(* ::Text:: *)
(*If the real parts of all eigenvalues of the Jacobian are negative, the point should be a stable equilibrium.*)
(*Make a list which has the stability of each equilibrium. (True here means stable.)*)


stability = Map[AllTrue[Re[Eigenvalues[#]], Negative]&, jacobiansList, {2}];


(* ::Text:: *)
(*Select all the equilibria that are stable.*)


stableEquils = Pick[duplicatesRemoved, stability, True];


(* ::Text:: *)
(*Check how many stable equilibria appear at each retention probability*)


MapThread[{#1, #2}&, {retentionValues,  Map[Dimensions, stableEquils][[All, 1]]}] (* Stuff after this is completely unnecessary, just prints the output nicely *) //
	TableForm[#, TableDirections -> Row, TableHeadings -> {None, {"Retention", "# Equilibria"}}]&


(* ::Subsubsection::Closed:: *)
(*Step 3: Make some temporary plots*)


(* ::Text:: *)
(*Make a list of colors that corresponds to whether the equilibrium is stable (blue) or unstable (red).*)


stabilityColors = Map[If[#, Directive[Blue, Opacity[1]], Directive[Red, Opacity[0.5]]]&, stability, {2}];


(* ::Text:: *)
(*Plot all equilibria, colored by stability*)


GraphicsRow[{ (* GraphicsRow is unnecessary, I just think it makes viewing these plots a bit nicer. *)

(* Equilibrium prey density *)
ListPlot[Table[MapThread[Style[{retentionValues[[i]], #1}, #2]&, {w /. duplicatesRemoved[[i]], stabilityColors[[i]]}], {i, Length[candidateStableEquils]}], PlotRange -> All,
	PlotLabel -> "Prey", AxesLabel -> {"Retention", "Prey"}, ImageSize -> Medium, PlotStyle -> {PointSize -> Medium}],

(* Equilibrium heterotroph density *)
ListPlot[Table[MapThread[Style[{retentionValues[[i]], #1}, #2]&, {h /. duplicatesRemoved[[i]], stabilityColors[[i]]}], {i, Length[candidateStableEquils]}], PlotRange -> All,
	PlotLabel -> "Heterotroph", AxesLabel -> {"Retention", "Heterotroph"}, ImageSize -> Medium, PlotStyle -> {PointSize -> Medium}],
	
(* Equilibrium phototroph density *)
ListPlot[Table[MapThread[Style[{retentionValues[[i]], #1}, #2]&, {p /. duplicatesRemoved[[i]], stabilityColors[[i]]}], {i, Length[candidateStableEquils]}], PlotRange -> All,
	PlotLabel -> "Phototroph", AxesLabel -> {"Retention", "Phototroph"}, ImageSize -> Medium, PlotStyle -> {PointSize -> Medium}],
	
(* Plot legend *)
PointLegend[{Blue, Directive[Red, Opacity[0.5]]}, {"Stable", "Unstable"}, LegendLabel -> "Equilibrium is...", LabelStyle -> 12, LegendMarkerSize -> Medium]

}]


(* ::Text:: *)
(*What we can see from the above plot is that above a retention probability of about 0.25, the stable equilibria seem to occur at negative population sizes. *)


(* ::Subsubsection::Closed:: *)
(*Step 4: Find nonnegative stable equilibria*)


(* ::Text:: *)
(*Determine whether the prey, heterotroph, and mixotroph population sizes are all nonnegative for each stable equilibrium.*)
(*We assume populations that are negative and very close to 0 (within 10^-10) are actually 0, and round them to 0 before we do the selection.*)


nonNegativePopulations = Map[AllTrue[NonNegative], {w, h, p} /. Chop[stableEquils], {2}];


(* ::Text:: *)
(*Select the stable equilibria with nonnegative populations.*)


nonNegativeEquils = Pick[stableEquils, nonNegativePopulations];


(* ::Text:: *)
(*Check how many nonnegative stable equilibria appear at each retention probability*)


MapThread[{#1, #2}&, {retentionValues,  Map[Dimensions, nonNegativeEquils][[All, 1]]}] //
	TableForm[#, TableDirections -> Row, TableHeadings -> {None, {"Retention", "# Equilibria"}}]&


MapThread[N[Chop[{#1, #2, #3, #4}]]&, 
{retentionValues,
Pick[w/.stableEquils, nonNegativePopulations],
Pick[h/.stableEquils, nonNegativePopulations],
Pick[p/.stableEquils, nonNegativePopulations]}] //TableForm[#, TableHeadings -> {None, {"Retention", "Prey", "Heterotroph", "Phototroph"}}]&


(* ::Subsubsection::Closed:: *)
(*Step 5: Investigate points where no positive stable equilibria were found*)


(* ::Text:: *)
(*Determine which stable equilibria have negative populations or all populations essentially 0*)
(*(We round anything within 10^-10 of 0 to 0.)*)


nonPositivePopulations = Map[AnyTrue[#, Negative] || AllTrue[#, NonPositive]&, {w, h, p} /. Chop[stableEquils], {2}];


(* ::Text:: *)
(*Determine which retention probabilities only have nonpositive stable equilibria*)


retentionNoPositiveEquil = Pick[retentionValues, Map[AllTrue[#, #1&]&, nonPositivePopulations]]


(* ::Text:: *)
(*For each of these retention probabilities, simulate population dynamics*)


simulatedPopDynamics = Table[residentEcoDynamics[retention, 2/10, 35, 4, 1, 1, attackTradeOff -> False, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6],
{retention, retentionNoPositiveEquil}];


(* ::Text:: *)
(*Check whether any of these populations go extinct*)
(*(They don't, so we'll move on to plotting them below)*)


MapThread[{#1, #2}&, {retentionNoPositiveEquil,  simulatedPopDynamics[[All, 2]]}] //
	TableForm[#, TableDirections -> Row, TableHeadings -> {None, {"Retention", "Extinction?"}}]&


(* ::Text:: *)
(*Plot population dynamics*)


Table[Plot[Evaluate[{w[t], h[t], p[t]}/.simulatedPopDynamics[[i, 1]]], {t, 3000, 4000},
PlotLabel -> "Retention = " <> ToString[N[retentionNoPositiveEquil[[i]]]], 
AxesLabel -> {"Time (days)", "Population size"}, ImageSize -> Medium],
{i, Length[simulatedPopDynamics]}];

(* This is completely unnecessary. It just arranges things nicely in a grid with a legend. *)
GraphicsGrid[Partition[
	Prepend[%, LineLegend[ColorData[97, "ColorList"][[1;;3]], {"Prey", "Heterotroph", "Phototroph"}, LegendLabel -> "Population", LabelStyle -> 12]], 
	UpTo[4]]]


(* ::Text:: *)
(*All of these plots show cyclical population dynamics.*)


(* ::Subsubsection::Closed:: *)
(*(Just to check, simulate population dynamics for the retention values with positive stable equilibria)*)


(* ::Text:: *)
(*Get the retention values with positive, stable equilibria.*)


retentionPositiveEquil = Pick[retentionValues, Map[AllTrue[#, #1&] == False&, nonPositivePopulations]]


(* ::Text:: *)
(*Simulate population dynamics.*)


simulatedPopDynamicsPositiveEquil = Table[residentEcoDynamics[retention, 2/10, 35, 4, 1, 1, attackTradeOff -> False, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6],
{retention, retentionPositiveEquil}];


(* ::Text:: *)
(*Plot population dynamics.*)
(*For retention < 0.25, everything reaches a stable equilibrium within 3000 time steps.*)
(*For retention = 0.25, we see damped oscillations, but they look like they're converging on the equilibrium value we found above (prey = 3.6, heterotroph = 3.3, phototroph = 1).*)


Table[Plot[Evaluate[{w[t], h[t], p[t]}/.simulatedPopDynamicsPositiveEquil[[i, 1]]], {t, 3000, 4000},
PlotLabel -> "Retention = " <> ToString[N[retentionPositiveEquil[[i]]]], 
AxesLabel -> {"Time (days)", "Population size"}, ImageSize -> Medium],
{i, Length[simulatedPopDynamicsPositiveEquil]}];

(* This is completely unnecessary. It just arranges things nicely in a grid with a legend. *)
GraphicsGrid[Partition[
	Prepend[%, LineLegend[ColorData[97, "ColorList"][[1;;3]], {"Prey", "Heterotroph", "Phototroph"}, LegendLabel -> "Population", LabelStyle -> 12]], 
	UpTo[4]]]


(* ::Subsubsection::Closed:: *)
(*Find the peak and trough of population cycles*)


(* ::Text:: *)
(*Find the minimum and maximum prey population sizes.*)


preyCycle = Table[
{NMinimize[{w[t]/.simulatedPopDynamics[[i, 1]], 3000 <= t <= 4000}, t, AccuracyGoal -> 6, PrecisionGoal -> 6][[1]],
NMaximize[{w[t]/.simulatedPopDynamics[[i, 1]], 3000 <= t <= 4000}, t, AccuracyGoal -> 6, PrecisionGoal -> 6][[1]]},
{i, Length[simulatedPopDynamics]}]


(* ::Text:: *)
(*Find the minimum and maximum total mixotroph population sizes.*)


mixotrophCycle = Table[
{NMinimize[{(h[t] + p[t])/.simulatedPopDynamics[[i, 1]], 3000 <= t <= 4000}, t, AccuracyGoal -> 6, PrecisionGoal -> 6][[1]],
NMaximize[{(h[t] + p[t])/.simulatedPopDynamics[[i, 1]], 3000 <= t <= 4000}, t, AccuracyGoal -> 6, PrecisionGoal -> 6][[1]]},
{i, Length[simulatedPopDynamics]}]


(* ::Text:: *)
(*Note that the phototroph and heterotroph cycles are slightly offset (one peaks slightly after the other), so the maximum is slightly lower than maximum heterotroph + maximum photroph.*)
(* Likewise, the minimum total mixotroph population is slightly higher than minimum heterotroph + minimum phototroph.*)


(* Calculate max & min pop. size for heterotroph and phototroph separately *)
Table[
population -> {NMinimize[{population[t]/.simulatedPopDynamics[[i, 1]], 3000 <= t <= 4000}, t, AccuracyGoal -> 6, PrecisionGoal -> 6][[1]],
NMaximize[{population[t]/.simulatedPopDynamics[[i, 1]], 3000 <= t <= 4000}, t, AccuracyGoal -> 6, PrecisionGoal -> 6][[1]]},
{i, Length[simulatedPopDynamics]}, {population, {h, p}}];

(* Compare to the max. and min. total mixotroph population *)
mixotrophCycle - ((h + p) /. %) //
	TableForm[#, TableDirections -> Row, TableHeadings -> {None, {"Mixotroph min. - (Heterotroph min. + Photroph min.)", "Mixotroph max. - (Heterotroph max. + Photroph max.)"}}]&


(* ::Subsubsection::Closed:: *)
(*Plot*)


(* ::Text:: *)
(*Get the positive, stable equilibria.*)


preyPositiveEquil = Pick[w /. stableEquils, Map[Not, nonPositivePopulations, {2}]];
mixotrophPositiveEquil = Pick[(h + p) /. stableEquils, Map[Not, nonPositivePopulations, {2}]];

(* Remove the entries that are empty because there is no positive, stable equilibrium at that retention value*)
preyPositiveEquil = Select[preyPositiveEquil, Length[#] > 0&];
mixotrophPositiveEquil = Select[mixotrophPositiveEquil, Length[#] > 0&];


(* ::Text:: *)
(*Pair the equilibria with their retention probabilities.*)


preyPositiveEquilVsRetention = MapThread[
	(* If there are multiple stable equilibria for a particular retention value, Map[{retention, #}&, equil] pairs each of the equilibria with the retention value. *)
	(* The Function[{retention, equil}] just allows us to refer to the entries of rententionPositiveEquil as "retention" and the entries of preyPositiveEquil as "equil" *)
	Function[{retention, equil}, Map[{retention, #}&, equil]],
	{retentionPositiveEquil, preyPositiveEquil}];
	
(* Make everything have nice dimensions *)
preyPositiveEquilVsRetention = Flatten[preyPositiveEquilVsRetention, 1];


mixotrophPositiveEquilVsRetention = MapThread[Function[{retention, equil}, Map[{retention, #}&, equil]], {retentionPositiveEquil, mixotrophPositiveEquil}];

(* Make everything have nice dimensions *)
mixotrophPositiveEquilVsRetention = Flatten[mixotrophPositiveEquilVsRetention, 1];


(* ::Text:: *)
(*Pair the cycle peak and troughs with their retention probabilities*)


preyCycleMinimaVsRetention = MapThread[{#1, #2}&, {retentionNoPositiveEquil, preyCycle[[All, 1]]}];

preyCycleMaximaVsRetention = MapThread[{#1, #2}&, {retentionNoPositiveEquil, preyCycle[[All, 2]]}];


mixotrophCycleMinimaVsRetention = MapThread[{#1, #2}&, {retentionNoPositiveEquil, mixotrophCycle[[All, 1]]}];

mixotrophCycleMaximaVsRetention = MapThread[{#1, #2}&, {retentionNoPositiveEquil, mixotrophCycle[[All, 2]]}];


(* ::Text:: *)
(*Join the stable equilibrium points to both cycle lists (they are both maxima and minima).*)
(*This is just for making the points connect nicely on the plot.*)


preyPlotPointsMinima = Join[preyPositiveEquilVsRetention, preyCycleMinimaVsRetention];
preyPlotPointsMaxima = Join[preyPositiveEquilVsRetention, preyCycleMaximaVsRetention];

mixotrophPlotPointsMinima = Join[mixotrophPositiveEquilVsRetention, mixotrophCycleMinimaVsRetention];
mixotrophPlotPointsMaxima = Join[mixotrophPositiveEquilVsRetention, mixotrophCycleMaximaVsRetention];


(* ::Text:: *)
(*Time to put these points into a plot!*)


tempBifurcationPlot = ListLinePlot[{preyPlotPointsMinima, preyPlotPointsMaxima,
mixotrophPlotPointsMinima, mixotrophPlotPointsMaxima}, PlotStyle -> {colorW, colorW, colorM, colorM},
PlotLabel -> "(A) Equilibrium population size vs. retention probability",
Frame -> {True, True, None, None},
FrameLabel -> {"Retention probability (\[Phi])", "Population size (cells/mL)"},
PlotLegends -> LineLegend[{colorW, colorM}, {"Prey", "Mixotroph"}],
AspectRatio -> 1/3, ImageSize -> 600, PlotRange -> {0, 26}]


(* ::Section:: *)
(*Make final plots*)


(* ::Subsubsection:: *)
(*Example population dynamics subpanels*)


(* ::Text:: *)
(*We'll make two example plots of the population dynamics: one at a retention probability that induces cycling, and one at a retention probability that does not.*)


(* ::Text:: *)
(*We'll use a retention probability of 0.1 for the plot without cycling and a retention probability of 0.6 for the plot with cycling.*)


retentionValues[[{3, 13}]]


(* ::Text:: *)
(*We'll mark these retention values on a plot so we can show them on the final figure.*)
(*For the y values of the points, we'll use the average of the population sizes.*)


pointsOfInterest = ListPlot[
Table[{{retentionValues[[i]], (* x value, retention probability *)
	Mean[{preyPlotPointsMinima[[i, 2]], preyPlotPointsMaxima[[i, 2]], mixotrophPlotPointsMinima[[i, 2]], mixotrophPlotPointsMaxima[[i, 2]]}]}}, (* y value, average of prey & mixotroph population sizes *)
	{i, {3, 13}}],
PlotMarkers -> {{"B", 15}, {"C", 15}}, PlotStyle -> Black];


(* ::Text:: *)
(*Add the marks to the bifurcation plot*)


bifurcationPlot = Show[tempBifurcationPlot, pointsOfInterest]


(* ::Text:: *)
(*Simulate the population dynamics for our two example cases.*)


examplePopDynamics = Table[residentEcoDynamics[retention, 2/10, 35, 4, 1, 1, attackTradeOff -> False, conversionTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6],
{retention, retentionValues[[{3, 13}]]}];


(* ::Text:: *)
(*Plot the population dynamics*)


examplePopPlots = Table[Show[
	Plot[Evaluate[{w[t], h[t] + p[t]} /. examplePopDynamics[[i, 1]]], {t, 0, 4000},
		PlotRange -> {{0, 4000}, {0, 26}},
		Frame -> {True, True, False, False}, FrameLabel -> {"Time (days)", "Population size"}, 
		PlotStyle -> {colorW, colorM}, 
		AspectRatio -> 1/4, ImageSize -> 800,
		
		(* The plot label. It's in a pane to get the title aligned properly :( *)
		PlotLabel -> Pane[{"(B) ", "(C) "}[[i]] <> "Population dynamics, retention = " <> ToString[N[retentionValues[[{3, 13}]][[i]]]], 
			(* Pane options. The image size here is what is needed to align the title with the bifurcation plot's title. Because the text is left-aligned, a smaller ImageSize moves the text more to the right. *)
			Alignment -> Left, ImageSize -> 531]],
		
	(* Add a gray rectangle to mark the part that is averaged for further analysis *)
	Graphics[{Gray, Opacity[0.3], Rectangle[{3000, 0}, {4000, 25}]}]],
{i, 2}]


(* ::Subsubsection:: *)
(*Assemble the plots into a figure*)


(* ::Text:: *)
(*Make a slightly nicer legend to replace the one we were using for the bifurcation plot*)


finalPlotLegend = Framed[ (* put a frame around the legend *)
	LineLegend[{colorW, colorM}, {"Prey", "Mixotroph"}, LegendLabel -> "Population",
		LabelStyle -> {14, "Arial", Black}], 
(* frame options *)
RoundingRadius -> 10, FrameStyle -> {Thin, Black}];


(* ::Text:: *)
(*Put the plots and new legend together in a grid*)


Grid[{{bifurcationPlot[[1]] (* bifurcationPlot[[1]] selects the plot without its legend *), finalPlotLegend},
{examplePopPlots[[1]] (* this just normally selects the first of two plots in examplePopPlots *), SpanFromLeft},
{examplePopPlots[[2]], SpanFromLeft}},
Spacings -> {0, 3}]
Export["bifurcation_diagram.pdf", %]



