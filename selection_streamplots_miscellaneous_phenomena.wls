#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Parameters*)


(* ::Text:: *)
(*Note on the parameter c_HP in the paper: In the code we have replaced c_HP with its reciprocal, gamma_C. This doesn't come up anywhere for the user, but just be aware that setting gamma_C = 0.7 in the functions below is equivalent to setting c_HP to 10/7. *)


(* ::Text:: *)
(*These parameters produce evolution to different ESSs at different light levels (25 vs. 1000) with conversion trade +/- attack trade off (same ESSs with either)*)
(*Removing conversion trade-off means we lose this, partly I think because of changes in ability of the mixotroph to survive/maintain decent population sizes at high retention rates*)


pW = 3;
gammaP = 2/3;(*2.9/3;*)
sigmaP = 0;
(*kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;*)
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


(* ::Section:: *)
(*Load model equations and analysis functions*)


Get[NotebookDirectory[] <> "model_and_functions.wls"]


(* ::Section:: *)
(*Effect of lower light absorptivities*)


kW = 0.2;
kH = kH2 = 0.1;
kP = kP2 = 0.3;


lightLevels = {100};


(* ::Subsection:: *)
(*No attack tradeoff*)


evolutionSimsCEhighLight = Map[(*Parallel*)Table[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6,
 conversionTradeOff -> True, attackTradeOff -> False], 1]&, lightLevels];


Table[averageEcoState[0.6 * i, i, 100, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True], {i, 0, 1, 1/10}]


Export["evolutionSimsCE_10e-6_light55_95.wdx", evolutionSimsCEhighLight];


evolutionSimsCEhighLight = Import["evolutionSimsCE_10e-6_light55_95.wdx"];


GraphicsGrid[Table[
{ListLinePlot[evolutionSimsCEhighLight[[i, All, 1]], AxesLabel -> {"Time (mutations)", "Retention"}, PlotLabel -> "Light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium],
ListLinePlot[evolutionSimsCEhighLight[[i, All, 2]], AxesLabel -> {"Time (mutations)", "Replication"}, PlotLabel -> "Light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium]},
{i, Length[lightLevels]}]]


(* ::Text:: *)
(*Simulate evolution with increased precision. (We see that extinction only occurs at higher light levels.)*)


evolutionSimsCEhighPrec = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 10, AccuracyGoal -> 10, conversionTradeOff -> True], 10]&, {65, 70, 75, 80, 85, 90, 95, 100}];


Export["evoltionSimsCE_10e-10_light65-100.wdx", evolutionSimsCEhighPrec]


evolutionSimsCEhighPrec[[All, All, All, 1001]] //MatrixForm
