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


(* ::Section::Closed:: *)
(*Equations*)


Get[NotebookDirectory[] <> "model_and_functions.wls"]


(* ::Section::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Light levels to investigate*)


(* ::Text:: *)
(*Light levels to simulate evolution from heterotrophy*)


simLightLevels = {22, 23, 24, 25, 26, 27, 28, 29, 30};
(*higherLightLevels = {55, 60, 65, 70, 75, 80, 85, 90};*)


(* ::Text:: *)
(*We also wanted to investigate slightly lower light levels. It turns out the heterotroph can't survive at these light levels, so we're putting them in a separate list.*)
(*(To avoid wasting simulation time, we'll just confirm that the heterotroph can't survive at these light levels, and thus there can be no evolution from strict heterotrophy.)*)


lowerLightLevels = {20, 21};


(* ::Text:: *)
(*Finally, we'll put all the light levels together for plotting.*)


lightLevels = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};


(* ::Text:: *)
(*A simple function to label outcomes of evolution according to the classification scheme in the paper*)


(* ::Subsection:: *)
(*Conversion efficiency trade-off*)


(* ::Subsubsection::Closed:: *)
(*Simulate evolution (light > 21 micromol photons/m^2/s)*)


evolutionSimsCE = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True], 100]&, simLightLevels];


Export["evolutionSimsCE_10e-6_light22-30_May15.wdx", evolutionSimsCE];


(*evolutionSimsCEmost = Import["evolutionSimsCE_10^-6_1-100.wdx"];
evolutionSimsCE29 = Import["evolutionSimsCE_10^-6_1-100_light29.wdx"];*)


(*evolutionSimsCE = Join[evolutionSimsCEmost[[1;;7]], evolutionSimsCE29, evolutionSimsCEmost[[{8}]]];
Export["evolutionSimsCE_10e-6_light22-30.wdx", evolutionSimsCE];*)


(*evolutionSimsCE = Import["evolutionSimsCE_10e-6_light22-30.wdx"];*)


(* ::Subsubsection::Closed:: *)
(*Simulation evolution (light > 30 micromol photons/m^2/s)*)


evolutionSimsCEhighLight = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True], 10]&, higherLightLevels];


Export["evolutionSimsCE_10e-6_light55_95.wdx", evolutionSimsCEhighLight];


evolutionSimsCEhighLight = Import["evolutionSimsCE_10e-6_light55_95.wdx"];


GraphicsGrid[Table[
{ListLinePlot[evolutionSimsCEhighLight[[i, All, 1]], AxesLabel -> {"Time (mutations)", "Retention"}, PlotLabel -> "Light = " <> ToString[higherLightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium],
ListLinePlot[evolutionSimsCEhighLight[[i, All, 2]], AxesLabel -> {"Time (mutations)", "Replication"}, PlotLabel -> "Light = " <> ToString[higherLightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium]},
{i, Length[higherLightLevels]}]]
Export["evolutionSimsCE_10e-6_light55_95.pdf", %]


(* ::Text:: *)
(*Simulate evolution with increased precision. (We see that extinction only occurs at higher light levels.)*)


evolutionSimsCEhighPrec = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 10, AccuracyGoal -> 10, conversionTradeOff -> True], 10]&, {65, 70, 75, 80, 85, 90, 95, 100}];


Export["evoltionSimsCE_10e-10_light65-100.wdx", evolutionSimsCEhighPrec]


evolutionSimsCEhighPrec[[All, All, All, 1001]] //MatrixForm


(* ::Subsubsection::Closed:: *)
(*Plot the trait values vs. time (just to see what's going on)*)


GraphicsGrid[Table[{
ListLinePlot[evolutionSimsCE[[i, All, 1]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Retention probability (phi)"}, PlotLabel -> "Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day", PlotRange -> {All, {0, 1}}],
ListLinePlot[evolutionSimsCE[[i, All, 2]], Frame -> True, FrameLabel -> {"Time (in mutations that have arisen)", "Replication probability (eta)"}, PlotLabel -> "Light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/day", PlotRange -> {All, {0, 1}}]}, 
{i, Length[simLightLevels]}], 
ImageSize -> 800]


(* ::Subsubsection::Closed:: *)
(*Simulate strict heterotroph at low light (light 20-21 micromol photons/m^2/s)*)


(* ::Text:: *)
(*Strict heterotroph never survives (extinction = True in all cases)*)


AllTrue[Table[residentEcoDynamics[0, 0, 20, 1, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False], 100][[All, 2]], #&]


AllTrue[Table[residentEcoDynamics[0, 0, 21, 1, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False], 100][[All, 2]], #&]


(* ::Text:: *)
(*Add these cases to the higher light simulated outcomes*)


evolutionSimsCE = Prepend[evolutionSimsCE, Table[{Table["NA", 1001], Table["NA", 1001]}, 100]];
evolutionSimsCE = Prepend[evolutionSimsCE, Table[{Table["NA", 1001], Table["NA", 1001]}, 100]];


(* ::Subsubsection::Closed:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


outcomesCE = Map[classifyOutcome[#[[All, -1]]]&, evolutionSimsCE, {2}];
outcomeCountsCE = Outer[Count, outcomesCE, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsCE = outcomeCountsCE / Dimensions[evolutionSimsCE][[2]];


(* ::Subsubsection::Closed:: *)
(*Plot a bar chart showing how light affects the endpoint of evolution ("outcome")*)


conversionTradeOffChart = BarChart[outcomeCountsCE, ChartLayout -> "Stacked", ChartLegends -> outcomeLabels, ChartStyle -> myColors, 
Frame -> {True, True, False, False},
(* Put the light levels on the x-axis *)
FrameTicks -> {Table[{i, lightLevels[[i]]}, {i, Length[lightLevels]}], Automatic},
FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "No trade-off"] (* This is the title because we're assuming everything has the CE trade-off, so we're just calling it normal *)


(* ::Subsection::Closed:: *)
(*Conversion efficiency and attack rate trade-offs*)


(* ::Subsubsection::Closed:: *)
(*Simulate evolution (light > 21 micromol photons/m^2/s)*)


evolutionSimsCEAR = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True], 100]&, simLightLevels];


Export["evolutionSimsCEAR_10e-6_light22-30_May15.wdx", evolutionSimsCEAR];


(*evolutionSimsCEARmost = Import["evolutionSimsCEAR_10e-6_sims1-100.wdx"];
evolutionSimsCEAR29 = Import["evolutionSimsCEAR_10e-6_sims1-100_light29.wdx"];*)


(*evolutionSimsCEAR = Join[evolutionSimsCEARmost[[1;;7]], evolutionSimsCEAR29, evolutionSimsCEARmost[[{8}]]];
Export["evolutionSimsCEAR_10e-6_light22-30.wdx", evolutionSimsCEAR];*)


(*evolutionSimsCEAR = Import["evolutionSimsCEAR_10e-6_light22-30.wdx"];*)


(* ::Subsubsection::Closed:: *)
(*Simulation evolution (light > 30 micromol photons/m^2/s)*)


evolutionSimsCEARhighLight = Map[ParallelTable[simulateEvolution[#, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True], 10]&, higherLightLevels];


Export["evolutionSimsCEAR_10e-6_light55_95.wdx", evolutionSimsCEARhighLight];


evolutionSimsCEARhighLight = Import["evolutionSimsCEAR_10e-6_light55_95.wdx"];


GraphicsGrid[Table[
{ListLinePlot[evolutionSimsCEARhighLight[[i, All, 1]], AxesLabel -> {"Time (mutations)", "Retention"}, PlotLabel -> "Light = " <> ToString[higherLightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium],
ListLinePlot[evolutionSimsCEARhighLight[[i, All, 2]], AxesLabel -> {"Time (mutations)", "Replication"}, PlotLabel -> "Light = " <> ToString[higherLightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium]},
{i, Length[higherLightLevels]}]]
Export["evolutionSimsCEAR_10e-6_light55_95.pdf", %]


(* ::Subsubsection::Closed:: *)
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


(* ::Section::Closed:: *)
(*Check whether the mixotroph population size dropped very low*)


(* ::Subsection::Closed:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection::Closed:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEAR = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEAR = Flatten[minPopSizesCEAR, {{1}, {2, 3}}];


extinctionPointsCEAR = Map[Select[#[[4]] <= 10^-6&], minPopSizesCEAR];


(* ::Text:: *)
(*Plot the points where the population goes extinct (green dots) and the evolutionary trajectories. The population doesn't go near the extinction points.*)


Transpose[evolutionSimsCEAR[[9]], {1, 3, 2}]


Table[Show[{

(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsCEAR[[i+2, All, {2, 3}]], 
PlotStyle -> {Green, PointSize[Large]},
FrameLabel -> {"Retention", "Replication"}, PlotRangePadding -> 0.02,
PlotLabel -> "(" <> ToUpperCase[Alphabet[]][[i]] <> ") Mixotroph population (cells/mL)\nSurface light = " <> ToString[simLightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", 
ImageSize -> 400, PlotRange -> {{0, 1}, {0, 1}}],

(* Plot the evolutionary trajectories*)
ListLinePlot[Transpose[evolutionSimsCEAR[[i+2]], {1, 3, 2}]]}],

{i, Length[simLightLevels]}]


(* ::Text:: *)
(*One extinction at light = 28*)


Table[{i, AnyTrue[#&][Map[AnyTrue[#&], Outer[((#2[[2]] - 0.05) <= #1[[2]] <= (#2[[2]] + 0.05)) && ((#2[[1]] - 0.05) <= #1[[1]] <= (#2[[1]] + 0.05)) &, Transpose[evolutionSimsCEAR[[9, i]]], extinctionPointsCEAR[[9, All, {2, 3}]], 1]]]},
{i, 100}]


Select[%, #[[2]]&]


(* ::Text:: *)
(*Simulation 56 goes extinct -- we'll reclassify that one*)


classifyOutcome[evolutionSimsCEAR[[9, 56, All, -1]]]


outcomeLabels


outcomeCountsCEAR


outcomeCountsCEAR[[9, 3]] = 1/100
outcomeCountsCEAR[[9, 2]] = outcomeCountsCEAR[[9, 2]] - 1/100


outcomeCountsCEAR


(* ::Text:: *)
(*Extinctions at light = 30*)


Table[{i, AnyTrue[#&][Map[AnyTrue[#&], Outer[((#2[[2]] - 0.05) <= #1[[2]] <= (#2[[2]] + 0.05)) && ((#2[[1]] - 0.05) <= #1[[1]] <= (#2[[1]] + 0.05)) &, Transpose[evolutionSimsCEAR[[11, i]]], extinctionPointsCEAR[[11, All, {2, 3}]], 1]]]},
{i, 100}]


Select[%, #[[2]]&]


(* ::Text:: *)
(*Simulations 89 and 91 go extinct: we'll reclassify*)


classifyOutcome[evolutionSimsCEAR[[11, 89, All, -1]]]
classifyOutcome[evolutionSimsCEAR[[11, 91, All, -1]]]


outcomeLabels


outcomeCountsCEAR


outcomeCountsCEAR[[11, 3]] = 2/100
outcomeCountsCEAR[[11, 2]] = outcomeCountsCEAR[[11, 2]] - 2/100


outcomeCountsCEAR


(* ::Subsection::Closed:: *)
(*Conversion trade-off only*)


(* ::Subsubsection::Closed:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCE = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCE = Flatten[minPopSizesCE, {{1}, {2, 3}}];


extinctionPointsCE = Map[Select[#[[4]] <= 10^-6&], minPopSizesCE];


(* ::Text:: *)
(*No extinction points for light at 22 or higher!*)


extinctionPointsCE


(* ::Section:: *)
(*Plot a line chart showing endpoint of evolution vs light*)


outcomesAndLightCEAR = MapThread[{#1, #2}&, {{lightLevels, lightLevels, lightLevels}, Transpose[outcomeCountsCEAR]}, 2];


linePlotCEAR = ListLinePlot[outcomesAndLightCEAR, PlotStyle -> {myColors[[1]], myColors[[2]], Black},
PlotMarkers -> Automatic,
PlotLegends -> {"Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"},
Frame -> True, FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "Simulated evolution from strict heterotrophy\n Attack rate tradeoff",
LabelStyle -> {Black, 12, "Arial"}, ImageSize -> 450]


outcomesAndLightCE = MapThread[{#1, #2}&, {{lightLevels, lightLevels, lightLevels}, Transpose[outcomeCountsCE]}, 2];


linePlotCE = ListLinePlot[outcomesAndLightCE, PlotStyle -> {myColors[[1]], myColors[[2]], Black},
PlotMarkers -> Automatic,
PlotLegends -> Placed[LineLegend[{"Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"}], {0.65, 0.5}],
Frame -> True, FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "Simulated evolution from strict heterotrophy\n No tradeoff", 
ImageSize -> 450, LabelStyle -> {Black, 12, "Arial"}]


(* ::Subsection:: *)
(*Plot all trade-off variations together*)


(* Add letter labels to plots *)
linePlots = MapThread[Show[#1, PlotLabel ->  #2 <> (PlotLabel /. Options[#1[[1]]])]&, {{linePlotCE, linePlotCEAR}, {"(A) ", "(B) "}}];

(* Remove legend from the second plot *)
linePlots = {linePlots[[1]], linePlots[[2, 1]]};

(* Put the plots in a row with a single, shared legend *)
Row[linePlots, Spacer[30]]

(* Export as pdf *)
Export["evo_from_heterotrophy_old_mutation_model.pdf", %];


(* ::Text:: *)
(*To make the final figure, we changed the plots to be arranged vertically instead of horizontally in an image editor. *)
