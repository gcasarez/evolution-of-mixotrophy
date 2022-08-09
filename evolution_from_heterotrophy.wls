#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Parameters*)


(* ::Text:: *)
(*Note on the parameter c_HP in the paper: In the code we have replaced c_HP with its reciprocal, gamma_C. This doesn't come up anywhere for the user, but just be aware that setting gammaC = 0.7 in the functions below is equivalent to setting c_HP to 10/7. *)


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
	0 < replication < 1, "No label"] (* We don't have a name for any state where replication is intermediate *) 


(* ::Text:: *)
(*A list of the outcomes we observe (we'll still be able to tell if things didn't get one of these labels, though)*)


(* outcomeLabels = {"Heterotroph", "Strict Kleptoplast", "Phototroph", "Kleptoplastidic Phototroph", "Occasional Kleptoplastidic Phototroph", "Extinction"}; *)
outcomeLabels = {(*"Strict kleptoplast", *)"Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"};


(* ::Text:: *)
(*Colors for plotting*)


(* myColors = "CandyColors"; *)
myColors = Append[(ColorData["CandyColors"] /@ Subdivide[3])[[{(*1,*) 2, 4}]], White];


(* ::Section:: *)
(*Run simulations*)


(* ::Subsection:: *)
(*Light levels to investigate*)


(* ::Text:: *)
(*Light levels to simulate evolution from heterotrophy*)


simLightLevels = {35, 55, 85};
(*higherLightLevels = {55, 60, 65, 70, 75, 80, 85, 90};*)


(* ::Text:: *)
(*We also wanted to investigate slightly lower light levels. It turns out the heterotroph can't survive at these light levels, so we're putting them in a separate list.*)
(*(To avoid wasting simulation time, we'll just confirm that the heterotroph can't survive at these light levels, and thus there can be no evolution from strict heterotrophy.)*)


lowerLightLevels = {20, 21};


(* ::Text:: *)
(*Finally, we'll put all the light levels together for plotting.*)


lightLevels = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};


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


(* ::Subsection:: *)
(*Conversion efficiency and attack rate trade-offs*)


(* ::Subsubsection:: *)
(*Simulate evolution (light > 21 micromol photons/m^2/s)*)


evolutionSimsCEAR = ParallelMap[{#, simulateEvolution[#, conversionTradeOff -> True, attackTradeOff -> True]}&, 
	Table[light, {light, simLightLevels}, 10]];


Export["evolutionSimsCEAR_10e-6_Aug9.wdx", evolutionSimsCEAR];


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


(* ::Section:: *)
(*Check whether the mixotroph population size dropped very low*)


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
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


(* ::Section::Closed:: *)
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


(* ::Section:: *)
(*Figures: evolutionary simulations*)


(* ::Subsection:: *)
(*Plot all trade-off variations together*)


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


(* ::Section:: *)
(*Figures: example ecology and evolution plots*)


(* ::Subsection:: *)
(*Colors and fonts*)


(* ::Text:: *)
(*Colors for plotting mixotroph population dynamics*)


colorW = RGBColor["#77aa77"];
colorH = colorH2 = RGBColor["#fdc552"];
colorP = colorP2 = RGBColor["#4573ba"];
colorM = RGBColor["#cd1aff"];
colorRetention = RGBColor["#a02b2b"];
colorReplication = RGBColor["#2c89a0"];


(* ::Text:: *)
(*Label font options (this saves so much trouble to specify here!*)


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


(* ::Subsubsection:: *)
(*Plot prey ecological dynamics*)


preyPlot = Plot[w[t-100] /. withMutantEcology, {t, 100, 4100}, PlotRange -> {{0, 4000}, {0, 20}}, PlotStyle -> colorW, Frame -> True, 
GridLines -> {{{100, Directive[Thick, Black, Dashed]}}, None},
FrameLabel -> {{"Prey (cells/mL)", None}, {"Time (days)", None}}, ImageSize -> 350, ImagePadding -> {{50, 50}, {50, 10}}, PlotRangePadding -> 0];


preyPlot = Show[{preyPlot, Plot[w[t+3900] /. residentEcology, {t, 0, 100}, PlotStyle -> colorW], Graphics[{Gray, Opacity[0.2], Rectangle[{3000, 0}, {4000, 20}]}]}]


(* ::Subsubsection:: *)
(*Plot mixotroph ecological dynamics*)


(* ::Text:: *)
(*Plot legend*)


(*mixotrophPlotLegend = LineLegend[{colorH, Directive[colorH2, Dashed], colorP, Directive[colorP2, Dashed]}, {"Resident heterotroph", "Mutant heterotroph", "Resident phototroph", "Mutant phototroph"}];*)


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


(* ::Text:: *)
(*This code plots the total mixotroph population (combined phototrophs and heterotrophs)*)


mixotrophPlot = Plot[{(h[t - 100] + p[t - 100]) /. withMutantEcology // Chop, (h2[t - 100] + p2[t - 100]) /. withMutantEcology // Chop}, {t, 100, 4100}, 
GridLines -> {{{100, Directive[Thick, Black, Dashed]}}, None},
PlotRange -> {{0, 4000}, {0, 0.6}}, PlotStyle -> {Directive[Thick, colorM], Directive[Thick, Dashed, colorM]}];


mixotrophPlot = Show[{mixotrophPlot, 
		Plot[(h[t + 3900] + p[t + 3900]) /. residentEcology // Chop, {t, 0, 100}, PlotRange -> {{0, 1000}, {0, 0.5}}, PlotStyle -> Directive[Thick, colorM]],
		Graphics[{Gray, Opacity[0.2], Rectangle[{3000, 0}, {4000, 20}]}]},
	Frame -> True,
	FrameLabel -> {{"Consumer biomass (cells/mL)", None}, {"Time (days)", None}}, 
	ImageSize -> 350, ImagePadding -> {{50, 50}, {50, 10}}, PlotRangePadding -> 0]


(* ::Subsubsection:: *)
(*Plot the proportion of the population that is the mutant*)


mutantProportionPlot = Plot[((h2[t - 100] + p2[t - 100])/ (h[t - 100] + p[t - 100] + h2[t - 100] + p2[t - 100]))/. withMutantEcology, {t, 100, 4100}, PlotRange -> {{0, 4000}, {0, 1}}, PlotStyle -> Black, Frame -> True,
GridLines -> {{{100, Directive[Thick, Black, Dashed]}}, None},
FrameLabel -> {"Time (days)", "Mutant proportion of population"}, ImageSize -> 350, ImagePadding -> {{50, 50}, {50, 10}}, PlotRangePadding -> 0];


mutantProportionPlot = Show[{mutantProportionPlot, Plot[0, {t, 0, 100}, PlotStyle -> Black], 
		Graphics[{Gray, Opacity[0.2], Rectangle[{3000, 0}, {4000, 20}]}]}];


(* ::Subsubsection:: *)
(*Combine plots into a list*)


ecologyPlots = {Labeled[Show[preyPlot, LabelStyle -> {Black, "Arial", 12}], "(A) Prey ecology", Top, Alignment -> {Left, Baseline}],
Labeled[Show[mixotrophPlot, LabelStyle -> {Black, "Arial", 12}],  "(B) Mixotroph ecology", Top, Alignment -> Left],
Labeled[Show[mutantProportionPlot, LabelStyle -> {Black, "Arial", 12}], "(C) Mutant spread", Top, Alignment -> Left]}


(* ::Subsection:: *)
(*Evolution plots*)


evolutionSimsCEAR = Import["evolutionSimsCEAR_10e-6_light22-30_May15.wdx"];


simLightLevels[[4]]


(* ::Text:: *)
(*Plot the plastid retention probability over evolutionary time*)


retentionPlot = ListLinePlot[evolutionSimsCEAR[[4, 5, 1]], Frame -> True, FrameLabel -> {"Time (mutations that have arisen)", "Retention"}, ImageSize -> 350, 
ImagePadding -> {{50, 50}, {50, 10}}, PlotStyle -> colorRetention, PlotRange -> {All, {0, 1}}, PlotRangePadding -> 0]


(* ::Text:: *)
(*Plot the plastid replication probability over evolutionary time*)


replicationPlot = ListLinePlot[evolutionSimsCEAR[[4, 5, 2]], Frame -> True, FrameLabel -> {"Time (mutations that have arisen)", "Replication"}, ImageSize -> 350, 
ImagePadding -> {{50, 50}, {50, 10}}, PlotStyle -> colorReplication, PlotRange -> {All, {0, 1}}, PlotRangePadding -> 0]


(* ::Text:: *)
(*Plot the evolutionary trajectory in trait space so we can see what all the different (retention, replication) combinations were*)


traitSpacePlot = ListLinePlot[Transpose[evolutionSimsCEAR[[4, 5]]], Frame -> True, FrameLabel -> {"Retention", "Replication"}, ImageSize -> 350, 
ImagePadding -> {{50, 50}, {50, 10}}, PlotStyle -> Black, PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0]


(* ::Text:: *)
(*Find and plot the "equilibrium" population sizes over evolutionary time*)


averagePopSizes = ParallelTable[Module[{residentEcology, extinction, retention = evolutionSimsCEAR[[4, 5, 1, i]], replication = evolutionSimsCEAR[[4, 5, 2, i]]}, 
	{residentEcology, extinction} = residentEcoDynamics[retention, replication, 25, 4, 1, 1, conversionTradeOff -> True, attackTradeOff -> True, PrecisionGoal -> 6, AccuracyGoal -> 6];
	Table[
		If[extinction, {"Unknown", 0},
			NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000)], 
	{stateVar, {w[t], (h[t] + p[t])}}]],
{i, Length[evolutionSimsCEAR[[4, 1, 1]]]}];


averagePopPlot = ListLinePlot[{averagePopSizes[[All, 1]], averagePopSizes[[All, 2]]}, PlotStyle -> {colorW, Directive[Thick, colorM]},
Frame -> True, FrameLabel -> {"Time (mutations that have arisen)", "Population size (cells/mL)"}, ImageSize -> 350,
ImagePadding -> {{50, 50}, {50, 10}}, PlotStyle -> Black, PlotRange -> {{0, 1000}, {0, 17}}, PlotRangePadding -> 0];


(* ::Subsubsection:: *)
(*Combine plots into a list*)


evolutionPlots = {Labeled[Show[averagePopPlot, LabelStyle -> {Black, "Arial", 12}], "(D) Average population size", Top, Alignment -> Left],
Labeled[Show[replicationPlot, LabelStyle -> {Black, "Arial", 12}], "(E) Replication evolution", Top, Alignment -> Left],
Labeled[Show[retentionPlot, LabelStyle -> {Black, "Arial", 12}], "(F) Retention evolution", Top, Alignment -> {Left, Baseline}],
Labeled[Show[traitSpacePlot, LabelStyle -> {Black, "Arial", 12}], "(G) Evolutionary trajectory in trait space", Top, Alignment -> Left]}


(* ::Subsection:: *)
(*Combine ecology and evolution plots into the figure*)


myArrow = Graphics[{Thick, Line[{{1/4, 1/8}, {1/4, 3/8}, {5/8, 0}, {1/4, -3/8}, {1/4, -1/8}}],
   Line[{{1/4, 1/8}, {-1/2, 1/8}, {-1/2, -1/8}, {1/4, -1/8}}]}];


mixotrophPlotLegend = Framed[ (* put a frame around the legend *)
	LineLegend[{colorW, Directive[Thick, colorM], Directive[Thick, Dashed, colorM]}, {"Prey", "Resident consumer", "Mutant consumer"}, LegendLabel -> "Population",
		LabelStyle -> {14, "Arial", Black}], 
(* frame options *)
RoundingRadius -> 10, FrameStyle -> {Thin, Black}];


Grid[Transpose[{Prepend[ecologyPlots, Style["Ecological Dynamics", FontFamily -> "Arial", FontColor -> Black, FontSize -> 12]],
{"", "", Show[myArrow, ImageSize -> 50], ""},
Prepend[evolutionPlots[[1;;3]], Style["Evolutionary Dynamics", FontFamily -> "Arial", FontColor -> Black, FontSize -> 12]], 
{"", "", Invisible[Show[myArrow, ImageSize -> 50]], ""},
{"", mixotrophPlotLegend, evolutionPlots[[4]], ""}}],
Alignment -> Center, Spacings -> 1]

Export["example_timeseries.pdf", %]



