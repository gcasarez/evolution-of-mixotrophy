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
lH = lH2 = 0.25;
lP = lP2 = 0.25;
dMax = 0.3;(*0.3;*)
sigmaD = 0; (* No more decay trade-off *)
a = a2 = 0.15;
gammaA = 1; (* Phototroph eats! *)
(*sigmaA = 0;*) (* 0 for no attack rate trade-off, 1 for attack rate trade-off *)
c = 0.1;
(*gammaC = 1;*) (* 1 for no conversion efficiency trade-off, 0.7 for conversion efficiency trade-off *)
sigmaC = 0; (* Added sigmaC because it looks like it's not defined here *)
(*iIn = 25; *)(* Incident light *)


(* ::Section::Closed:: *)
(*Equations*)


Get[NotebookDirectory[] <> "model_and_functions.wls"]


(* ::Section:: *)
(*Phototroph eats!*)


lightLevels = {25, 35, 50};


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEAR = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEAR = Flatten[ecoDynamicsCEAR, {{1}, {2,3}, {4}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEAR = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEAR = Flatten[minPopSizesCEAR, {{1}, {2, 3}}];


extinctionPointsCEAR = Table[Pick[ecoDynamicsCEAR[[i]], Map[#[[4]] <= 10^-6&, minPopSizesCEAR[[i]]]], {i, Length[ecoDynamicsCEAR]}];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEAR = Table[Pick[ecoDynamicsCEAR[[i]], Map[#[[4]] > 10^-6&, minPopSizesCEAR[[i]]]], {i, Length[ecoDynamicsCEAR]}];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCEAR = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEAR, {2}];


Export["phototroph_eats_CEAR.wdx", derivsForPlotsCEAR]


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCE = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
conversionTradeOff -> True, attackTradeOff -> False]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCE = Flatten[ecoDynamicsCE, {{1}, {2,3}, {4}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCE = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCE = Flatten[minPopSizesCE, {{1}, {2, 3}}];


extinctionPointsCE = Table[Pick[ecoDynamicsCE[[i]], Map[#[[4]] <= 10^-6&, minPopSizesCE[[i]]]], {i, Length[ecoDynamicsCE]}];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCE = Table[Pick[ecoDynamicsCE[[i]], Map[#[[4]] > 10^-6&, minPopSizesCE[[i]]]], {i, Length[ecoDynamicsCE]}];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCE = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCE, {2}];


Export["phototroph_eats_CE.wdx", derivsForPlotsCE]


(* ::Section:: *)
(*Plots*)


(* ::Text:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions -> {"OperatorSubstitution"->False}}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions->{"OperatorSubstitution"->False}}];


(* ::Text:: *)
(*Color scheme (prey density) and legend*)


maxPrey = Max[{Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]],
	Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]],
	40 (* This is the maximum in the other figures, add it here to put everything on the same scale across figures*)}];


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


preyDensityLegend = BarLegend[{preyDensityColors[Exp[#] - 1]&, {Log[1], Log[1 + maxPrey]}}, 
Ticks -> Table[{Log[i + 1], ToString[i]}, {i, {0, 1, 2.5, 5, 10, 15, 25, 40}}],
LegendLayout -> "Column", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {25, 450}];


(* ::Text:: *)
(*Make plots*)


streamPlotsCE = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") No tradeoff\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/" <> 
		ToString[Superscript["m", "2"], StandardForm] <> "/s", ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCE[[i]]], StreamStyle -> Black, StreamColorFunction -> None(*, StreamPoints \[Rule] noExtinctionPointsCE[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsCE]}]

streamPlotsCEAR = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1,2]], #[[3, 1]]}&, derivsForPlotsCEAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i + Length[lightLevels]]]] <> ") Attack tradeoff\nSurface light = " <> 
		ToString[lightLevels[[i]]] <> "\[Mu]mol photons/" <> ToString[Superscript["m", "2"], StandardForm] <> "/s", 
	ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCEAR[[i]]], StreamStyle -> Black, StreamColorFunction -> None (*, StreamPoints \[Rule] noExtinctionPointsCEAR[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsCEAR[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsCEAR]}]


(* ::Text:: *)
(*Combine plots and information for export*)


Column[{Style["Phototrophic state can also engage in predation\n", {FontSize -> 16, FontFamily -> "Arial", FontColor -> Black}],
Row[{Grid[
	(* Some extremely fancy rasterization of the plots without rasterizing their labels (removes the labels, raterizes, and re-adds the labels *)
	Map[
		Labeled[Rasterize[Show[#, FrameLabel -> None, PlotLabel -> None](*, RasterSize \[Rule] 500*), ImageSize -> Medium], 
		Append[(FrameLabel /. Options[#])[[All, 1]], Style[(PlotLabel /. Options[#]), 14, TextAlignment -> Center]], {Left, Bottom, Top},
		LabelStyle -> Append[(LabelStyle /. Options[#]), (TextAlignment -> Center)], RotateLabel -> True]&,
	{streamPlotsCE, streamPlotsCEAR}, {2}], 
Spacings -> {2, 3}], 

	(* Legend *)
	Column[{
		BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Column", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {25, 450}],
		PointLegend[{Black, Black}, {"Direction of selection", "Consumer extinction"}, 
			LegendMarkers -> {Graphics[{Arrowheads[0.5], Arrow[{{0, 0}, {1, 0}}]}], Graphics[Rectangle[{0, 0}]]},
			LegendMarkerSize -> {{75, 75}, {20, 20}}]},
	Frame -> True, Alignment -> {Center, Center}]}, 
Spacer[50]],
Graphics[{Thickness[0.0025], Arrowheads[0.025], Arrow[{{0, 0}, {100, 0}}], Inset[Style["Increasing light", 20], Background -> White]}, ImageSize -> {1200, 50}]},
Alignment -> Center]

Export["stream_plots_phototroph_eats.pdf", %]



