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


(* ::Section:: *)
(*Equations*)


Get[NotebookDirectory[] <> "model_and_functions.wls"]


(* ::Section:: *)
(*Find ecological equilibrium & derivative of invasion fitness*)


lightLevels = {25, 35, 50};


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEAR = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEAR = Flatten[ecoDynamicsCEAR, {{1}, {2,3}, {4}}];


(* ::Subsubsection:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


minPopSizesCEAR = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEAR = Flatten[minPopSizesCEAR, {{1}, {2, 3}}];


extinctionPointsCEAR = Map[Select[#[[4]] <= 10^-6&], minPopSizesCEAR];


(* ::Text:: *)
(*Where does the population not go extinct?*)


noExtinctionPointsCEAR = Pick[ecoDynamicsCEAR, Map[#[[4]] > 10^-6&, minPopSizesCEAR, {2}]];


(* ::Subsubsection:: *)
(*Plot the points where the population goes extinct (green dots, overlaid on colors showing mixotroph population sizes)*)


Table[Show[{
(* Plot of the mixotroph population sizes *)
ListDensityPlot[Map[{#[[2]], #[[3]], #[[6]] + #[[7]]}&, ecoDynamicsCEAR[[i]]], 
FrameLabel -> {"Retention", "Replication"}, PlotRangePadding -> 0.02,
PlotLabel -> "(" <> ToUpperCase[Alphabet[]][[i]] <> ") Mixotroph population (cells/mL)\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", 
ImageSize -> 400, PlotLegends -> Automatic, PlotRange -> All, ClippingStyle -> Black, ColorFunction -> "ValentineTones"],

(* Plot points where the mixotroph goes extinct *)
ListPlot[extinctionPointsCEAR[[i, All, {2, 3}]], PlotStyle -> {Green, PointSize[Large]}]}],
{i, Length[ecoDynamicsCEAR]}]


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness at the points where the resident is not extinct*)


derivsForPlotsCEAR = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEAR, {2}];


(* ::Text:: *)
(*Quick vector plot to get an idea of the magnitude/direction of derivative of invasion fitness. Colors indicate magnitude; arrows indicate direction*)


Table[
Show[ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsCEAR[[i]]], 
	VectorPoints -> noExtinctionPointsCEAR[[i, All, 2;;3]],
	ColorFunction -> "Pastel", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium, VectorScale -> {Automatic, Automatic, None}], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCEAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCEAR]}]


(* ::Subsubsection:: *)
(*Now we'll make some stream plots.*)


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Pastel"(*"Aquamarine"*)][Sqrt[(p/maxPrey)](**7/10 + 3/10*)]


Table[
Show[ListStreamDensityPlot[Map[{#[[1]], {#[[2]], #[[3,1]]}}&, derivsForPlotsCEAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full},
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False, 
	StreamStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCEAR[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCEAR]}]


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCE = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCE = Flatten[ecoDynamicsCE, {{1}, {2,3}, {4}}];


(* ::Subsubsection:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


minPopSizesCE = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCE = Flatten[minPopSizesCE, {{1}, {2, 3}}];


extinctionPointsCE = Map[Select[#[[4]] <= 10^-6&], minPopSizesCE];


(* ::Text:: *)
(*Where does the population not go extinct?*)


noExtinctionPointsCE = Pick[ecoDynamicsCE, Map[#[[4]] > 10^-6&, minPopSizesCE, {2}]];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCE = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCE, {2}];


(* ::Text:: *)
(*Quick vector plot to get an idea of the magnitude/direction of derivative of invasion fitness. Colors indicate magnitude; arrows indicate direction*)


Table[
Show[ListVectorDensityPlot[Map[{#[[1]], {#[[2]], Norm[#[[2]]]}}&, derivsForPlotsCE[[i]]], 
	VectorPoints -> noExtinctionPointsCE[[i, All, 2;;3]],
	ColorFunction -> "Pastel", VectorStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium, VectorScale -> {Automatic, Automatic, None}], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCE]}]


(* ::Subsubsection:: *)
(*Make some stream plots.*)


maxPrey = Max[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]];
preyDensityColors[p_] := ColorData["Pastel"][Sqrt[p/maxPrey]]


Table[
Show[ListStreamDensityPlot[Map[{#[[1]], {#[[2]], #[[3,1]]}}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 0,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False, 
	StreamStyle -> Black, PlotLegends -> Automatic, FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") Surface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/20, 0], Max[#[[2]] - 1/20, 0]}, {Min[#[[1]] + 1/20, 1], Min[#[[2]] + 1/20, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]],
{i, Length[derivsForPlotsCE]}]


(* ::Section:: *)
(*Put the plots together*)


(* ::Text:: *)
(*Color scheme (prey density)*)


Histogram[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]]]


Histogram[Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]]


maxPrey = Max[{Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]],
	Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]}];


ClearAll[preyDensityColors]


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


(* ::Text:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}];
(*SetOptions[SwatchLegend, LabelStyle \[Rule] {FontSize -> 12, FontFamily -> "Arial"}];
SetOptions[PointLegend, LabelStyle \[Rule] {FontSize -> 12, FontFamily -> "Arial"}];*)


streamPlotsCE = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i]]] <> ") No tradeoff\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCE[[i]]], StreamStyle -> Black(*, StreamPoints \[Rule] noExtinctionPointsCE[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsCE[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsCE]}]

streamPlotsCEAR = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1,2]], #[[3, 1]]}&, derivsForPlotsCEAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> "(" <> Capitalize[Alphabet[][[i + Length[lightLevels]]]] <> ") Attack tradeoff\nSurface light = " <> ToString[lightLevels[[i]]] <> " \[Mu]mol photons/m^2/s", ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCEAR[[i]]], StreamStyle -> Black(*, StreamPoints \[Rule] noExtinctionPointsCEAR[[i, All, 1]]*)],
	
	Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPointsCEAR[[i, All, 2;;3]]]]}],
{i, Length[derivsForPlotsCEAR]}]


Column[{Row[{Grid[
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
Graphics[{Thickness[0.0025], Arrowheads[0.025], Arrow[{{0, 0}, {100, 0}}], Inset[Style["Increasing light", 20], Background -> White]}, ImageSize -> {1200, 50}]}]

Export["stream_plots.pdf", %]



