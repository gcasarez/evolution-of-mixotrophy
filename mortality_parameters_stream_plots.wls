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
gammaP = 2/3;
sigmaP = 0;
kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;
hW = 50;
hP = hP2 = 50;
lW = 0.5;
(*lH = lH2 = 0.25;
lP = lP2 = 0.25;*)
dMax = 0.3;
sigmaD = 0; (* No more decay trade-off *)
a = a2 = 0.15;
gammaA = 0; (* We want this to be 0 consistently, because we've decided the phototroph doesn't eat *)
c = 0.1;
sigmaC = 0; (* This doesn't do anything anymore *)


(* ::Section::Closed:: *)
(*Equations*)


Get[NotebookDirectory[] <> "model_and_functions.wls"]


(* ::Section:: *)
(*Phototroph mortality rate is higher*)


lH = lH2 = 0.15;
lP = lP2 = 0.35;


light = 50;


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARhighlP = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARhighlP = Flatten[ecoDynamicsCEARhighlP, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEARhighlP = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEARhighlP = Flatten[minPopSizesCEARhighlP, 1];


extinctionPointsCEARhighlP = Pick[ecoDynamicsCEARhighlP, Map[#[[4]] <= 10^-6&, minPopSizesCEARhighlP]];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEARhighlP = Pick[ecoDynamicsCEARhighlP, Map[#[[4]] > 10^-6&, minPopSizesCEARhighlP]];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEARhighlP = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEARhighlP];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and down*)


Select[derivsForPlotsCEARhighlP, #[[1]] == {1, 0}&][[All, 2]]


Export["higherPhototrophMortalityCEAR.wdx", derivsForPlotsCEARhighlP]


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEhighlP = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEhighlP = Flatten[ecoDynamicsCEhighlP, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEhighlP = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEhighlP = Flatten[minPopSizesCEhighlP, 1];


extinctionPointsCEhighlP = Pick[ecoDynamicsCEhighlP, Map[#[[4]] <= 10^-6&, minPopSizesCEhighlP]];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEhighlP = Pick[ecoDynamicsCEhighlP, Map[#[[4]] > 10^-6&, minPopSizesCEhighlP]];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEhighlP = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEhighlP];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and down*)


Select[derivsForPlotsCEhighlP, #[[1]] == {1, 0}&][[All, 2]]


Export["higherPhototrophMortalityCE.wdx", derivsForPlotsCEhighlP]


(* ::Section:: *)
(*Heterotroph mortality rate is higher*)


lH = lH2 = 0.35;
lP = lP2 = 0.15;


(*light = 50;*)


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARhighlH = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARhighlH = Flatten[ecoDynamicsCEARhighlH, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEARhighlH = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEARhighlH = Flatten[minPopSizesCEARhighlH, 1];


extinctionPointsCEARhighlH = Pick[ecoDynamicsCEARhighlH, Map[#[[4]] <= 10^-6&, minPopSizesCEARhighlH]];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEARhighlH = Pick[ecoDynamicsCEARhighlH, Map[#[[4]] > 10^-6&, minPopSizesCEARhighlH]];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEARhighlH = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEARhighlH];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and up.*)


Select[derivsForPlotsCEARhighlH, #[[1]] == {1, 0}&][[All, 2]]


Export["higherHeterotrophMortalityCEAR.wdx", derivsForPlotsCEARhighlH]


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEhighlH = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEhighlH = Flatten[ecoDynamicsCEhighlH, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEhighlH = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEhighlH = Flatten[minPopSizesCEhighlH, 1];


extinctionPointsCEhighlH = Pick[ecoDynamicsCEhighlH, Map[#[[4]] <= 10^-6&, minPopSizesCEhighlH]];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEhighlH = Pick[ecoDynamicsCEhighlH, Map[#[[4]] > 10^-6&, minPopSizesCEhighlH]];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEhighlH = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEhighlH];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and up*)


Select[derivsForPlotsCEhighlH, #[[1]] == {1, 0}&][[All, 2]]


Export["higherHeterotrophMortalityCE.wdx", derivsForPlotsCEhighlH]


(* ::Section:: *)
(*Default mortality rates (lH = lP)*)


lH = lH2 = 0.25;
lP = lP2 = 0.25;


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARdefault = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARdefault = Flatten[ecoDynamicsCEARdefault, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEARdefault = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEARdefault = Flatten[minPopSizesCEARdefault, 1];


extinctionPointsCEARdefault = Pick[ecoDynamicsCEARdefault, Map[#[[4]] <= 10^-6&, minPopSizesCEARdefault]];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEARdefault = Pick[ecoDynamicsCEARdefault, Map[#[[4]] > 10^-6&, minPopSizesCEARdefault]];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEARdefault = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEARdefault];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and down.*)


Select[derivsForPlotsCEARdefault, #[[1]] == {1, 0}&][[All, 2]]


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEdefault = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEdefault = Flatten[ecoDynamicsCEdefault, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEdefault = Table[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEdefault = Flatten[minPopSizesCEdefault, 1];


extinctionPointsCEdefault = Pick[ecoDynamicsCEdefault, Map[#[[4]] <= 10^-6&, minPopSizesCEdefault]];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEdefault = Pick[ecoDynamicsCEdefault, Map[#[[4]] > 10^-6&, minPopSizesCEdefault]];


(* ::Subsubsection:: *)
(*Find the derivative of invasion fitness (with respect to retention and replication) at each (retention, replication) combination where the mixotroph does not go extinct*)


derivsForPlotsCEdefault = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
	conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEdefault];


(* ::Text:: *)
(*Check the derivatives of invasion fitness at obligately kleptoplastidic mixotrophy. They point right and up*)


Select[derivsForPlotsCEdefault, #[[1]] == {1, 0}&][[All, 2]]


(* ::Section:: *)
(*Plots*)


(* ::Subsection:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions -> {"OperatorSubstitution"->False}}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions->{"OperatorSubstitution"->False}}];


(* ::Subsection:: *)
(*Color scheme (prey density) and legend*)


maxPrey = Max[{noExtinctionPointsCEhighlP[[All, 5]], noExtinctionPointsCEdefault[[All, 5]], noExtinctionPointsCEhighlH[[All, 5]],
noExtinctionPointsCEARhighlP[[All, 5]], noExtinctionPointsCEARdefault[[All, 5]], noExtinctionPointsCEARhighlH[[All, 5]]}];


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


(*BarLegend[{preyDensityColors[#]&, {0, maxPrey}}]*)
preyDensityLegend = BarLegend[{preyDensityColors[Exp[#] - 1]&, {Log[1], Log[1 + maxPrey]}}, 
Ticks -> Table[{Log[i + 1], ToString[i]}, {i, {0, 1, 2.5, 5, 10, 15, 25, 40}}],
LegendLayout -> "Column", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {25, 450}];


(* ::Subsection:: *)
(*Set up plot labels*)


plotLabels = (*{"(A) High phototroph mortality\nLow heterotroph mortality\nNo tradeoff", 
"(B) Equal heterotroph and\nphototroph mortality\nNo tradeoff",
"(C) Low phototroph mortality\nHigh heterotroph mortality\nNo tradeoff",
"(D) High phototroph mortality\nLow heterotroph mortality\nAttack tradeoff", 
"(E) Equal heterotroph and\nphototroph mortality\nAttack tradeoff",
"(F) Low phototroph mortality\nHigh heterotroph mortality\nAttack tradeoff"};*)
Table[Pane[Style[letter, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}], ImageSize -> 300],
{letter, {"(A)", "(B)", "(C)", "(D)", "(E)", "(F)"}}];


(* ::Subsection:: *)
(*Mark points where the population went extinct with black boxes*)


extinctionGraphics = 
Table[Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPoints[[All, 2;;3]]]], 
	{extinctionPoints, 
	{extinctionPointsCEhighlP, extinctionPointsCEdefault, extinctionPointsCEhighlH, 
	extinctionPointsCEARhighlP, extinctionPointsCEARdefault, extinctionPointsCEARhighlH}}];


(* ::Subsection:: *)
(*Collect derivatives in a list*)


derivsForPlots = {derivsForPlotsCEhighlP, derivsForPlotsCEdefault, derivsForPlotsCEhighlH,
	derivsForPlotsCEARhighlP, derivsForPlotsCEARdefault, derivsForPlotsCEARhighlH};


(* ::Subsection:: *)
(*Make plots*)


plots = Table[
Show[{ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlots[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	FrameLabel -> {"Retention Probability (\[Phi])", "Replication Probability (\[Eta])"},
	PlotLabel -> plotLabels[[i]], ImageSize -> Medium], 
	
	ListStreamPlot[Map[{#[[1]], #[[2]]}&, derivsForPlots[[i]]], StreamStyle -> Black, StreamColorFunction -> None],
	
	extinctionGraphics[[i]]}],
{i, Length[derivsForPlots]}]


(* ::Text:: *)
(*Combine plots and information for export*)


Row[{Grid[

	(* Row labels *)
	Join[Map[Rotate[Style[#, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}], 90Degree]&, 
		{{""}, {"No tradeoff"}, {"Attack tradeoff"}}, {2}],
		
	(* Column labels *)
	Join[Map[Style[#, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}]&, 
		{{"High phototroph mortality, low heterotroph mortality", "Equal phototroph and heterotroph mortality", "Low phototroph mortality, high heterotroph mortality"}},
		{2}],
	 
	(* Some extremely fancy rasterization of the plots without rasterizing their labels (removes the labels, raterizes, and re-adds the labels *)
	Map[
		Labeled[Rasterize[Show[#, FrameLabel -> None, PlotLabel -> None, ImagePadding -> 1](*, RasterSize \[Rule] 500*), ImageSize -> Medium], 
		Append[(FrameLabel /. Options[#])[[All, 1]], Style[(PlotLabel /. Options[#]), 14, TextAlignment -> Center]], {Left, Bottom, Top},
		LabelStyle -> Append[(LabelStyle /. Options[#]), (TextAlignment -> Center)], RotateLabel -> True]&,
	{plots[[1;;3]], plots[[4;;6]]}, {2}]
		(* close out the two joinings of row and column labels *)], 2], 
Spacings -> {2, 3}], 

	(* Legend *)
	Column[{
		BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Column", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {25, 450}],
		PointLegend[{Black, Black}, {"Direction of selection", "Consumer extinction"}, 
			LegendMarkers -> {Graphics[{Arrowheads[0.5], Arrow[{{0, 0}, {1, 0}}]}], Graphics[Rectangle[{0, 0}]]},
			LegendMarkerSize -> {{75, 75}, {20, 20}}]},
	Frame -> True, Alignment -> {Center, Center}]}, 
Spacer[50]]

Export["stream_plots_varying_mortality.pdf", %]



