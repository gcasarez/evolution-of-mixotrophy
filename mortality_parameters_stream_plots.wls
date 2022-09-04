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


ecoDynamicsCEARhighlP = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARhighlP = Flatten[ecoDynamicsCEARhighlP, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEARhighlP = ParallelTable[{light, retention, replication, 
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


Export["higherPhototrophMortalityCEAR.wdx", derivsForPlotsCEARhighlP];


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEhighlP = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEhighlP = Flatten[ecoDynamicsCEhighlP, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEhighlP = ParallelTable[{light, retention, replication, 
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


Export["higherPhototrophMortalityCE.wdx", derivsForPlotsCEhighlP];


(* ::Section:: *)
(*Heterotroph mortality rate is higher*)


lH = lH2 = 0.35;
lP = lP2 = 0.15;


(*light = 50;*)


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARhighlH = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARhighlH = Flatten[ecoDynamicsCEARhighlH, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEARhighlH = ParallelTable[{light, retention, replication, 
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


Export["higherHeterotrophMortalityCEAR.wdx", derivsForPlotsCEARhighlH];


(* ::Subsection:: *)
(*Conversion trade-off only*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEhighlH = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEhighlH = Flatten[ecoDynamicsCEhighlH, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEhighlH = ParallelTable[{light, retention, replication, 
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


Export["higherHeterotrophMortalityCE.wdx", derivsForPlotsCEhighlH];


(* ::Section:: *)
(*Default mortality rates (lH = lP)*)


lH = lH2 = 0.25;
lP = lP2 = 0.25;


(* ::Subsection:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEARdefault = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6,
	conversionTradeOff -> True, attackTradeOff -> True]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEARdefault = Flatten[ecoDynamicsCEARdefault, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEARdefault = ParallelTable[{light, retention, replication, 
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


ecoDynamicsCEdefault = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]],
 {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEdefault = Flatten[ecoDynamicsCEdefault, {{1,2}, {3}}];


(* ::Subsubsection:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEdefault = ParallelTable[{light, retention, replication, 
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


(* ::Subsection::Closed:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions -> {"OperatorSubstitution"->False}}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions->{"OperatorSubstitution"->False}}];


(* ::Subsection::Closed:: *)
(*Color scheme (prey density) and legend*)


maxPrey = Max[{noExtinctionPointsCEhighlP[[All, 5]], noExtinctionPointsCEdefault[[All, 5]], noExtinctionPointsCEhighlH[[All, 5]],
noExtinctionPointsCEARhighlP[[All, 5]], noExtinctionPointsCEARdefault[[All, 5]], noExtinctionPointsCEARhighlH[[All, 5]]}];


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


(*BarLegend[{preyDensityColors[#]&, {0, maxPrey}}]*)
preyDensityLegend = BarLegend[{preyDensityColors[Exp[#] - 1]&, {Log[1], Log[1 + maxPrey]}}, 
Ticks -> Table[{Log[i + 1], ToString[i]}, {i, {0, 1, 2.5, 5, 10, 15, 25, 40}}],
LegendLayout -> "Column", LegendLabel -> "Prey density (cells/mL)", LegendMarkerSize -> {25, 450}];


(* ::Subsection::Closed:: *)
(*Set up plot labels*)


plotLabels = 
Table[Pane[Style[letter, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}], ImageSize -> 72 * 2.25],
{letter, {"(A)", "(B)", "(C)", "(D)", "(E)", "(F)"}}];


(* ::Subsection::Closed:: *)
(*Mark points where the population went extinct with black boxes*)


extinctionGraphics = 
Table[Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, extinctionPoints[[All, 2;;3]]]], 
	{extinctionPoints, 
	{extinctionPointsCEhighlP, extinctionPointsCEdefault, extinctionPointsCEhighlH, 
	extinctionPointsCEARhighlP, extinctionPointsCEARdefault, extinctionPointsCEARhighlH}}];


(* ::Subsection::Closed:: *)
(*Collect points where the population didn't go extinct*)


noExtinctionPoints = {noExtinctionPointsCEhighlP, noExtinctionPointsCEdefault, noExtinctionPointsCEhighlH, 
	noExtinctionPointsCEARhighlP, noExtinctionPointsCEARdefault, noExtinctionPointsCEARhighlH};


(* ::Subsection::Closed:: *)
(*Collect derivatives in a list*)


derivsForPlots = {derivsForPlotsCEhighlP, derivsForPlotsCEdefault, derivsForPlotsCEhighlH,
	derivsForPlotsCEARhighlP, derivsForPlotsCEARdefault, derivsForPlotsCEARhighlH};


(* ::Subsection::Closed:: *)
(*Make plots*)


(* ::Subsubsection::Closed:: *)
(*Prey density*)


preyDensityPlots = Table[Show[Graphics[Inset[
	Rasterize[ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlots[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1, 
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	PlotRangePadding -> 0, ImagePadding -> 0, ImageMargins -> 0, FrameLabel -> None, Frame -> None]],
	
	{0, 0}, {0, 0}, Scaled[1/1.1]], PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.05, ImagePadding -> Full], 
	Frame -> True, FrameStyle -> Black, LabelStyle -> {Black, "Arial", 12}, ImageSize -> 72*2.25],

{i, Length[derivsForPlots]}]


(* ::Subsubsection::Closed:: *)
(*Selection gradient*)


vectorPlots = Table[ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivsForPlots[[i]]], 
	VectorStyle -> Darker[Gray], VectorColorFunction -> None, VectorPoints -> noExtinctionPoints[[i, All, 2;;3]],
	VectorScale -> {Automatic, Automatic, None}],
{i, Length[derivsForPlots]}]


(* ::Subsubsection::Closed:: *)
(*Combine plots with extinction graphics*)


plots = MapThread[Show, {preyDensityPlots, vectorPlots, extinctionGraphics}];


(* ::Subsection:: *)
(*Combine plots and information for export*)


(* ::Subsubsection::Closed:: *)
(*Axis labels*)


ylabel = Rotate[Pane[Style["Replication (\[Eta])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12]], Pi/2]
xlabel = Pane[Style["Retention (\[Phi])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], 
	ImageSize -> 2.5 * 72, Alignment -> {0.2, Center}]


(* ::Text:: *)
(*Row & column labels*)


rowLabels = Map[Rotate[Style[#, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}], 90Degree]&, 
		{"No tradeoff", "Attack tradeoff"}];


columnLabels = Map[Style[#, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}]&, 
		{"", "",  "High phototroph mortality,\nlow heterotroph mortality", "Equal phototroph and\nheterotroph mortality", 
		"Low phototroph mortality,\nhigh heterotroph mortality"}];


(* ::Subsubsection:: *)
(*Legend*)


(* ::Text:: *)
(*arrowhead for the legend*)


myArrowhead = Graphics[Triangle[{{0, 0}, {-1, -1/4}, {-1, 1/4}}]];


legend = Framed[Column[{
	Pane[Style["Legend", FontSize -> 12, FontFamily -> "Arial", FontWeight -> Bold, FontColor -> Black], 
		ImageSize -> {72 * 1, 30}, Alignment -> Center],
		
	BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Column", 
	"TickSide" -> Left,
	LegendLabel -> Placed[Pane["Prey density (cells/\!\(\*SuperscriptBox[\(m\), \(2\)]\))", ImageSize -> 72 * 1], Top],
	LegendMarkerSize -> {24, 144}, LabelStyle -> {Black, "Arial", 12}],
	
	PointLegend[{Darker[Gray], Black, Black, Black}, 
	{Pane["Direction of selection", ImageSize -> 72 * 1],
		Pane["Consumer extinction", ImageSize -> 72 * 1]}, 
	LegendMarkers -> {Graphics[{Arrowheads[{{0.35, 1, myArrowhead}}], Arrow[{{0, 0}, {16, 0}}]}, AspectRatio -> 1/8], 
		Graphics[Rectangle[{0, 0}]]},
	LegendMarkerSize -> {{32, 16}, 16},
	LabelStyle -> {Black, "Arial", 12}]},
			
	Alignment -> {Center, Center}, Spacings -> 0],
RoundingRadius -> 10]


(* ::Subsubsection:: *)
(*Plot compilation*)


compiledPlots = Grid[{{Grid[ (* We do a grid of a grid to get around an export glitch *)
	{columnLabels,
	Join[{rowLabels[[1]], ylabel}, plots[[1;;3]]], 
	Join[{rowLabels[[2]], ylabel}, plots[[4;;6]]],
	{"", "", xlabel, xlabel, xlabel}},
	Alignment -> {Center, Center},
	Spacings -> {{0 (* leftmost space, before grid starts *), Scaled[1/100] (* between row labels and ylabels *),
			0 (* between ylabels and plots *), Scaled[1/150] (* between plots *), Scaled[1/150] (* between plots *)}, 
		{0 (* uppermost space, before grid starts *), Scaled[1/100] (* between column labels and plots *),
		Scaled[5/144] (* between plots *), 
		Scaled[1/150] (* between xlabels and plots *)}}],
		
	legend}}, Spacings -> 1]

Export["vector_plots_varying_mortality.pdf", %]



