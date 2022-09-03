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


(* ::Subsection::Closed:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection::Closed:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEAR = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEAR = Flatten[ecoDynamicsCEAR, {{1}, {2,3}, {4}}];


(* ::Subsubsection::Closed:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCEAR = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEAR = Flatten[minPopSizesCEAR, {{1}, {2, 3}}];


extinctionPointsCEAR = Table[Pick[ecoDynamicsCEAR[[i]], Map[#[[4]] <= 10^-6&, minPopSizesCEAR[[i]]]], {i, Length[ecoDynamicsCEAR]}];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCEAR = Table[Pick[ecoDynamicsCEAR[[i]], Map[#[[4]] > 10^-6&, minPopSizesCEAR[[i]]]], {i, Length[ecoDynamicsCEAR]}];


(* ::Subsubsection::Closed:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCEAR = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEAR, {2}];


Export["phototroph_eats_CEAR.wdx", derivsForPlotsCEAR]


(* ::Subsection::Closed:: *)
(*Conversion trade-off only*)


(* ::Subsubsection::Closed:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCE = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
conversionTradeOff -> True, attackTradeOff -> False]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCE = Flatten[ecoDynamicsCE, {{1}, {2,3}, {4}}];


(* ::Subsubsection::Closed:: *)
(*What trait combinations cause the mixotroph population to go extinct, or drop below 10^-6?*)


minPopSizesCE = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCE = Flatten[minPopSizesCE, {{1}, {2, 3}}];


extinctionPointsCE = Table[Pick[ecoDynamicsCE[[i]], Map[#[[4]] <= 10^-6&, minPopSizesCE[[i]]]], {i, Length[ecoDynamicsCE]}];


(* ::Text:: *)
(*Also find the points where the mixotroph population does not go extinct*)


noExtinctionPointsCE = Table[Pick[ecoDynamicsCE[[i]], Map[#[[4]] > 10^-6&, minPopSizesCE[[i]]]], {i, Length[ecoDynamicsCE]}];


(* ::Subsubsection::Closed:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCE = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCE, {2}];


Export["phototroph_eats_CE.wdx", derivsForPlotsCE]


(* ::Section:: *)
(*Plots*)


(* ::Subsection::Closed:: *)
(*Set font and colors*)


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


(* ::Subsection::Closed:: *)
(*Make plots*)


(* ::Subsubsection::Closed:: *)
(*Prey density*)


{preyDensityPlotsCE, preyDensityPlotsCEAR} = Table[Show[Graphics[Inset[
	Rasterize[ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlots[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1, 
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	PlotRangePadding -> 0, ImagePadding -> 0, ImageMargins -> 0, FrameLabel -> None, Frame -> None]],
	
	{0, 0}, {0, 0}, Scaled[1/1.1]], PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.05, ImagePadding -> Full], 
	Frame -> True, FrameStyle -> Black, LabelStyle -> {Black, "Arial", 12}, ImageSize -> 72*2.25],

{derivsForPlots, {derivsForPlotsCE, derivsForPlotsCEAR}}, {i, Length[lightLevels]}]


(* ::Subsubsection::Closed:: *)
(*Selection gradient*)


vectorPlotsCE = Table[ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCE[[i]]], 
	VectorStyle -> Darker[Gray], VectorColorFunction -> None, VectorPoints -> noExtinctionPointsCE[[i, All, 2;;3]],
	VectorScale -> {Automatic, Automatic, None}],
{i, Length[lightLevels]}]

vectorPlotsCEAR = Table[ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCEAR[[i]]], 
	VectorStyle -> Darker[Gray], VectorColorFunction -> None, VectorPoints -> noExtinctionPointsCEAR[[i, All, 2;;3]],
	VectorScale -> {Automatic, Automatic, None}],
{i, Length[lightLevels]}]


(* ::Subsubsection::Closed:: *)
(*Extinction points*)


{extinctionGraphicsCE, extinctionGraphicsCEAR} = 
Table[Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, 
	extinctionPoints[[i, All, 2;;3]]]],
	{extinctionPoints, {extinctionPointsCE, extinctionPointsCEAR}}, {i, Length[lightLevels]}];


(* ::Subsubsection::Closed:: *)
(*Combine plots with extinction graphics*)


plotsCE = MapThread[Show, {preyDensityPlotsCE, vectorPlotsCE, extinctionGraphicsCE}];
plotsCEAR = MapThread[Show, {preyDensityPlotsCEAR, vectorPlotsCEAR, extinctionGraphicsCEAR}];


(* ::Text:: *)
(*Label combined plots*)


plotLabels = Map[Pane[Style[#, FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], ImageSize -> 160, Alignment -> Left]&,
	{"(A)", "(B)", "(C)", "(D)", "(E)", "(F)"}];


plotsCE = MapThread[Labeled[#1, #2, {Top, Left}]&, {plotsCE, plotLabels[[1;;3]]}];
plotsCEAR = MapThread[Labeled[#1, #2, {Top, Left}]&, {plotsCEAR, plotLabels[[4;;6]]}];


(* ::Subsubsection::Closed:: *)
(*Labels for plots*)


ylabel = Rotate[Pane[Style["Replication (\[Eta])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12]], Pi/2];
xlabel = Pane[Style["Retention (\[Phi])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], 
	ImageSize -> 2.5 * 72, Alignment -> {0.2, Center}];


(* ::Text:: *)
(*Row labels*)


rowLabels = Map[Rotate[Style[#, {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}], 90Degree]&, 
		{"No tradeoff", "Attack tradeoff"}];


(* ::Text:: *)
(*Arrow showing direction of increasing light*)


lightLabel = Graphics[{Thickness[0.0035], Arrowheads[0.025], Line[{{0, 0}, {200, 0}}], Arrow[{{344, 0}, {544, 0}}], 
Inset[Pane[Style["Increasing light", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], Alignment -> Center]]}];


lightLabel = Pane[Show[lightLabel, ImageSize -> 500], Alignment -> 0.38, ImageSize -> 550];


(* ::Subsubsection::Closed:: *)
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
RoundingRadius -> 10];


(* ::Subsubsection::Closed:: *)
(*Plot compilation*)


compiledPlots = Grid[{{Grid[ (* We do a grid of a grid to get around an export glitch *)
	{Join[{rowLabels[[1]], ylabel}, plotsCE], 
	Join[{rowLabels[[2]], ylabel}, plotsCEAR],
	{"", "", xlabel, xlabel, xlabel},
	{"", "",  lightLabel, SpanFromLeft, SpanFromLeft}},
	Alignment -> {Center, Center},
	Spacings -> {{0 (* leftmost space, before grid starts *), Scaled[1/100] (* between row labels and ylabels *),
			0 (* between ylabels and plots *), Scaled[1/150] (* between plots *), Scaled[1/150] (* between plots *)}, 
		{0 (* uppermost space, before grid starts *), Scaled[5/144] (* between plots *), 
		Scaled[1/150] (* between xlabels and plots *), Scaled[1/150] (* between arrow and xlabels *)}}],
		
	legend}}, Spacings -> 1]

Export["vector_plots_phototroph_eats.pdf", %]



