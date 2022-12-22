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
(*Find ecological equilibrium & derivative of invasion fitness*)


lightLevels = {25, 35, 50};


(* ::Subsection::Closed:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection::Closed:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCEAR = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> True]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCEAR = Flatten[ecoDynamicsCEAR, {{1}, {2,3}, {4}}];


(* ::Subsubsection::Closed:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


minPopSizesCEAR = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> True]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCEAR = Flatten[minPopSizesCEAR, {{1}, {2, 3}}];


extinctionPointsCEAR = Map[Select[#[[4]] <= 10^-6&], minPopSizesCEAR];


(* ::Text:: *)
(*Where does the population not go extinct?*)


noExtinctionPointsCEAR = Pick[ecoDynamicsCEAR, Map[#[[4]] > 10^-6&, minPopSizesCEAR, {2}]];


(* ::Subsubsection::Closed:: *)
(*Find the derivative of invasion fitness at the points where the resident is not extinct*)


derivsForPlotsCEAR = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> True],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCEAR, {2}];


(* ::Subsection::Closed:: *)
(*Conversion trade-off only*)


(* ::Subsubsection::Closed:: *)
(*Find the average prey, heterotroph, and phototroph densities across a range of retention and replication levels*)


ecoDynamicsCE = ParallelTable[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, conversionTradeOff -> True, attackTradeOff -> False]],
{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


ecoDynamicsCE = Flatten[ecoDynamicsCE, {{1}, {2,3}, {4}}];


(* ::Subsubsection::Closed:: *)
(*What light/trait combinations cause the mixotroph population to go extinct?*)


minPopSizesCE = ParallelTable[{light, retention, replication, 
	minimumMixotroph[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
	conversionTradeOff -> True, attackTradeOff -> False]},
 {light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];


minPopSizesCE = Flatten[minPopSizesCE, {{1}, {2, 3}}];


extinctionPointsCE = Map[Select[#[[4]] <= 10^-6&], minPopSizesCE];


(* ::Text:: *)
(*Where does the population not go extinct?*)


noExtinctionPointsCE = Pick[ecoDynamicsCE, Map[#[[4]] > 10^-6&, minPopSizesCE, {2}]];


(* ::Subsubsection::Closed:: *)
(*Find the derivative of invasion fitness*)


derivsForPlotsCE = 
ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, conversionTradeOff -> True, attackTradeOff -> False],
{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPointsCE, {2}];


(* ::Section::Closed:: *)
(*Plot prey density and selection gradient*)


(* ::Subsection::Closed:: *)
(*Colors and font*)


(* ::Subsubsection::Closed:: *)
(*Color scheme (prey density)*)


maxPrey = Max[{Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCEAR][[All, All, 5]]],
	Flatten[Map[Select[#[[4]] == False&], ecoDynamicsCE][[All, All, 5]]]}];


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


(* ::Subsubsection::Closed:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black}];


(* ::Subsection::Closed:: *)
(*Prey density plots*)


(* ::Subsubsection::Closed:: *)
(*No attack trade-off*)


preyDensityPlotsCE = Table[Show[Graphics[Inset[
	Rasterize[ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlotsCE[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1, 
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	PlotRangePadding -> 0, ImagePadding -> 0, ImageMargins -> 0, FrameLabel -> None, Frame -> None]],
	
	{0, 0}, {0, 0}, Scaled[1/1.1]], PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.05, ImagePadding -> Full], 
	Frame -> True, FrameStyle -> Black, LabelStyle -> {Black, "Arial", 12}, ImageSize -> 72*2.25],

{i, Length[derivsForPlotsCE]}]


(* ::Subsubsection::Closed:: *)
(*Attack trade-off*)


preyDensityPlotsCEAR = Table[Show[Graphics[Inset[
	Rasterize[ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivsForPlotsCEAR[[i]]], 
	PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1, 
	ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
	PlotRangePadding -> 0, ImagePadding -> 0, ImageMargins -> 0, FrameLabel -> None, Frame -> None]],
	
	{0, 0}, {0, 0}, Scaled[1/1.1]], PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.05, ImagePadding -> Full], 
	Frame -> True, FrameStyle -> Black, LabelStyle -> {Black, "Arial", 12}, ImageSize -> 72*2.25],

{i, Length[derivsForPlotsCEAR]}]


(* ::Subsection::Closed:: *)
(*Selection vector plots*)


(* ::Subsubsection::Closed:: *)
(*No attack trade-off*)


vectorPlotsCE = Table[ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCE[[i]]], 
	VectorStyle -> Darker[Gray], VectorColorFunction -> None, VectorPoints -> noExtinctionPointsCE[[i, All, 2;;3]],
	VectorScale -> {Automatic, Automatic, None}],
{i, Length[derivsForPlotsCE]}]


(* ::Subsubsection::Closed:: *)
(*Attack trade-off*)


vectorPlotsCEAR = Table[ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCEAR[[i]]], 
	VectorStyle -> Darker[Gray], VectorColorFunction -> None, VectorPoints -> noExtinctionPointsCEAR[[i, All, 2;;3]],
	VectorScale -> {Automatic, Automatic, None}],
{i, Length[derivsForPlotsCEAR]}]


(* ::Subsection::Closed:: *)
(*Plots showing extinction points*)


(* ::Subsubsection::Closed:: *)
(*No attack trade-off*)


extinctionPlotCE = Table[
Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, 
extinctionPointsCE[[i, All, 2;;3]]]],

{i, Length[derivsForPlotsCE]}]


(* ::Subsubsection::Closed:: *)
(*Attack trade-off*)


extinctionPlotCEAR = Table[
Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, 
extinctionPointsCEAR[[i, All, 2;;3]]]],

{i, Length[derivsForPlotsCEAR]}]


(* ::Subsection::Closed:: *)
(*Markers for the ESSs*)


(* ::Subsubsection::Closed:: *)
(*No attack trade-off*)


(* ::Text:: *)
(*Find the retention probabilities where replication = 1 is an ESS*)


fkpESSpointsCE = Table[
	Select[derivsForPlotsCE[[i]], (#[[1, 2]] == 1) (* replication = 1*) && 
		(Chop[#[[2, 1]], 10^-1] == 0) (* no selection on retention *) && 
		(#[[2, 2]] > 0) (* positive selection on replication *)&][[All, 1, 1]],
	{i, Length[derivsForPlotsCE]}];


(* ::Text:: *)
(*Get the minimum & maximum retention probabilities where replication = 1 is an ESS*)


fkpESSrangeCE = Map[MinMax, fkpESSpointsCE];


(* ::Text:: *)
(*Draw rectangles marking the region where relication = 1 is an ESS*)


fkpRegionsCE = Map[Graphics[{FaceForm[], EdgeForm[{Thickness[0.03], Black}], 
	Rectangle[{Max[#[[1]] - 1/20, 0], 1.015}, {Min[#[[2]] + 1/20, 1], 1.085}, RoundingRadius -> 0.035]}]&, fkpESSrangeCE];


(* ::Subsubsection::Closed:: *)
(*Attack trade-off, obligately kleptoplastidic mixotrophy*)


(* ::Text:: *)
(*Marker for ESS where retention > 0, replication = 0*)


(* ::Text:: *)
(*Use interpolation to find a function for d fitness / d retention when replication = 0*)


dFitdPhiInterpolation = Table[Interpolation[Select[derivsForPlotsCEAR[[2]], #[[1, 2]] == 0&][[All, 1;;2, 1]]],
	{i, Length[derivsForPlotsCEAR]}];


(* ::Text:: *)
(*ESS is at the point where d fitness / d retention = 0*)


okpPhiValue = Map[phi /. FindRoot[#[phi], {phi, 0.9}]&, dFitdPhiInterpolation];


(* ::Text:: *)
(*Draw a circle at the ESS*)


okpPoint = Map[Graphics[{Thickness[0.03], Black, Circle[{#, 0}, 0.07]}]&, okpPhiValue];


(* ::Subsubsection::Closed:: *)
(*Attack tradeoff, facultatively kleptoplastidic mixotrophy*)


(* ::Text:: *)
(*Find the retention probabilities where replication = 1 is an ESS*)


fkpESSpointsCEAR = Table[
	Select[derivsForPlotsCEAR[[i]], (#[[1, 2]] == 1) (* replication = 1*) && 
		(Chop[#[[2, 1]], 10^-1] == 0) (* no selection on retention *) && 
		(#[[2, 2]] > 0) (* positive selection on replication *)&][[All, 1, 1]],
	{i, Length[derivsForPlotsCEAR]}];


(* ::Text:: *)
(*Find the points where replication = 1 and the consumer goes extinct*)


extinctionAtReplication1 = Table[Select[extinctionPointsCEAR[[i]], #[[3]] == 1&][[All, 2]], {i, Length[extinctionPointsCEAR]}];


(* ::Text:: *)
(*Remove extinction points from FKP ESS points*)


fkpESSpointsCEAR = MapThread[Complement, {fkpESSpointsCEAR, extinctionAtReplication1}];


(* ::Text:: *)
(*Check that the ESS exists -- i.e. whether there are any retention probabilities where replication = 1 is an ESS*)


fkpESSexists = Map[Length[#] > 0&, fkpESSpointsCEAR];


(* ::Text:: *)
(*Get the minimum & maximum retention probabilities where replication = 1 is an ESS*)


fkpESSrangeCEAR = Map[MinMax, fkpESSpointsCEAR];


(* ::Text:: *)
(*Also get the maximum retention probability where consumers with replication = 1 can survive*)


maxSurvivableRetention = Table[Min[Select[extinctionPointsCEAR[[i]], #[[3]] == 1&][[All, 2]]] - 1/10, {i, Length[extinctionPointsCEAR]}];


(* ::Text:: *)
(*Draw rectangles marking the region where relication = 1 is an ESS*)


fkpRegionsCEAR = MapThread[If[#2, 
	Graphics[{FaceForm[], EdgeForm[{Thickness[0.03], Black}], 
		Rectangle[{Max[#1[[1]] - 1/20, 0], 1.015}, {Min[#1[[2]] + 1/20, 1], 1.085}, RoundingRadius -> 0.035]}],
	Graphics[] (* Draw nothing if fkp is not an ESS*)
	]&, 
	{fkpESSrangeCEAR, fkpESSexists}];


(* ::Subsection::Closed:: *)
(*Combine plots, extinction markers, and ESS markers*)


annotatedPlotsCE = MapThread[Show[{#1, #2, #3, #4}, 
	ImagePadding -> {{Automatic, 5}, {Automatic, 8}}, Background -> None]&, 
	{preyDensityPlotsCE, vectorPlotsCE, fkpRegionsCE, extinctionPlotCE}];
	
annotatedPlotsCE = MapThread[
	Labeled[#1, 
		Pane[Style[#2, FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], ImageSize -> 160, Alignment -> Left],
		{Top, Left}]&,
	{annotatedPlotsCE, {"(A)", "(B)", "(C)"}}]


annotatedPlotsCEAR = MapThread[Show[{#1, #2, #3, okpPoint, #4}, 
	ImagePadding -> {{Automatic, 5}, {Automatic, 8}}, Background -> None]&, 
	{preyDensityPlotsCEAR, vectorPlotsCEAR, fkpRegionsCEAR, extinctionPlotCEAR}];
	
annotatedPlotsCEAR = MapThread[
	Labeled[#1, 
		Pane[Style[#2, FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], ImageSize -> 160, Alignment -> Left],
		{Top, Left}]&,
	{annotatedPlotsCEAR, {"(D)", "(E)", "(F)"}}]


(* ::Subsection::Closed:: *)
(*Legend*)


myArrowhead = Graphics[Triangle[{{0, 0}, {-1, -1/4}, {-1, 1/4}}]];


legend = Framed[Column[{
	Pane[Style["Legend", FontSize -> 12, FontFamily -> "Arial", FontWeight -> Bold, FontColor -> Black], 
		ImageSize -> 72 * 1.4, Alignment -> -0.9],
		
	BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Column", 
	"TickSide" -> Left,
	LegendLabel -> Placed[Pane["Prey density (cells/\!\(\*SuperscriptBox[\(m\), \(2\)]\))", ImageSize -> 72 * 1.2], Right],
	LegendMarkerSize -> {24, 144}, LabelStyle -> {Black, "Arial", 12}],
	
	PointLegend[{Darker[Gray], Black, Black, Black}, 
	{Pane["Direction of selection", ImageSize -> 72 * 1.4],
		Pane["Consumer extinction", ImageSize -> 72 * 1.4],
		Pane["ESS: Obligate kleptoplasty", ImageSize -> 72 * 1.4],
		Pane["ESS: Facultative kleptoplasty", ImageSize -> 72 * 1.4]}, 
	LegendMarkers -> {Graphics[{Arrowheads[{{0.35, 1, myArrowhead}}], Arrow[{{0, 0}, {16, 0}}]}, AspectRatio -> 1/8], 
		Graphics[Rectangle[{0, 0}]],
		Graphics[{FaceForm[], Disk[]}], 
		Graphics[{FaceForm[], Rectangle[{0, 0}, {2, 1}, RoundingRadius -> 0.5]}]},
	LegendMarkerSize -> {{32, 16}, 16, 20, {32, 16}},
	LabelStyle -> {Black, "Arial", 12},
	BaseStyle -> EdgeForm[{Thick, Opacity[1]}]]}, (* This makes the edges of the rectangle and circle look right *)
			
	Alignment -> {Right, Center}, Spacings -> 0],
RoundingRadius -> 10]


(* ::Subsection::Closed:: *)
(*Axis labels and grid label*)


lightLabel = Graphics[{Thickness[0.0035], Arrowheads[0.025], Line[{{0, 0}, {200, 0}}], Arrow[{{344, 0}, {544, 0}}], 
Inset[Pane[Style["Increasing light", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], Alignment -> Center]]}]


ylabel = Rotate[Pane[Style["Replication (\[Eta])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12]], Pi/2]
xlabel = Pane[Style["Retention (\[Phi])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], 
	ImageSize -> 2.5 * 72, Alignment -> {0.2, Center}]


(* ::Subsection::Closed:: *)
(*Compile plots into a grid*)


compiledPlots = Grid[{{Grid[ (* We do a grid of a grid to get around an export glitch *)
	{Join[{ylabel}, annotatedPlotsCE], 
	Join[{ylabel}, annotatedPlotsCEAR],
	{, xlabel, xlabel, xlabel},
	{, Pane[Show[lightLabel, ImageSize -> 500], Alignment -> 0.38, ImageSize -> 550], SpanFromLeft, SpanFromLeft}},
	Alignment -> {Center, Center},
	Spacings -> {{0(* leftmost space, before grid starts *), 0 (* between ylabels and plots *), 
		Scaled[1/150] (* between plots *), Scaled[1/150] (* between plots *)}, 
		{0 (* uppermost space, before grid starts *), Scaled[5/144] (* between plots *), Scaled[1/150] (* between xlabels and plots *),
		Scaled[1/150] (* between arrow and xlabels *)}}],
	legend}}, Spacings -> 1]


Export["marked_vectorplots.pdf", compiledPlots]


(* ::Section::Closed:: *)
(*Supplementary plot showing magnitude of selection gradient*)


(* ::Subsubsection::Closed:: *)
(*Selection gradient colors*)


(* ::Text:: *)
(*Find the minimum and maximum vector magnitudes*)


vectorMagnitudesCE = Map[Norm, derivsForPlotsCE[[All, All, 2]], {2}];
vectorMagnitudesCEAR = Map[Norm,
								(* Remove derivatives that are "NA" before finding magnitudes *)
								DeleteCases[derivsForPlotsCEAR[[All, All, 2]], {"NA", "NA"}, 2], 
							{2}];


vectorMagnitudeRange = MinMax[{vectorMagnitudesCE, vectorMagnitudesCEAR}];


(* ::Text:: *)
(*Although the range is broad, most vector magnitudes are on the small end of the range*)


Quantile[Flatten[{vectorMagnitudesCE, vectorMagnitudesCEAR}], {0.9, 0.925, 0.95, 0.99}]


(* ::Text:: *)
(*Create a color scheme based on vector magnitudes*)


vectorMagnitudeColors[m_] := ColorData[{"ValentineTones", "Reverse"}][If[m <= 60, m/65, 1]];


(* ::Subsubsection::Closed:: *)
(*Selection vector plots*)


vectorMagnitudePlotsCE = Table[ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCE[[i]]], 
	VectorColorFunction -> vectorMagnitudeColors, VectorColorFunctionScaling -> False,
	VectorPoints -> noExtinctionPointsCE[[i, All, 2;;3]],
	VectorScale -> {Automatic, Automatic, None},
	VectorRange -> {0, 60}, ClippingStyle -> vectorMagnitudeColors[65],
	PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.05,
	LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.25],
{i, Length[derivsForPlotsCE]}];


vectorMagnitudePlotsCEAR = Table[ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivsForPlotsCEAR[[i]]],
	VectorColorFunction -> vectorMagnitudeColors, VectorColorFunctionScaling -> False,
	VectorPoints -> noExtinctionPointsCEAR[[i, All, 2;;3]],
	VectorScale -> {Automatic, Automatic, None},
	VectorRange -> {0, 60}, ClippingStyle -> vectorMagnitudeColors[65],
	PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.05,
	LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.25],
{i, Length[derivsForPlotsCEAR]}];


(* ::Subsubsection::Closed:: *)
(*Combine plots with markers showing extinction points*)


extinctionPlotCE = Table[
Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, 
extinctionPointsCE[[i, All, 2;;3]]]],

{i, Length[derivsForPlotsCE]}];


extinctionPlotCEAR = Table[
Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, 
extinctionPointsCEAR[[i, All, 2;;3]]]],

{i, Length[derivsForPlotsCEAR]}];


vectorMagnitudeMarkedPlotsCE = MapThread[Show, {vectorMagnitudePlotsCE, extinctionPlotCE}];
vectorMagnitudeMarkedPlotsCEAR = MapThread[Show, {vectorMagnitudePlotsCEAR, extinctionPlotCEAR}];


(* ::Subsubsection::Closed:: *)
(*Label plots*)


vectorMagnitudeLabeledPlotsCE = MapThread[
	Labeled[#1, Pane[Style[#2, 14, Black, "Arial"], ImageSize -> 72 * 2.25], Top]&,
	{vectorMagnitudeMarkedPlotsCE, {"(A)", "(B)", "(C)"}}];


vectorMagnitudeLabeledPlotsCEAR = MapThread[
	Labeled[#1, Pane[Style[#2, 14, Black, "Arial"], ImageSize -> 72 * 2.25], Top]&,
	{vectorMagnitudeMarkedPlotsCEAR, {"(D)", "(E)", "(F)"}}];


(* ::Subsubsection::Closed:: *)
(*Legend*)


(* ::Text:: *)
(*Make a legend to go with the color scheme*)


vectorMagnitudeLegend = BarLegend[{vectorMagnitudeColors, {0, 70}},
  "Ticks" -> Table[i, {i, 0, 70, 10}],
 "TickLabels" -> Append[Table[If[Mod[i, 10] == 0, ToString[i], ""], {i, 0, 60, 10}], ">60"],
 LegendLabel -> "Gradient\nmagnitude",
 LegendMarkerSize -> {20, 72*3},
 LabelStyle -> {12, Black, "Arial"}];


(* ::Text:: *)
(*An arrowhead to represent the direction of selection*)


myArrowhead = Graphics[Triangle[{{0, 0}, {-1, -1/4}, {-1, 1/4}}]];


(* ::Text:: *)
(*Full legend*)


vectorMagnitudeFullLegend = Framed[Column[{
	Pane[Style["Legend", FontSize -> 12, FontFamily -> "Arial", FontWeight -> Bold, FontColor -> Black], 
		ImageSize -> 72 * 1, Alignment -> Center],
		
	vectorMagnitudeLegend,
	
	PointLegend[{Darker[Gray], Black, Black, Black}, 
	{Pane["Direction of selection", ImageSize -> 72 * 1],
		Pane["Consumer extinction", ImageSize -> 72 * 1]}, 
	LegendMarkers -> {Graphics[{Arrowheads[{{0.35, 1, myArrowhead}}], Arrow[{{0, 0}, {16, 0}}]}, AspectRatio -> 1/8], 
		Graphics[Rectangle[{0, 0}]]},
	LegendMarkerSize -> {{32, 16}, 16},
	LabelStyle -> {Black, "Arial", 12},
	BaseStyle -> EdgeForm[{Thick, Opacity[1]}]]}, (* This makes the edges of the rectangle and circle look right *)
			
	Alignment -> {Center, Center}, Spacings -> 0],
RoundingRadius -> 10]


(* ::Subsubsection::Closed:: *)
(*Axis labels and grid label*)


lightLabel = Graphics[{Thickness[0.0035], Arrowheads[0.025], Line[{{0, 0}, {200, 0}}], Arrow[{{344, 0}, {544, 0}}], 
Inset[Pane[Style["Increasing light", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], Alignment -> Center]]}]


ylabel = Rotate[Pane[Style["Replication (\[Eta])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12]], Pi/2];
xlabel = Pane[Style["Retention (\[Phi])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], 
	ImageSize -> 2.5 * 72, Alignment -> {0.2, Center}];


(* ::Subsubsection::Closed:: *)
(*Compile plots into a grid*)


compiledVectorMagnitudePlots = Grid[{{Grid[ (* We do a grid of a grid to get around an export glitch *)
	{Join[{ylabel}, vectorMagnitudeLabeledPlotsCE], 
	Join[{ylabel}, vectorMagnitudeLabeledPlotsCEAR],
	{, xlabel, xlabel, xlabel},
	{, Pane[Show[lightLabel, ImageSize -> 500], Alignment -> 0.38, ImageSize -> 550], SpanFromLeft, SpanFromLeft}},
	Alignment -> {Center, Center},
	Spacings -> {{0(* leftmost space, before grid starts *), 0 (* between ylabels and plots *), 
		Scaled[1/150] (* between plots *), Scaled[1/150] (* between plots *)}, 
		{0 (* uppermost space, before grid starts *), Scaled[1/50] (* between plots *), Scaled[1/150] (* between xlabels and plots *),
		Scaled[1/150] (* between arrow and xlabels *)}}],
	vectorMagnitudeFullLegend}}, Spacings -> 1]


Export["vector_magnitude_plots.pdf", compiledVectorMagnitudePlots]


(* ::Section::Closed:: *)
(*Extinction at high retention & replication is due to presence of prey*)


(* ::Text:: *)
(*If we approach extinction points from neighboring points, prey are present and the mixotroph population goes extinct.*)
(*If we approach via a "jump" (large mutation), prey might not be present, and the mixotroph population can survive.*)


neighborPt = ecoDynamicsCEAR[[3, 120, 5;;7]]
jumpPt = ecoDynamicsCEAR[[3, 66, 5;;7]]


neighborDynamics = residentEcoDynamics[1, 1, 50, neighborPt[[1]], neighborPt[[2]], neighborPt[[3]], attackTradeOff -> True, conversionTradeOff -> True]
jumpDynamics = residentEcoDynamics[1, 1, 50, jumpPt[[1]], jumpPt[[2]], jumpPt[[3]], attackTradeOff -> True, conversionTradeOff -> True]


(* ::Text:: *)
(*Extinction when prey are present*)


Plot[If[t <= 32.9, ((h[t] + p[t])/. neighborDynamics[[1]]), 0], {t, 0, 4000}]


(* ::Text:: *)
(*Survival when prey are not present*)


Plot[((h[t] + p[t])/. jumpDynamics[[1]]), {t, 0, 4000}, PlotRange -> {{0, 4000}, {0, 25}}]


(* ::Text:: *)
(*Introduce prey into prey-free equilibrium. We can see that the prey-free equilibrium is unstable, and the mixotroph population goes extinct when prey are introduced*)


preyFreeEquil = {w[4000], h[4000], p[4000]} /. jumpDynamics[[1]]
preyIntroductionPt = preyFreeEquil + {0.000001, 0, 0}


preyIntroductionDynamics = residentEcoDynamics[1, 1, 50, preyIntroductionPt[[1]], preyIntroductionPt[[2]], preyIntroductionPt[[3]], 
	attackTradeOff -> True, conversionTradeOff -> True]
