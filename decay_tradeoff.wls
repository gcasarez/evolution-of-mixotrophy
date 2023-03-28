#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
(*Parameters*)


(* ::Text:: *)
(*Note on the parameter c_HP in the paper: In the code we have replaced c_HP with its reciprocal, gamma_C. This doesn't come up anywhere for the user, but just be aware that setting gamma_C = 0.7 in the functions below is equivalent to setting c_HP to 10/7. *)


pW = 3;
gammaP = 2/3;(*2/3;*)(*2.9/3;*)
(*sigmaP = 0;*)
kW = 0.1;
kH = kH2 = 0.05;
kP = kP2 = 0.15;
hW = 50;(*50;*)
hP = hP2 = 50;(*10;*)(*50;*)
lW = 0.5;
lH = lH2 = 0.25;
lP = lP2 = 0.25;
dMax = 0.3;(*0.3;*)
(*sigmaD = 1; (* Decay trade-off *)*)
a = a2 = 0.15;
gammaA = 0; (* Phototroph does not eat *)
(*sigmaA = 0;*) (* 0 for no attack rate trade-off, 1 for attack rate trade-off *)
c = 0.1;
(*gammaC = 1;*) (* 1 for no conversion efficiency trade-off, 0.7 for conversion efficiency trade-off *)
sigmaC = 0; (* Added sigmaC because it looks like it's not defined here *)
(*iIn = 25; *)(* Incident light *)


(* ::Section::Closed:: *)
(*Equations*)


Get[NotebookDirectory[] <> "model_and_functions.wls"]


lightLevels = {35};


(* ::Text:: *)
(*Function to find the selection gradient, given the attack trade-off*)


selectionGradient[attackTO_] :=
Module[{ecoDynamics, extinctionPoints, noExtinctionPoints, derivsForPlots},

	(* Ecological dynamics *)
	ecoDynamics = Table[Join[{light, retention, replication}, averageEcoState[retention, replication, light, 4, 1, 1, PrecisionGoal -> 6, AccuracyGoal -> 6, 
			conversionTradeOff -> True, attackTradeOff -> attackTO]],
		{light, lightLevels}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}];
		
	ecoDynamics = Flatten[ecoDynamics, {{1}, {2,3}, {4}}];
	
	(* Where does the mixotroph go extinct? *)
	extinctionPoints = Map[Select[#[[4]] == True&], ecoDynamics];
	(* Also find the points where the mixotroph population does not go extinct *)
	noExtinctionPoints = Map[Select[#[[4]] == False&], ecoDynamics];
	
	(* Selection gradient *)
	derivsForPlots = 
		ParallelMap[{{#[[2]], #[[3]]}, invasionDerivatives[#[[2]], #[[3]], #[[1]], #[[5]], #[[6]], #[[7]], PrecisionGoal -> 6, AccuracyGoal -> 6, MaxSteps -> 40000, 
			conversionTradeOff -> True, attackTradeOff -> attackTO],
		{#[[5]], #[[6]], #[[7]]}}&, noExtinctionPoints, {2}];
		
	(* Return selection gradient & extinction points *)
	Return[{derivsForPlots, extinctionPoints}]
]


(* ::Section::Closed:: *)
(*Photosynthesis and decay not affected by retention or replication*)


(* ::Text:: *)
(*This just repeats what the "main" model is, for comparison*)


sigmaP = 0; (* no photosynthesis trade-off *)
sigmaD = 0; (* no decay trade-off *)


(* ::Text:: *)
(*Attack trade-off*)


gradientNormalCEAR = selectionGradient[True];


(* ::Text:: *)
(*No attack trade-off*)


gradientNormalCE = selectionGradient[False];


{ListVectorPlot[gradientNormalCEAR[[1, 1, All, 1;;2]]], ListVectorPlot[gradientNormalCE[[1, 1, All, 1;;2]]]}


(* ::Section::Closed:: *)
(*Replication improves photosynthesis or reduces plastid decay*)


(* ::Text:: *)
(*This trade-off is built into the equations in "model_and_functions.wls." We just have to set sigmaP and sigmaD to 1 to "turn on" the relationship.*)


(* ::Subsection::Closed:: *)
(*Replication improves photosynthesis*)


sigmaP = 1; (* yes photosynthesis trade-off *)
sigmaD = 0; (* no decay trade-off *)


(* ::Text:: *)
(*Attack trade-off*)


gradientReplicationPhotoCEAR = selectionGradient[True];


(* ::Text:: *)
(*No attack trade-off*)


gradientReplicationPhotoCE = selectionGradient[False];


{ListVectorPlot[gradientReplicationPhotoCEAR[[1, 1, All, 1;;2]]], ListVectorPlot[gradientReplicationPhotoCE[[1, 1, All, 1;;2]]]}


(* ::Subsection::Closed:: *)
(*Replication reduces plastid decay*)


sigmaP = 0; (* no photosynthesis trade-off *)
sigmaD = 1; (* yes decay trade-off *)


(* ::Text:: *)
(*Attack trade-off*)


gradientReplicationDecayCEAR = selectionGradient[True];


(* ::Text:: *)
(*No attack trade-off*)


gradientReplicationDecayCE = selectionGradient[False];


{ListVectorPlot[gradientReplicationDecayCEAR[[1, 1, All, 1;;2]]], ListVectorPlot[gradientReplicationDecayCE[[1, 1, All, 1;;2]]]}


(* ::Section::Closed:: *)
(*Retention improves photosynthesis or reduces plastid decay*)


(* ::Text:: *)
(*Now we instead model rentention (phi) as improving consumer photosynthesis and/or reducing plastid decay.*)
(*Again we can use sigmaD = 1 or sigmaP = 1 to "turn on" the relationship between retention and decay or photosynthesis, respectively.*)


(* ::Subsection::Closed:: *)
(*Model equations*)


(* ::Text:: *)
(*Prey*)


wDot[w_, hW_, lW_, a_, pW_, kappa_, iIn_, gammaA_, eta_, phi_, eta2_, phi2_, sigmaA_, p_, h_, a2_, h2_, p2_]:= 
(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-(a(1-sigmaA (phi + eta)/2)) w h-(gammaA a(1-sigmaA (phi + eta)/2)) w p-(a2(1-sigmaA (phi2 + eta2)/2)) w h2 - (gammaA a2(1-sigmaA (phi2 + eta2)/2)) w p2 


(* ::Text:: *)
(*Resident heterotroph*)


hDot[eta_, c_, gammaC_, kappa_, gammaP_, pW_, sigmaP_, p_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, lH_, dMax_, sigmaD_]:=
(1-eta) (c/(gammaC c))((gammaP pW + sigmaP phi(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi (a(1-sigmaA (phi + eta)/2)) w h + (1-phi)c (a(1-sigmaA (phi + eta)/2)) w h+ (1-eta) c (gammaA a(1-sigmaA (phi + eta)/2)) w p -lH h + (dMax (1-sigmaD phi))p


(* ::Text:: *)
(*Resident phototroph*)


pDot[eta_, gammaP_, sigmaP_, pW_, p_, kappa_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, gammaC_, sigmaC_, c_, dMax_, sigmaD_, lP_] :=
eta ((gammaP pW + sigmaP phi(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi (a(1-sigmaA (phi + eta)/2)) w h +eta (gammaA a(1-sigmaA (phi + eta)/2)) (gammaC c) p w - (dMax (1-sigmaD phi)) p - lP p


(* ::Text:: *)
(*Mutant heterotroph*)


h2Dot[eta2_, c_, gammaC_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, lH2_, dMax_, sigmaD_] :=
(1-eta2) (c/(gammaC c))((gammaP pW + sigmaP phi2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 (a2(1-sigmaA (phi2 + eta2)/2)) w h2 + (1-phi2) c (a2(1-sigmaA (phi2 + eta2)/2)) w h2 +(1-eta2) c (gammaA a2(1-sigmaA (phi2 + eta2)/2)) w p2 - lH2 h2+ (dMax(1-sigmaD phi2)) p2


(* ::Text:: *)
(*Mutant phototroph*)


p2Dot[eta2_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, c_, gammaC_, sigmaC_, dMax_, sigmaD_, lP2_] :=
eta2 ((gammaP pW + sigmaP phi2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 (a2(1-sigmaA (phi2 + eta2)/2)) w h2+eta2 (gammaA a2(1-sigmaA (phi2 + eta2)/2)) (gammaC c) p2 w-(dMax(1-sigmaD phi2)) p2-lP2 p2


(* ::Subsection::Closed:: *)
(*Retention improves photosynthesis*)


sigmaP = 1; (* yes photosynthesis trade-off *)
sigmaD = 0; (* no decay trade-off *)


(* ::Text:: *)
(*Attack trade-off*)


gradientRetentionPhotoCEAR = selectionGradient[True];


(* ::Text:: *)
(*No attack trade-off*)


gradientRetentionPhotoCE = selectionGradient[False];


{ListVectorPlot[gradientRetentionPhotoCEAR[[1, 1, All, 1;;2]]], ListVectorPlot[gradientRetentionPhotoCE[[1, 1, All, 1;;2]]]}


(* ::Subsection::Closed:: *)
(*Retention reduces plastid decay*)


sigmaP = 0; (* no photosynthesis trade-off *)
sigmaD = 1; (* yes decay trade-off *)


(* ::Text:: *)
(*Attack trade-off*)


gradientRetentionDecayCEAR = selectionGradient[True];


(* ::Text:: *)
(*No attack trade-off*)


gradientRetentionDecayCE = selectionGradient[False];


{ListVectorPlot[gradientRetentionDecayCEAR[[1, 1, All, 1;;2]]], ListVectorPlot[gradientRetentionDecayCE[[1, 1, All, 1;;2]]]}


(* ::Section::Closed:: *)
(*Retention and replication improve photosynthesis or reduces plastid decay*)


(* ::Text:: *)
(*Here the sum of retention and replication improves consumer photosynthesis and reduces plastid decay.*)


(* ::Subsection::Closed:: *)
(*Model functions*)


(* ::Text:: *)
(*Prey*)


wDot[w_, hW_, lW_, a_, pW_, kappa_, iIn_, gammaA_, eta_, phi_, eta2_, phi2_, sigmaA_, p_, h_, a2_, h2_, p2_]:= 
(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-(a(1-sigmaA (phi + eta)/2)) w h-(gammaA a(1-sigmaA (phi + eta)/2)) w p-(a2(1-sigmaA (phi2 + eta2)/2)) w h2 - (gammaA a2(1-sigmaA (phi2 + eta2)/2)) w p2 


(* ::Text:: *)
(*Resident heterotroph*)


hDot[eta_, c_, gammaC_, kappa_, gammaP_, pW_, sigmaP_, p_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, lH_, dMax_, sigmaD_]:=
(1-eta) (c/(gammaC c))((gammaP pW + sigmaP ((phi + eta)/2) (pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi (a(1-sigmaA (phi + eta)/2)) w h + (1-phi)c (a(1-sigmaA (phi + eta)/2)) w h+ (1-eta) c (gammaA a(1-sigmaA (phi + eta)/2)) w p -lH h + (dMax (1-sigmaD ((phi + eta)/2)))p


(* ::Text:: *)
(*Resident phototroph*)


pDot[eta_, gammaP_, sigmaP_, pW_, p_, kappa_, hP_, iIn_, phi_, a_, w_, h_, gammaA_, sigmaA_, gammaC_, sigmaC_, c_, dMax_, sigmaD_, lP_] :=
eta ((gammaP pW + sigmaP ((phi + eta)/2) (pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi (a(1-sigmaA (phi + eta)/2)) w h +eta (gammaA a(1-sigmaA (phi + eta)/2)) (gammaC c) p w - (dMax (1-sigmaD ((phi + eta)/2))) p - lP p


(* ::Text:: *)
(*Mutant heterotroph*)


h2Dot[eta2_, c_, gammaC_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, lH2_, dMax_, sigmaD_] :=
(1-eta2) (c/(gammaC c))((gammaP pW + sigmaP ((phi2 + eta2)/2) (pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 (a2(1-sigmaA (phi2 + eta2)/2)) w h2 + (1-phi2) c (a2(1-sigmaA (phi2 + eta2)/2)) w h2 +(1-eta2) c (gammaA a2(1-sigmaA (phi2 + eta2)/2)) w p2 - lH2 h2+ (dMax(1-sigmaD ((phi2 + eta2)/2))) p2


(* ::Text:: *)
(*Mutant phototroph*)


p2Dot[eta2_, gammaP_, sigmaP_, pW_, p2_, kappa_, hP2_, iIn_, phi2_, a2_, w_, h2_, gammaA_, sigmaA_, c_, gammaC_, sigmaC_, dMax_, sigmaD_, lP2_] :=
eta2 ((gammaP pW + sigmaP ((phi2 + eta2)/2) (pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 (a2(1-sigmaA (phi2 + eta2)/2)) w h2+eta2 (gammaA a2(1-sigmaA (phi2 + eta2)/2)) (gammaC c) p2 w-(dMax(1-sigmaD ((phi2 + eta2)/2))) p2-lP2 p2


(* ::Subsection::Closed:: *)
(*Retention and replication improve photosynthesis*)


sigmaP = 1; (* yes photosynthesis trade-off *)
sigmaD = 0; (* no decay trade-off *)


(* ::Text:: *)
(*Attack trade-off*)


gradientBothPhotoCEAR = selectionGradient[True];


(* ::Text:: *)
(*No attack trade-off*)


gradientBothPhotoCE = selectionGradient[False];


{ListVectorPlot[gradientBothPhotoCEAR[[1, 1, All, 1;;2]]], ListVectorPlot[gradientBothPhotoCE[[1, 1, All, 1;;2]]]}


(* ::Subsection::Closed:: *)
(*Retention and replication reduce plastid decay*)


sigmaP = 0; (* no photosynthesis trade-off *)
sigmaD = 1; (* yes decay trade-off *)


(* ::Text:: *)
(*Attack trade-off*)


gradientBothDecayCEAR = selectionGradient[True];


(* ::Text:: *)
(*No attack trade-off*)


gradientBothDecayCE = selectionGradient[False];


{ListVectorPlot[gradientBothDecayCEAR[[1, 1, All, 1;;2]]], ListVectorPlot[gradientBothDecayCE[[1, 1, All, 1;;2]]]}


(* ::Section::Closed:: *)
(*Plots*)


(* ::Subsection::Closed:: *)
(*Colors and legend*)


(* ::Subsubsection::Closed:: *)
(*Set font*)


SetOptions[ListDensityPlot, LabelStyle -> {FontSize -> 12, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions -> {"OperatorSubstitution"->False}}];


SetOptions[{SwatchLegend, PointLegend, BarLegend}, LabelStyle -> {FontSize -> 14, FontFamily -> "Arial", FontColor -> Black, PrivateFontOptions->{"OperatorSubstitution"->False}}];


(* ::Subsubsection::Closed:: *)
(*Color scheme (prey density)*)


Map[Function[gradient, Flatten[gradient[[1, All, All, 3, 1]]]],
	{gradientNormalCEAR, gradientNormalCE, gradientReplicationPhotoCEAR, gradientReplicationPhotoCE,
	gradientReplicationDecayCEAR, gradientReplicationDecayCE, gradientRetentionPhotoCEAR,
	gradientRetentionPhotoCE, gradientRetentionDecayCEAR, gradientRetentionDecayCE,
	gradientBothPhotoCEAR, gradientBothPhotoCE, gradientBothDecayCEAR, gradientBothDecayCE}];
	
maxPrey = Max[%]


maxPrey = Max[maxPrey, 40] (* To match other plots *)


preyDensityColors[p_] := ColorData[{"BeachColors", "Reverse"}][Log[(p + 1)]/Log[maxPrey + 1]]


(* ::Subsubsection::Closed:: *)
(*Legend*)


legend = Framed[Column[{
	Pane[Style["Legend", FontSize -> 12, FontFamily -> "Arial", FontWeight -> Bold, FontColor -> Black], 
		ImageSize -> 72 * 1, Alignment -> Center],
		
	BarLegend[{preyDensityColors[#]&, {0, maxPrey}}, LegendLayout -> "Column", 
	LegendLabel -> Placed[Pane["Prey density (cells/\!\(\*SuperscriptBox[\(m\), \(2\)]\))", ImageSize -> 72 * 1], Top],
	LegendMarkerSize -> {24, 72 * 2}, LabelStyle -> {Black, "Arial", 12}],
	
	PointLegend[{Darker[Gray], Black}, 
	{Pane["Direction of selection", ImageSize -> 72 * 1],
	Pane["Consumer extinction", ImageSize -> 72 * 1]}, 
	LegendMarkers -> {Graphics[{Arrowheads[{{0.35, 1, myArrowhead}}], Arrow[{{0, 0}, {16, 0}}]}, AspectRatio -> 1/8],
		Graphics[Rectangle[{0, 0}]]},
	LegendMarkerSize -> {16},
	LabelStyle -> {Black, "Arial", 12}]},
	
	Alignment -> {Center, Center}, Spacings -> 0],
RoundingRadius -> 10]


(* ::Subsection::Closed:: *)
(*Make plots*)


plots = Map[Function[gradient, 
	Module[{preyPlot, vectorPlot, extinctionPlot, noExtinctionPoints, derivatives = gradient[[1, 1]], extinction = gradient[[2, 1]]},
	
	(* Plot prey population size *)
	preyPlot = Show[Graphics[Inset[
		Rasterize[ListDensityPlot[Map[{#[[1, 1]], #[[1, 2]], #[[3, 1]]}&, derivatives],
			PlotRange -> {{0, 1}, {0, 1}, Full}, InterpolationOrder -> 1,
			ColorFunction -> preyDensityColors, ColorFunctionScaling -> False,
			PlotRangePadding -> 0, ImagePadding -> 0, ImageMargins -> 0, FrameLabel -> None, Frame -> None]],
		{0, 0}, {0, 0}, Scaled[1/1.1]], 
		PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.05, ImagePadding -> Full], 
		Frame -> True, FrameStyle -> Black, LabelStyle -> {Black, "Arial", 12}, ImageSize -> 72*2.25];
		
	(* Points where consumer extinction does not occur *)
	noExtinctionPoints = Complement[Flatten[Table[{retention, replication}, {retention, 0, 1, 1/10}, {replication, 0, 1, 1/10}], 1],
		extinction[[All, 2;;3]]];
	
	(* Plot selection gradient *)
	vectorPlot = ListVectorPlot[Map[{#[[1]], #[[2]]}&, derivatives], 
	VectorStyle -> Darker[Gray], VectorColorFunction -> None, VectorPoints -> noExtinctionPoints,
	VectorScale -> {Automatic, Automatic, None}];
	
	(* Markers for points where consumer extinction occurs *)
	extinctionPlot = Graphics[Map[Rectangle[{Max[#[[1]] - 1/10, 0], Max[#[[2]] - 1/10, 0]}, {Min[#[[1]] + 1/10, 1], Min[#[[2]] + 1/10, 1]}]&, 
		extinction[[All, 2;;3]]]];
		
	(* Combine everything into one plot *)
	Show[{preyPlot, vectorPlot, extinctionPlot}, ImagePadding -> {{Automatic, 5}, {Automatic, 8}}]

]],

{
 {gradientNormalCEAR, gradientReplicationDecayCEAR, gradientReplicationPhotoCEAR},
 {gradientNormalCEAR, gradientBothDecayCEAR, gradientBothPhotoCEAR}}, {2}]


 (* Full list of all plots that can be made. Just substitute for the list of gradients above to plot everything *)
 (*{{gradientNormalCE,   gradientReplicationDecayCE,   gradientReplicationPhotoCE},
 {gradientNormalCEAR, gradientReplicationDecayCEAR, gradientReplicationPhotoCEAR},
 {gradientNormalCE,   gradientRetentionDecayCE,   gradientReplicationPhotoCE},
 {gradientNormalCEAR, gradientRetentionDecayCEAR, gradientRetentionPhotoCEAR},
 {gradientNormalCE,   gradientBothDecayCE,   gradientBothPhotoCE},
 {gradientNormalCEAR, gradientBothDecayCEAR, gradientBothPhotoCEAR}}*)


(* ::Subsubsection::Closed:: *)
(*Label plots*)


labeledPlots = MapThread[
	Labeled[#1, 
		Pane[Style[#2, FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], ImageSize -> 160, Alignment -> Left],
		{Top, Left}]&,
	{plots, {{"(A)", "(B)", "(C)"}, {"(D)", "(E)", "(F)"}}}, 2];


ylabel = Rotate[Pane[Style["Replication (\[Eta])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12]], Pi/2]
xlabel = Pane[Style["Retention (\[Phi])", FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], 
	ImageSize -> 2.5 * 72, Alignment -> {0.2, Center}]


colLabel = Map[Rotate[Pane[Style[#, FontColor -> Black, FontFamily -> "Arial", FontSize -> 12]], Pi/2]&,
	{"Decreases plastid decay", "Increases photosynthesis"}]

rowLabel = Map[Pane[Style[#, FontColor -> Black, FontFamily -> "Arial", FontSize -> 12], 
	ImageSize -> 2 * 72, Alignment -> {0.2, Center}]&, 
	{"No effect", "Replication...", "Retention + replication..."}]


(* ::Subsubsection::Closed:: *)
(*Combine plots with legend and export*)


Grid[{Join[{"", ""}, rowLabel, {""}],
Join[{colLabel[[1]], ylabel}, labeledPlots[[1]], {legend}],
Join[{colLabel[[2]], ylabel}, labeledPlots[[2]], {SpanFromAbove}],
{"", "", xlabel, xlabel, xlabel, ""}},
Alignment -> {Center, Center}, Spacings -> {0.2, {0, 0, Scaled[1/50], Scaled[1/150]}}]


Export["decay_and_photosynthesis_tradeoffs.pdf", %]
