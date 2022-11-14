#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section::Closed:: *)
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
	0 < retention < 1, "Mainly kleptoplastidic mixotroph"],
	replication == 1, (* If replication evolves to 1, the mixotroph is either a phototroph (does not retain plastid), a kleptoplastidic phototroph (always retains plastid) or an occasional kleptoplastidic phototroph (sometimes retains plastid *)
	Which[retention == 0, "Strict phototroph",
	retention > 0, "Facultatively kleptoplastidic phototroph"],
	0 < replication < 1, "No label"] (* We don't have a name for any state where replication is intermediate *) 


(* ::Text:: *)
(*A list of the outcomes we observe (we'll still be able to tell if things didn't get one of these labels, though)*)


(* outcomeLabels = {"Heterotroph", "Strict Kleptoplast", "Phototroph", "Kleptoplastidic Phototroph", "Occasional Kleptoplastidic Phototroph", "Extinction"}; *)
outcomeLabels = {(*"Strict kleptoplast", "Obligately kleptoplastidic mixotroph",*) "Mainly kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph", "Extinction"};


(* ::Text:: *)
(*Colors for plotting*)


(* myColors = "CandyColors"; *)
myColors = Append[(ColorData["CandyColors"] /@ Subdivide[3])[[{(*1,*) 2, 4}]], White];


(* ::Section:: *)
(*Some of the simulations look like they reached an ESS (or something?) in the middle of phi-eta space*)


(* ::Subsection:: *)
(*Import saved simulations*)


evolutionSimsCEAR = Import[NotebookDirectory[] <> "evolutionSimsCEAR.wdx"];


(* ::Text:: *)
(*Sort simulations by light level*)


gEvoSimsCEAR = GatherBy[SortBy[evolutionSimsCEAR, #[[1]]&], #[[1]]&];
simEndpointsCEAR = gEvoSimsCEAR[[All, All, 2, 1;;2, -1]];


(* ::Subsection::Closed:: *)
(*Classify simulations by endpoint*)


(* ::Text:: *)
(*Classify simulations*)


outcomesCEAR = Map[classifyOutcome, simEndpointsCEAR, {2}];


(* ::Text:: *)
(*Count the classifications*)


outcomeCountsCEAR = Outer[Count, outcomesCEAR, Append[outcomeLabels, "No label"], 1];


(* ::Text:: *)
(*Turn the counts into fraction of simulation*)


outcomeCountsCEAR = outcomeCountsCEAR / Map[Length, outcomesCEAR];


(* ::Subsection:: *)
(*Investigate the simulations that don't fall into our previous known ESSs*)


(* ::Text:: *)
(*I called these "stopped in the middle sims" because they stop evolving in the middle of parameter space 0 < phi, eta < 1*)


stoppedInTheMiddleSims = Pick[gEvoSimsCEAR, outcomesCEAR, "No label"];


stoppedInTheMiddleLight = stoppedInTheMiddleSims[[All, All, 1]];


stoppedInTheMiddlePhiEta = stoppedInTheMiddleSims[[All, All, 2, 1;;2]];


(* ::Subsubsection:: *)
(*Plots of the evolutionary trajectories*)


Map[ListLinePlot[#, 
	Frame -> {True, True, False, False}, FrameLabel -> {"Mutations", "Trait"},
	PlotRange -> {{0, 500}, {0, 1}},
	PlotStyle -> {Blue, Orange},
	LabelStyle -> {FontSize -> 12, FontColor -> Black, FontFamily -> "Arial"}]&, 
	stoppedInTheMiddlePhiEta, {2}] // Flatten
LineLegend[{Blue, Orange}, {"Retention", "Replication"}, LabelStyle -> {FontSize -> 12, FontColor -> Black, FontFamily -> "Arial"}]


(* ::Subsubsection:: *)
(*Example ecological dynamics*)


(* ::Text:: *)
(*Light*)


stoppedInTheMiddleLight[[5, 1]]


(* ::Text:: *)
(*Final retention and replication*)


stoppedInTheMiddlePhiEta[[5, 1, All, -1]]


(* ::Text:: *)
(*Evo. trajectory*)


stoppedInTheMiddleExampleEvoPlot = ListLinePlot[stoppedInTheMiddlePhiEta[[5, 1]],
	PlotRange -> {{0, 500}, {0, 1}},
	Frame -> {True, True, False, False},
	FrameLabel -> {"Mutations", "Trait value"},
	PlotLabel -> "Example of possible new ESS",
	PlotLegends -> {"Retention", "Replication"},
	LabelStyle -> {FontSize -> 12, FontColor -> Black, FontFamily -> "Arial"},
	ImageSize -> 300]


(* ::Text:: *)
(*Eco dynamics at end of trajectory*)


stoppedInTheMiddleExampleEco = residentEcoDynamics[21/25, 19/100, 40, 4, 1, 1, conversionTradeOff -> True, attackTradeOff -> True];


stoppedInTheMiddleExampleEcoPlot  = Plot[Evaluate[{w[t], h[t] + p[t]} /. stoppedInTheMiddleExampleEco[[1]]], {t, 3000, 4000},
Frame -> {True, True, False, False}, FrameLabel -> {"Time (days)", "Population (cells m^-2)"},
PlotLegends -> {"Prey", "Consumer"},
PlotLabel -> "Ecological dynamics at end of simulation",
LabelStyle -> {FontSize -> 12, FontColor -> "Black", FontFamily -> "Arial"},
ImageSize -> 300]


(* ::Text:: *)
(*Combined plot of eco/evo dynamics*)


Column[{stoppedInTheMiddleExampleEvoPlot, stoppedInTheMiddleExampleEcoPlot}]
Export["example_possible_new_ess.pdf", %]



