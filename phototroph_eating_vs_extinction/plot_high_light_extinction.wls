#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
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
myColors = Append[(ColorData["CandyColors"] /@ Subdivide[3])[[{(*1,*) 2, 4}]], Black];


(* ::Section:: *)
(*Import*)


phototrophDoesntEatCE = Import[NotebookDirectory[] <> "evolution_phototroph_doesnt_eat_CE_2.wdx"];


phototrophDoesntEatCEAR = Import[NotebookDirectory[] <> "evolution_phototroph_doesnt_eat_CEAR_2.wdx"];


phototrophEatsCE = Import[NotebookDirectory[] <> "evolution_phototroph_eats_CE_2.wdx"];


phototrophEatsCEAR = Import[NotebookDirectory[] <> "evolution_phototroph_eats_CEAR_2.wdx"];


Map[Dimensions, {phototrophDoesntEatCE, phototrophDoesntEatCEAR, phototrophEatsCE, phototrophEatsCEAR}]


(* ::Text:: *)
(*Light levels*)


allLightLevels = Map[DeleteDuplicates[#[[All, 1]]]&, 
	{phototrophDoesntEatCE, phototrophDoesntEatCEAR, phototrophEatsCE, phototrophEatsCEAR}];


(* ::Text:: *)
(*Check that all simulations have the same light levels*)


allLightLevels[[1]] == allLightLevels[[2]] == allLightLevels[[3]] == allLightLevels[[4]]

If[%, lightLevels = allLightLevels[[1]]]


(* ::Section:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


{gEvoSimsNoEatingCE, gEvoSimsNoEatingCEAR, gEvoSimsEatingCE, gEvoSimsEatingCEAR} = 
	Map[Function[{rawEvoSims}, GatherBy[SortBy[rawEvoSims, #[[1]]&], #[[1]]&]], 
		{phototrophDoesntEatCE, phototrophDoesntEatCEAR, phototrophEatsCE, phototrophEatsCEAR}];


{simEndpointsNoEatingCE, simEndpointsNoEatingCEAR, simEndpointsEatingCE, simEndpointsEatingCEAR} = 
	Map[Function[{gEvoSims}, gEvoSims[[All, All, 2, 1;;2, -1]]],
		{gEvoSimsNoEatingCE, gEvoSimsNoEatingCEAR, gEvoSimsEatingCE, gEvoSimsEatingCEAR}];


{outcomesNoEatingCE, outcomesNoEatingCEAR, outcomesEatingCE, outcomesEatingCEAR} = 
	Map[classifyOutcome, {simEndpointsNoEatingCE, simEndpointsNoEatingCEAR, simEndpointsEatingCE, simEndpointsEatingCEAR}, {3}];


{outcomeCountsNoEatingCE, outcomeCountsNoEatingCEAR, outcomeCountsEatingCE, outcomeCountsEatingCEAR} = 
	Map[Function[{outcomes}, Outer[Count, outcomes, outcomeLabels, 1] / Map[Length, outcomes]],
		{outcomesNoEatingCE, outcomesNoEatingCEAR, outcomesEatingCE, outcomesEatingCEAR}];


Map[Total, {outcomeCountsNoEatingCE, outcomeCountsNoEatingCEAR, outcomeCountsEatingCE, outcomeCountsEatingCEAR}, {2}]


(* ::Section:: *)
(*Plots*)


{noEatingColor, eatingColor} = Map[ColorData["BrownCyanTones"], {0, 1}]


linePlots = MapThread[
	ListPlot[#1[[All, All, 3]], 
		DataRange -> MinMax[lightLevels],
		PlotRange -> {MinMax[lightLevels], {0, 1}},
		PlotRangePadding -> {10, 0.05},
		PlotLabel -> "Evolution to extinction:\n" <> #2,
		Frame -> {True, True, False, False},
		FrameLabel -> {"Light (\[Mu]mol photons \!\(\*SuperscriptBox[\(m\), \(-2\)]\) \!\(\*SuperscriptBox[\(s\), \(-1\)]\))", If[#2 == "No trade-off", "Fraction of simulations", None]},
		PlotLegends -> If[#2 == "Attack trade-off", {"Phototroph does not eat", "Phototroph eats"}, None],
		LabelStyle -> {FontSize -> 12, FontColor -> Black},
		PlotStyle -> {noEatingColor, Directive[eatingColor, Dashed]},
		PlotMarkers -> Automatic,
		Joined -> True,
		ImageSize -> 72 * 3]&,
	{{{outcomeCountsNoEatingCE, outcomeCountsEatingCE}, {outcomeCountsNoEatingCEAR, outcomeCountsEatingCEAR}}, {"No trade-off", "Attack trade-off"}}]


Row[linePlots, Spacer[36]]
Export["phototroph_eating_increases_survival.pdf", %]



