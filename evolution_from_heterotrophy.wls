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


(* ::Section::Closed:: *)
(*Run simulations*)


(* ::Subsection::Closed:: *)
(*Light levels to investigate*)


(* ::Text:: *)
(*Light levels to simulate evolution from heterotrophy*)


simLightLevels = Table[i, {i, 20, 90, 5}];


(* ::Subsection::Closed:: *)
(*Conversion efficiency trade-off*)


(* ::Subsubsection::Closed:: *)
(*Simulate evolution*)


evolutionSimsCE = ParallelMap[{#, simulateEvolution[#, conversionTradeOff -> True, attackTradeOff -> False]}&, 
	Flatten[Table[light, {light, simLightLevels}, 50]]];


Export[NotebookDirectory[] <> "evolution_simsCE.wdx", evolutionSimsCE];


evolutionSimsCE1 = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCE_Aug20_peak_gr_method.wdx"];


evolutionSimsCE2 = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCE_Aug29_peak_gr_method_light85-90.wdx"];


evolutionSimsCE = Join[evolutionSimsCE1, evolutionSimsCE2];


gEvoSimsCE = GatherBy[SortBy[evolutionSimsCE, #[[1]]&], #[[1]]&];
simEndpointsCE = gEvoSimsCE[[All, All, 2, 1;;2, -1]];


(* ::Subsubsection::Closed:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


outcomesCE = Map[classifyOutcome, simEndpointsCE, {2}];


outcomeCountsCE = Outer[Count, outcomesCE, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsCE = outcomeCountsCE / Dimensions[outcomesCE][[2]];


Map[Total, outcomeCountsCE]


(* ::Subsubsection::Closed:: *)
(*Plot a bar chart showing how light affects the endpoint of evolution ("outcome")*)


lightLevels = DeleteDuplicates[Flatten[gEvoSimsCE[[All, All, 1]]]];


conversionTradeOffChart = BarChart[outcomeCountsCE, ChartLayout -> "Stacked", ChartLegends -> outcomeLabels, ChartStyle -> myColors, 
Frame -> {True, True, False, False},
(* Put the light levels on the x-axis *)
FrameTicks -> {Table[{i, lightLevels[[i]]}, {i, Length[lightLevels]}], Automatic},
FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "No trade-off"] (* This is the title because we're assuming everything has the CE trade-off, so we're just calling it normal *)


(* ::Subsection::Closed:: *)
(*Conversion efficiency and attack rate trade-offs*)


(* ::Subsubsection:: *)
(*Simulate evolution*)


evolutionSimsCEAR = ParallelMap[{#, simulateEvolution[#, conversionTradeOff -> True, attackTradeOff -> True]}&, 
	Flatten[Table[light, {light, simLightLevels}, 150]]];


Export["evolutionSimsCEAR_Sept2_light150_more_sims_peak_gr_method.wdx", evolutionSimsCEAR];


(* ::Text:: *)
(*We ran 50 extra simulations at light levels 50 and 55*)


evolutionSimsCEARextraSims = ParallelMap[{#, simulateEvolution[#, conversionTradeOff -> True, attackTradeOff -> True]}&, 
	Flatten[Table[light, {light, {50, 55}}, 50]]];


evolutionSimsCEAR0 = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCEAR_Aug16_peak_gr_method.wdx"];


evolutionSimsCEAR1 = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCEAR_Aug18_peak_gr_method.wdx"];


evolutionSimsCEAR2 = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCEAR_Aug20_peak_gr_method.wdx"];


evolutionSimsCEAR3 = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCEAR_Aug21_light85-90_peak_gr_method.wdx"];


evolutionSimsCEAR4 = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCEAR_Aug30_light50-55_peak_gr_method.wdx"];


evolutionSimsCEAR = Join[evolutionSimsCEAR0, evolutionSimsCEAR1, evolutionSimsCEAR2, evolutionSimsCEAR3, evolutionSimsCEAR4];


gEvoSimsCEAR = GatherBy[SortBy[evolutionSimsCEAR, #[[1]]&], #[[1]]&];
simEndpointsCEAR = gEvoSimsCEAR[[All, All, 2, 1;;2, -1]];


(* ::Subsubsection:: *)
(*Classify and count the outcomes of the simulations (endpoints of evolution)*)


outcomesCEAR = Map[classifyOutcome, simEndpointsCEAR, {2}];


outcomeCountsCEAR = Outer[Count, outcomesCEAR, outcomeLabels, 1];

(* Turn the counts into fraction of simulations *)
outcomeCountsCEAR = outcomeCountsCEAR / Map[Length, outcomesCEAR];


Map[Total, outcomeCountsCEAR]


(* ::Subsubsection::Closed:: *)
(*Investigate the simulations that ended up in unexpected places*)


(* ::Text:: *)
(*Find the simulations that weren't able to be classified (those are the ones that ended in a weird place)*)


unclassifiedSimsCEAR = Position[outcomesCEAR, "No label", {2}]


(* ::Text:: *)
(*Find their endpoints*)


Map[simEndpointsCEAR[[#[[1]], #[[2]]]]&, unclassifiedSimsCEAR]


(* ::Text:: *)
(*Plot those simulations!*)


Row[Append[
	Map[
		ListLinePlot[gEvoSimsCEAR[[#[[1]], #[[2]], 2, 1;;2]], PlotLabel -> "Light = " <> ToString[gEvoSimsCEAR[[#[[1]], #[[2]], 1]]],
			Frame -> {True, True, False, False}, FrameLabel -> {"Mutations", "Trait"}, PlotRange -> {{0, 500}, {0, 1}},
			PlotStyle -> {Blue, Orange}, LabelStyle -> {14, Black, "Arial"}, ImageSize -> 300]&, 
		unclassifiedSimsCEAR],
	LineLegend[{Blue, Orange}, {"Retention", "Replication"}, LabelStyle -> {14, Black, "Arial"}]]]


(* ::Text:: *)
(*Restart the simulations from their endpoints. They, uh, just stay there. Let's try manually introducing some mutations.*)


restartSimsCEAR = ParallelMap[
	simulateEvolution[lightLevels[[#[[1]]]], initialRetention -> simEndpointsCEAR[[#[[1]], #[[2]], 1]],
	initialReplication -> simEndpointsCEAR[[#[[1]], #[[2]], 2]], attackTradeOff -> True, conversionTradeOff -> True]&,
unclassifiedSimsCEAR];


restartSimsCEAR[[All, 1;;2, -1]]


restartSimsCEAR[[All, 1;;2, -1]] == Map[simEndpointsCEAR[[#[[1]], #[[2]]]]&, unclassifiedSimsCEAR]


Row[Append[
	Table[
		ListLinePlot[{MapThread[{#1, #2}&, {Range[0, 500], gEvoSimsCEAR[[unclassifiedSimsCEAR[[i, 1]], unclassifiedSimsCEAR[[i, 2]], 2, 1]]}],
				MapThread[{#1, #2}&, {Range[0, 500], gEvoSimsCEAR[[unclassifiedSimsCEAR[[i, 1]], unclassifiedSimsCEAR[[i, 2]], 2, 2]]}],
				MapThread[{#1, #2}&, {Range[500, 1000], restartSimsCEAR[[i, 1]]}],
				MapThread[{#1, #2}&, {Range[500, 1000], restartSimsCEAR[[i, 2]]}]}, 
			PlotLabel -> "Light = " <> ToString[gEvoSimsCEAR[[unclassifiedSimsCEAR[[i, 1]], unclassifiedSimsCEAR[[i, 2]], 1]]],
			Frame -> {True, True, False, False}, FrameLabel -> {"Mutations", "Trait"}, PlotRange -> {{0, 1000}, {0, 1}},
			PlotStyle -> {Directive[{Dashed, Blue}], Directive[{Dashed, Orange}], Blue, Orange}, LabelStyle -> {14, Black, "Arial"}, ImageSize -> 300], 
		{i, Length[restartSimsCEAR]}],
	LineLegend[{Blue, Orange}, {"Retention", "Replication"}, LabelStyle -> {14, Black, "Arial"}]]]


Transpose[restartSimsCEAR[[1, 1;;2, All;;-2]]]


MapThread[If[#1 == #2, #1]&, {Transpose[restartSimsCEAR[[1, 1;;2, 3;;All]]], Transpose[restartSimsCEAR[[1, 1;;2, All;;-3]]]}]


(* ::Text:: *)
(*I'm pretty sure the problem is that by accident we've gotten things is a slightly different part of the cycle or something. I suppose I should check.*)


residentEcoUnclassifiedCEAR = Map[
residentEcoDynamics[simEndpointsCEAR[[#[[1]], #[[2]], 1]], simEndpointsCEAR[[#[[1]], #[[2]], 2]],
	lightLevels[[#[[1]]]], 1, 1, 1, maxTime -> 8000, attackTradeOff -> True, conversionTradeOff -> True]&,
unclassifiedSimsCEAR];


(* ::Text:: *)
(*Find cycle period*)


cyclePeriods = Map[interpolatingFunctionPeriod[(h /. #)]&, residentEcoUnclassifiedCEAR[[All, 1]]];


(* ::Text:: *)
(*Plot resident ecology*)


Map[Plot[((h[t] + p[t]) /. #), {t, 7000, 8000}]&, residentEcoUnclassifiedCEAR[[All, 1]]]


(* ::Text:: *)
(*Try increasing retention*)


retentionMutants = MapThread[
	rareMutantEcoDynamics[(w /. #1), (h /. #1), (p /. #1), 
		simEndpointsCEAR[[#2[[1]], #2[[2]], 1]] + 1/100, simEndpointsCEAR[[#2[[1]], #2[[2]], 2]],
		lightLevels[[#2[[1]]]], (h[4000] /. #1), (p[4000] /. #1), 4000, 8000,
		attackTradeOff -> True, conversionTradeOff -> True]&,
	{residentEcoUnclassifiedCEAR[[All, 1]], unclassifiedSimsCEAR}];


(* ::Text:: *)
(*Well, these all look like they can invade, so what gives?*)


Map[Plot[((h2[t] + p2[t]) /. #), {t, 7000, 8000}]&, retentionMutants[[All, 1]]]


(* ::Text:: *)
(*Estimate mutant growth rate*)


estGrowthRates = MapThread[
	(((h2[8000] + p2[8000]) / (h2[8000 - #2] + p2[8000 - #2])) /. #1)&,
	{retentionMutants[[All, 1]], cyclePeriods}]


(* ::Text:: *)
(*Yep, it's a problem with estimation of the growth rate. Try starting from another time point*)


estGrowthRates2 = MapThread[
	(((h2[4000 + #2 * 20] + p2[4000 + #2 * 20]) / (h2[4000 + #2 * 19] + p2[4000 + #2 * 19])) /. #1)&,
	{retentionMutants[[All, 1]], cyclePeriods}]


estGrowthRates3 = MapThread[
	(((h2[7999] + p2[7999]) / (h2[7999 - #2] + p2[7999 - #2])) /. #1)&,
	{retentionMutants[[All, 1]], cyclePeriods}]


estGrowthRates3 = MapThread[
	(((h2[7900] + p2[7900]) / (h2[7900 - #2] + p2[7900 - #2])) /. #1)&,
	{retentionMutants[[All, 1]], cyclePeriods}]


(* ::Text:: *)
(*Ugh. I guess there's something going on at t = 8000?*)


estGrowthRatesLots = MapThread[
	Table[(((h2[endpt] + p2[endpt]) / (h2[endpt - #2] + p2[endpt - #2])) /. #1),
		{endpt, {4100, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 7900, 7950, 7990, 8000}}]&,
	{retentionMutants[[All, 1]], cyclePeriods}];


TableForm[Map[If[# < 1, Style[#, Bold], #]&, Transpose[estGrowthRatesLots], {2}], 
TableHeadings -> {{4100, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 7900, 7950, 7990, 8000}, {"Sim 1", "Sim 2", "Sim 3"}}]


estGrowthRatesLots2 = MapThread[
	Table[(((h2[endpt] + p2[endpt]) / (h2[endpt - 5 * #2] + p2[endpt - 5 * #2])) /. #1),
		{endpt, {4300, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 7900, 7950, 7990, 8000}}]&,
	{retentionMutants[[All, 1]], cyclePeriods}];


TableForm[Map[If[# < 1, Style[#, Bold], #]&, Transpose[estGrowthRatesLots2], {2}], 
TableHeadings -> {{4300, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 7900, 7950, 7990, 8000}, {"Sim 1", "Sim 2", "Sim 3"}}]


estGrowthRatesLots3 = MapThread[
	Table[(((h2[endpt] + p2[endpt]) / (h2[endpt - 9 * #2] + p2[endpt - 9 * #2])) /. #1),
		{endpt, {4500, 5000, 5500, 6000, 6500, 7000, 7500, 7900, 7950, 7990, 8000}}]&,
	{retentionMutants[[All, 1]], cyclePeriods}];


TableForm[Map[If[# < 1, Style[#, Bold], #]&, Transpose[estGrowthRatesLots3], {2}], 
TableHeadings -> {{4300, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 7900, 7950, 7990, 8000}, {"Sim 1", "Sim 2", "Sim 3"}}]


estGrowthRatesLots4 = MapThread[
	Table[(((h2[endpt] + p2[endpt]) / (h2[endpt - 75 * #2] + p2[endpt - 75 * #2])) /. #1),
		{endpt, {7900, 7950, 7990, 8000}}]&,
	{retentionMutants[[All, 1]], cyclePeriods}];


TableForm[Map[If[# < 1, Style[#, Bold], #]&, Transpose[estGrowthRatesLots4], {2}], 
TableHeadings -> {{7900, 7950, 7990, 8000}, {"Sim 1", "Sim 2", "Sim 3"}}]


estGrowthRatesLots5 = MapThread[
	Table[(((h2[endpt] + p2[endpt]) / (h2[endpt - #2] + p2[endpt - #2])) /. #1),
		{endpt, (4000 + #2 * {10, 15, 20, 30, 40, 50})}]&,
	{retentionMutants[[All, 1]], cyclePeriods}];


TableForm[Map[If[# < 1, Style[#, Bold], #]&, Transpose[estGrowthRatesLots5], {2}], 
TableHeadings -> {{"10x", "15x", "20x", "30x", "40x", "50x"}, {"Sim 1", "Sim 2", "Sim 3"}}]


(* ::Subsubsection::Closed:: *)
(*Plot a bar chart showing how light affects the endpoint of evolution ("outcome")*)


lightLevels = DeleteDuplicates[Flatten[gEvoSimsCEAR[[All, All, 1]]]];


conversionAndAttackTradeOffChart = BarChart[outcomeCountsCEAR, ChartLayout -> "Stacked", ChartLegends -> outcomeLabels, ChartStyle -> myColors, 
Frame -> {True, True, False, False},
(* Put the light levels on the x-axis *)
FrameTicks -> {Table[{i, lightLevels[[i]]}, {i, Length[lightLevels]}], Automatic},
FrameLabel -> {"Surface light (\[Mu]mol photons/m^2/s)", "Fraction of simulations"},
PlotLabel -> "Attack rate trade-off"] (* This is the title because we're assuming everything has the CE trade-off, so we're just calling it normal *)


(* ::Section::Closed:: *)
(*Plot endpoint of evolution vs light*)


(* ::Subsection::Closed:: *)
(*Legend*)


(* ::Code::Initialization::Bold:: *)
linePlotLegend = Framed[
	LineLegend[{myColors[[1]], myColors[[2]], Black}, 
		{"Obligately klepto. mixotroph", "Facultatively klepto. phototroph", "Extinction"},
		(*LegendLabel -> Placed["Endpoint of evolution", Above],*)
		LegendMarkers -> Automatic, (*LegendMarkerSize -> {15, 6}, *)
		LegendMargins -> 0, Spacings -> 0.25,
		LabelStyle -> {Black, 12, "Arial"}],
	RoundingRadius -> 10, FrameStyle -> Thickness[0.5]];


(* ::Subsection::Closed:: *)
(*Conversion trade-off*)


outcomesAndLightCE = MapThread[{#1, #2}&, {{lightLevels, lightLevels, lightLevels}, Transpose[outcomeCountsCE]}, 2];


linePlotCE = ListLinePlot[outcomesAndLightCE, PlotStyle -> {myColors[[1]], myColors[[2]], Black},
PlotMarkers -> Automatic,
PlotLegends -> Placed[linePlotLegend, {0.55, 0.5}],
Frame -> True, FrameLabel -> {"Surface light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\))", "Fraction of simulations"},
PlotLabel -> "(A) No tradeoff", 
ImageSize -> 4.4 * 72, LabelStyle -> {Black, 12, "Arial"}]


(* ::Subsection::Closed:: *)
(*Attack and conversion trade-offs*)


outcomesAndLightCEAR = MapThread[{#1, #2}&, {{lightLevels, lightLevels, lightLevels}, Transpose[outcomeCountsCEAR]}, 2];


linePlotCEAR = ListLinePlot[outcomesAndLightCEAR, PlotStyle -> {myColors[[1]], myColors[[2]], Black},
PlotMarkers -> Automatic,
Frame -> True, FrameLabel -> {"Surface light (\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s)", ""},
PlotLabel -> "(B) Attack rate tradeoff",
LabelStyle -> {Black, 12, "Arial"}, ImageSize -> 4.4 * 72]


(* ::Subsection::Closed:: *)
(*Plot both trade-off variations together*)


(* Add letter labels to plots *)
linePlots = {linePlotCE, linePlotCEAR};

(* Put the plots in a column with a single, shared legend *)
Row[linePlots, Spacer[1]]

(* Export as pdf *)
Export["evo_from_heterotrophy.pdf", %];


(* ::Section::Closed:: *)
(*Plot evolutionary trajectories*)


(* ::Text:: *)
(*A list of letters for labeling things that might have more than 26 parts*)


letterList = Flatten[Outer[#1 <> #2&, Prepend[ToUpperCase[Alphabet[]], ""], ToUpperCase[Alphabet[]]]];


(* ::Subsection::Closed:: *)
(*Conversion trade-off*)


(* ::Subsubsection::Closed:: *)
(*Color the evolutionary trajectories based on the ESS the simulations end up at*)


lineColorsCE = 
Map[If[# == "Facultatively kleptoplastidic phototroph", myColors[[2]], 
	If[# == "Extinction", Black, 
	Green]]&, 
outcomesCE, {2}];


(* ::Subsubsection::Closed:: *)
(*Plot the evolutionary trajectories in retention-replication space*)


trajectoryPlotsCE = Table[
ListLinePlot[Transpose[gEvoSimsCE[[i, All, 2, 1;;2]], 2 <-> 3],  
PlotStyle -> Map[Directive[Thin, #]&, lineColorsCE[[i]]],
PlotLabel -> "(" <> letterList[[i - 1]] <> ") Light = " <> ToString[lightLevels[[i]]] (*<> " \[Mu]mol photons/m^2/day"*),
Frame -> True,(* FrameLabel -> {"Retention prob. (\[Phi])", "Replication prob. (\[Eta])"},*)
PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.025,
LabelStyle -> {11, Black, "Arial"},
ImageSize -> 72 * 2.1, AspectRatio -> 1], 
{i, 2, Length[lightLevels]}];


(* ::Subsubsection::Closed:: *)
(*Add axis labels to the plots*)


xlabel = Style["Retention prob. (\[Phi])", FontSize -> 12, FontFamily -> "Arial", FontColor -> Black];
ylabel = Rotate[Pane[
			Style["Replication prob. (\[Eta])", FontSize -> 12, FontFamily -> "Arial", FontColor -> Black], 
		ImageSize -> 72 * 2.1, Alignment -> -0.33], 
	90 Degree];


(* ::Text:: *)
(*Put the plots in a grid with axis labels on the left and bottom*)


trajPlotRowsCE = Partition[trajectoryPlotsCE, UpTo[4]];
trajPlotRowsCE = Map[Prepend[#, ylabel]&, trajPlotRowsCE]; (* ylabels at the left of each row *)
trajPlotRowsCE = Insert[trajPlotRowsCE, {"", "", "", xlabel, xlabel}, -2]; (* xlabels at the bottom of the right two rows *)
trajPlotRowsCE = Append[trajPlotRowsCE, {"", xlabel, xlabel}]; (* xlabels at the bottom of the left two rows *)


(* ::Subsubsection::Closed:: *)
(*Add a legend to the plots*)


trajPlotLegend = SwatchLegend[{myColors[[2]], Black}, {"Facultatively kleptoplastidic phototroph", "Extinction"},
LegendLabel -> Style["Endpoint of evolution", 14],
LabelStyle -> {12, "Arial", Black}];

trajPlotLegend = Framed[trajPlotLegend, RoundingRadius -> 5];
trajPlotLegend = Pane[trajPlotLegend, ImageSize -> {72 * 4.2, 72 * 2.1}, Alignment -> Center];


trajPlotRowsCE[[-2]] = Join[trajPlotRowsCE[[-2]], {trajPlotLegend, SpanFromLeft}]; (* add the legend to the list of plots *)


(* ::Subsubsection::Closed:: *)
(*Show the plots and legend together in a grid*)


Grid[trajPlotRowsCE, Alignment -> {Center, Top}, Spacings -> {1, {1, 1, 1, 0, 0.25, 0}}];
Export["trajectories_ce.pdf", %]


(* ::Subsection::Closed:: *)
(*Attack & conversion trade-offs*)


(* ::Subsubsection::Closed:: *)
(*Color the evolutionary trajectories based on the ESS the simulations end up at*)


lineColorsCEAR = 
Map[If[# == "Obligately kleptoplastidic mixotroph", myColors[[1]], 
	If[# == "Facultatively kleptoplastidic phototroph", myColors[[2]], 
	Black]]&, 
outcomesCEAR, {2}];


(* ::Subsubsection::Closed:: *)
(*Plot the evolutionary trajectories in retention-replication space*)


trajectoryPlotsCEAR = Table[
ListLinePlot[Transpose[gEvoSimsCEAR[[i, All, 2, 1;;2]], 2 <-> 3],  
PlotStyle -> Map[Directive[Thin, #]&, lineColorsCEAR[[i]]],
PlotLabel -> "(" <> letterList[[i - 1]] <> ") Light = " <> ToString[lightLevels[[i]]] (*<> " \[Mu]mol photons/m^2/day"*),
Frame -> True,(* FrameLabel -> {"Retention prob. (\[Phi])", "Replication prob. (\[Eta])"},*)
PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.025,
LabelStyle -> {11, Black, "Arial"},
ImageSize -> 72 * 2.1, AspectRatio -> 1], 
{i, 2, Length[lightLevels]}];


(* ::Subsubsection::Closed:: *)
(*Add axis labels to the plots*)


xlabel = Style["Retention prob. (\[Phi])", FontSize -> 12, FontFamily -> "Arial", FontColor -> Black];
ylabel = Rotate[Pane[
			Style["Replication prob. (\[Eta])", FontSize -> 12, FontFamily -> "Arial", FontColor -> Black], 
		ImageSize -> 72 * 2.1, Alignment -> -0.33], 
	90 Degree];


(* ::Text:: *)
(*Put the plots in a grid with axis labels on the left and bottom*)


trajPlotRowsCEAR = Partition[trajectoryPlotsCEAR, UpTo[4]];
trajPlotRowsCEAR = Map[Prepend[#, ylabel]&, trajPlotRowsCEAR]; (* ylabels at the left of each row *)
trajPlotRowsCEAR = Insert[trajPlotRowsCEAR, {"", "", "", xlabel, xlabel}, -2]; (* xlabels at the bottom of the right two rows *)
trajPlotRowsCEAR = Append[trajPlotRowsCEAR, {"", xlabel, xlabel}]; (* xlabels at the bottom of the left two rows *)


(* ::Subsubsection::Closed:: *)
(*Add a legend to the plots*)


trajPlotLegend = SwatchLegend[myColors[[1;;2]], {"Obligately kleptoplastidic mixotroph", "Facultatively kleptoplastidic phototroph"},
LegendLabel -> Style["Endpoint of evolution", 14],
LabelStyle -> {12, "Arial", Black}];

trajPlotLegend = Framed[trajPlotLegend, RoundingRadius -> 5];
trajPlotLegend = Pane[trajPlotLegend, ImageSize -> {72 * 4.2, 72 * 2.1}, Alignment -> Center];


trajPlotRowsCEAR[[-2]] = Join[trajPlotRowsCEAR[[-2]], {trajPlotLegend, SpanFromLeft}]; (* add the legend to the list of plots *)


(* ::Text:: *)
(*Show the plots and legend together in a grid*)


Grid[trajPlotRowsCEAR, Alignment -> {Center, Top}, Spacings -> {1, {1, 1, 1, 0, 0.25, 0}}];
Export["trajectories_cear.pdf", %];


(* ::Section::Closed:: *)
(*Conversion trade-off: investigate drift in retention probability*)


(* ::Text:: *)
(*We want to check that drift only happens when replication = 1*)


(* ::Subsection::Closed:: *)
(*Plot retention vs. time for each simulation*)


MapThread[ListLinePlot[#1[[All, 2, 1]], PlotRange -> {{0, 500}, {0, 1}},
Frame -> True, FrameLabel -> {"Time (mutations)", "Retention prob."},
PlotLabel -> "Light = " <> ToString[#2]]&, 
{gEvoSimsCE[[2;;All]], lightLevels[[2;;All]]}]


(* ::Text:: *)
(*Most drift looks like it occurs at light = 25 and light = 55*)


MapThread[ListLinePlot[#1[[All, 2, 1]], PlotRange -> {{0, 500}, {0, 1}},
Frame -> True, FrameLabel -> {"Time (mutations)", "Retention prob."},
PlotLabel -> "Light = " <> ToString[#2],
ImageSize -> 400]&, 
{gEvoSimsCE[[{2, 8}]], lightLevels[[{2, 8}]]}]


(* ::Subsection::Closed:: *)
(*Take a look at trait values that are revisited (these should occur in regions where retention is drifting)*)


talliedTraitVisits = Map[Tally, Transpose[gEvoSimsCE[[2;;All, All, 2, 1;;2]], 3 <-> 4], {2}];
repeatTraitVisits = Map[Select[#[[2]] > 1&], talliedTraitVisits, {2}][[All, All, All, 1]];


(* ::Text:: *)
(*All repeated trait values occur where replication = 1*)


Sort[DeleteDuplicates[Flatten[repeatTraitVisits, 2]]]


(* ::Subsection::Closed:: *)
(*Additionally investigate individual simulations to see if we can visually observe drift only after replication reaches 1*)


(* ::Text:: *)
(*Make a marker that we can plot to show the time point at which replication = 1*)


replication1Marker = Map[If[NumberQ[#], {{#, 0}, {#, 1}}, {{0, 0}}]&,
	Map[FirstPosition[#, 1]&, gEvoSimsCE[[All, All, 2, 2]], {2}], {3}];


(* ::Text:: *)
(*Set the light level of the simulations to plot (I looked at all simulated values between 25 and 60)*)


investigateLight = 25; (* Change this to change the light level of the simulations plotted below *)
lightListPosition = FirstPosition[lightLevels, investigateLight][[1]];


(* ::Text:: *)
(*Legend for plots below*)


LineLegend[{Red, Blue, Black}, {"Retention", "Replication", "Point where replication reaches 1"}]


(* ::Text:: *)
(*Plot simulations: we can see that drift in retention occurs only after replication = 1*)


MapThread[Show[
	ListLinePlot[#1, PlotStyle -> {Red, Blue}], (* traits vs. time *)
	ListLinePlot[#2, PlotStyle -> Black], (* marker for where replication reaches 1*)
	Frame -> True, FrameLabel -> {"Time (mutations)", "Trait"}]&, 
{gEvoSimsCE[[lightListPosition, All, 2, 1;;2]], replication1Marker[[lightListPosition]]}]


(* ::Section::Closed:: *)
(*Example ecology and evolution plots*)


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Ecology plots*)


(* ::Subsubsection::Closed:: *)
(*Run the simulation*)


(* ::Text:: *)
(*First simulate a heterotroph alone (resident) to get the initialization for a simulation with mutants*)


{residentEcology, extinction} = residentEcoDynamics[0, 0, 35, 4, 1, 1, conversionTradeOff -> True, attackTradeOff -> False]


Plot[{w[t]/.residentEcology, h[t]/.residentEcology, p[t]/.residentEcology}, {t, 0, 300}, PlotRange -> All]


{avgPrey, avgHeterotroph, avgPhototroph} = Table[NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000), {stateVar, {w[t], h[t], p[t]}}];


(* ::Text:: *)
(*Simulate a resident heterotrophic population being invaded by a mutant that can sometimes replicate the plastid*)


{withMutantEcology, extinction} = mutantAndResidentEcoDynamics[0, 0, 1/100, 0, 35, 
	avgPrey, (99/100) * avgHeterotroph, (99/100) * avgPhototroph,  (1/100) * avgHeterotroph, (1/100) * avgPhototroph, 
	conversionTradeOff -> True, attackTradeOff -> False]


(* ::Subsubsection::Closed:: *)
(*Plot prey ecological dynamics*)


(* ::Text:: *)
(*x-axis plot ticks*)


ecologicalTimePlotTicks = {Charting`FindScaledTicks[Charting`FindTicks[{0, 1},{0, 1}],{Identity,Identity}][0, 4000, {3, 9}], Automatic}; 
	(* replace the above with Charting`ScaledTicks["Linear"][0, 1, {3, 9}] for Mathematica v < 13 *)


preyPlot = Plot[w[t-100] /. withMutantEcology, {t, 100, 4100}, PlotRange -> {{0, 4000}, {0, 20}}, PlotStyle -> colorW, 
	Frame -> True,  PlotRangePadding -> 0,
	GridLines -> {{{100, Directive[Black, Dashed]}}, None},
	FrameLabel -> {"Time (days)", "Prey (cells/\!\(\*SuperscriptBox[\(m\), \(2\)]\))"}, 
	FrameTicks -> {Automatic, ecologicalTimePlotTicks},
	ImageSize -> 72 * 2.8, AspectRatio -> 3.25/4,
	LabelStyle -> {12, Black, "Arial"}];


preyPlot = Show[{preyPlot, Plot[w[t+3900] /. residentEcology, {t, 0, 100}, PlotStyle -> colorW], Graphics[{Gray, Opacity[0.2], Rectangle[{3000, 0}, {4000, 20}]}]}]


(* ::Subsubsection::Closed:: *)
(*Plot mixotroph ecological dynamics*)


(* ::Text:: *)
(*y-axis plot ticks*)


consumerPlotTicks = {Charting`FindScaledTicks[Charting`FindTicks[{0, 1},{0, 1}],{Identity,Identity}][0, 2, {4, 6}], Automatic};


(* ::Text:: *)
(*This code plots the total mixotroph population (combined phototrophs and heterotrophs)*)


mixotrophPlot = Plot[{(h[t - 100] + p[t - 100]) /. withMutantEcology // Chop, (h2[t - 100] + p2[t - 100]) /. withMutantEcology // Chop}, {t, 100, 4100}, 
	PlotStyle -> {Directive[Thick, colorM], Directive[Thick, Dashed, colorM]},
	GridLines -> {{{100, Directive[Black, Dashed]}}, None},
	PlotRange -> {{0, 4000}, {0, 2}}, PlotRangePadding -> 0,
	Frame -> True, FrameTicks -> {consumerPlotTicks, ecologicalTimePlotTicks},
	FrameLabel -> {{"Consumer (cells/\!\(\*SuperscriptBox[\(m\), \(2\)]\))", None}, {"Time (days)", None}},
	LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.8, AspectRatio -> 3.25/4];


mixotrophPlot = Show[{mixotrophPlot, 
		Plot[(h[t + 3900] + p[t + 3900]) /. residentEcology // Chop, {t, 0, 100}, 
			PlotStyle -> Directive[Thick, colorM]],
			Graphics[{Gray, Opacity[0.2], Rectangle[{3000, 0}, {4000, 20}]}]}]


(* ::Subsubsection::Closed:: *)
(*Plot the proportion of the population that is the mutant*)


mutantProportionPlot = Plot[((h2[t - 100] + p2[t - 100])/ (h[t - 100] + p[t - 100] + h2[t - 100] + p2[t - 100]))/. withMutantEcology, {t, 100, 4100}, 
	PlotRange -> {{0, 4000}, {0, 1}},  PlotRangePadding -> 0,
	PlotStyle -> Black, 
	Frame -> True, FrameTicks -> {Automatic, ecologicalTimePlotTicks},
	FrameLabel -> {"Time (days)", "Proportion mutants"},
	GridLines -> {{{100, Directive[Black, Dashed]}}, None},
	LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.8, AspectRatio -> 3.25/4];


mutantProportionPlot = Show[{mutantProportionPlot, Plot[0, {t, 0, 100}, PlotStyle -> Black], 
		Graphics[{Gray, Opacity[0.2], Rectangle[{3000, 0}, {4000, 20}]}]}]


(* ::Subsubsection::Closed:: *)
(*Combine plots into a list*)


ecologyPlots = {Labeled[Show[preyPlot, LabelStyle -> {Black, "Arial", 12}], "(A) Prey ecology", Top, Alignment -> {Left, Baseline}],
Labeled[Show[mixotrophPlot, LabelStyle -> {Black, "Arial", 12}],  "(B) Mixotroph ecology", Top, Alignment -> Left],
Labeled[Show[mutantProportionPlot, LabelStyle -> {Black, "Arial", 12}], "(C) Mutant spread", Top, Alignment -> Left]}


(* ::Subsection::Closed:: *)
(*Evolution plots*)


evolutionSimsCEexample = Import[NotebookDirectory[] <> "Saved_analyses/evolutionSimsCE_Aug20_peak_gr_method.wdx"];


gEvoSimsCEexample = GatherBy[SortBy[evolutionSimsCEexample, #[[1]]&], #[[1]]&];


exampleSim = gEvoSimsCEexample[[4, 3]];


(* ::Subsubsection::Closed:: *)
(*Plot the plastid retention probability over evolutionary time*)


retentionPlot = ListLinePlot[exampleSim[[2, 1]], 
	DataRange -> {0, Length[exampleSim[[2, 1]]] - 1},
	PlotStyle -> colorRetention, PlotRange -> {{0, 500}, {0, 1}}, PlotRangePadding -> 0,
	Frame -> True, FrameLabel -> {"Time (mutations)", "Retention"}, 
	LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.8, AspectRatio -> 3.25/4]


(* ::Subsubsection::Closed:: *)
(*Plot the plastid replication probability over evolutionary time*)


replicationPlot = ListLinePlot[exampleSim[[2, 2]], 
	DataRange -> {0, Length[exampleSim[[2, 2]]] - 1},
	PlotStyle -> colorReplication, PlotRange -> {{0, 500}, {0, 1}}, PlotRangePadding -> 0,
	Frame -> True, FrameLabel -> {"Time (mutations)", "Replication"}, 
	LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.8, AspectRatio -> 3.25/4]


(* ::Subsubsection::Closed:: *)
(*Plot the evolutionary trajectory in trait space so we can see what all the different (retention, replication) combinations were*)


traitSpacePlot = ListLinePlot[Transpose[exampleSim[[2, 1;;2]]], 
	PlotStyle -> Black, PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0,
	Frame -> True, FrameLabel -> {"Retention", "Replication"}, 
	LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.8, AspectRatio -> 3.25/4]


(* ::Subsubsection::Closed:: *)
(*Find and plot the "equilibrium" population sizes over evolutionary time*)


averagePopSizes = ParallelTable[Module[{residentEcology, extinction, retention = exampleSim[[2, 1, i]], replication = exampleSim[[2, 1, i]]}, 
	{residentEcology, extinction} = residentEcoDynamics[retention, replication, 35, 4, 1, 1, 
		conversionTradeOff -> True, attackTradeOff -> False, PrecisionGoal -> 6, AccuracyGoal -> 6];
	Table[
		If[extinction, {"Unknown", 0},
			NIntegrate[stateVar/. residentEcology, {t, 3000, 4000}]/(4000 - 3000)], 
	{stateVar, {w[t], (h[t] + p[t])}}]],
{i, Length[exampleSim[[2, 1]]]}];


averagePopPlot = ListLinePlot[{averagePopSizes[[All, 1]], averagePopSizes[[All, 2]]}, 
DataRange -> {0, Length[averagePopSizes] - 1},
PlotStyle -> {colorW, Directive[Thick, colorM]},
PlotRange -> {{0, 500}, {0, 17}}, PlotRangePadding -> 0,
Frame -> True, FrameLabel -> {"Time (mutations)", "Pop. size (cells/\!\(\*SuperscriptBox[\(m\), \(2\)]\))"},
LabelStyle -> {12, Black, "Arial"},
	ImageSize -> 72 * 2.8, AspectRatio -> 3.25/4]


(* ::Text:: *)
(*Prey go extinct when mutation 111 is introduced*)


FirstPosition[Chop[averagePopSizes[[All, 1]], 10^-6], 0] - 1 (* subtract 1 so we don't count the initial resident as one of the mutants *)


(* ::Subsubsection::Closed:: *)
(*Combine plots into a list*)


evolutionPlots = {Labeled[Show[averagePopPlot, LabelStyle -> {Black, "Arial", 12}], "(D) Average population size", Top, Alignment -> Left],
Labeled[Show[replicationPlot, LabelStyle -> {Black, "Arial", 12}], "(E) Replication evolution", Top, Alignment -> Left],
Labeled[Show[retentionPlot, LabelStyle -> {Black, "Arial", 12}], "(F) Retention evolution", Top, Alignment -> Left],
Labeled[Show[traitSpacePlot, LabelStyle -> {Black, "Arial", 12}], "(G) Evo. trajectory in trait space", Top, Alignment -> Left]}


(* ::Subsection::Closed:: *)
(*Combine ecology and evolution plots into the figure*)


myArrow = Graphics[{Thick, Line[{{1/4, 1/8}, {1/4, 3/8}, {5/8, 0}, {1/4, -3/8}, {1/4, -1/8}}],
   Line[{{1/4, 1/8}, {-1/2, 1/8}, {-1/2, -1/8}, {1/4, -1/8}}]}];


mixotrophPlotLegend = Framed[ (* put a frame around the legend *)
	LineLegend[{colorW, Directive[Thick, colorM], Directive[Thick, Dashed, colorM]}, {"Prey", "Resident consumer", "Mutant consumer"}, LegendLabel -> "Population",
		LabelStyle -> {12, "Arial", Black}], 
(* frame options *)
RoundingRadius -> 10];


Grid[Transpose[{Prepend[ecologyPlots, Style["Ecological Dynamics", FontFamily -> "Arial", FontColor -> Black, FontSize -> 12]],
Prepend[evolutionPlots[[1;;3]], Style["Evolutionary Dynamics", FontFamily -> "Arial", FontColor -> Black, FontSize -> 12]], 
{"", mixotrophPlotLegend, evolutionPlots[[4]], ""}}],
Alignment -> {Center, Center}, Spacings -> {2.5, 2}]

Export["example_timeseries.pdf", %]


(* ::Text:: *)
(*The final step is to open this in an image editor and add arrows: between ecology and evolution columns, and pointing from E and F to G. I also moved the legend a bit. This is doable in Mathematica, but sometimes it's not worth the annoyance.*)
