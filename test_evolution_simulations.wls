#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Test rare mutant ecological dynamics at (retention = 1, replication = 0.1)*)


(* ::Subsection::Closed:: *)
(*Resident ecological dynamics*)


ecoEx4 = residentEcoDynamics[1, 1/10, 50, 4, 1, maxTime -> 8000, attackTradeOff -> False, conversionTradeOff -> True];


(* ::Text:: *)
(*Period of cyclical dynamics*)


estPeriod = interpolatingFunctionPeriod[(w /. ecoEx4[[1]])];


(* ::Text:: *)
(*Time to start simulating mutant dynamics*)


tmin = 4000;


(* ::Text:: *)
(*Check that the resident is neutral (s = 0) at resident equilibrium*)


((h[8000] + p[8000]) / (h[8000 - estPeriod] + p[8000 - estPeriod]) /. ecoEx4[[1]]) - 1


(* ::Subsection::Closed:: *)
(*Check that if we introduce a mutant with the same traits as the resident, it is also neutral at equilibrium*)


mutDynamics4 = rareMutantEcoDynamics[(w /. ecoEx4[[1]]),  (h /. ecoEx4[[1]]), (p /. ecoEx4[[1]]), 1, 1/10, 50, 
(*(h[tmin]/(h[tmin] + p[tmin]) /. ecoEx4[[1]])*) (h[tmin] /. ecoEx4[[1]]), (p[tmin] /. ecoEx4[[1]]),
tmin, 8000, conversionTradeOff -> True, attackTradeOff -> False];


Plot[(h2[t] /. mutDynamics4[[1]]) / (h[t] /. ecoEx4[[1]]), {t, 4000, 8000}, PlotRange -> Full]


s4 = (((h2[8000] + p2[8000]) / (h2[8000 - estPeriod] + p2[8000 - estPeriod])) /. mutDynamics4[[1]]) - 1


((h[8000] + p[8000]) / (h[8000 - estPeriod] + p[8000 - estPeriod]) /. ecoEx4[[1]]) - 1


Max[1 - Exp[-2 * s4], 0]


Plot[{(h[t] /. ecoEx4[[1]]), (h2[t] /. mutDynamics4[[1]])}, {t, 7000, 8000}]


(* ::Subsection::Closed:: *)
(*Check that a mutant with higher replication than the resident has s > 0 (and thus > 0 fixation probability)*)


mutDynamics5 = rareMutantEcoDynamics[(w /. ecoEx4[[1]]),  (h /. ecoEx4[[1]]), (p /. ecoEx4[[1]]), 1, 1/10 + 1/100, 50, 
(*(h[tmin]/(h[tmin] + p[tmin]) /. ecoEx4[[1]])*) (h[tmin] /. ecoEx4[[1]]), (p[tmin] /. ecoEx4[[1]]),
tmin, 8000, conversionTradeOff -> True, attackTradeOff -> False];


s5 = (((h2[8000] + p2[8000]) / (h2[8000 - estPeriod] + p2[8000 - estPeriod])) /. mutDynamics5[[1]]) - 1


((h2[4000 + 20 * estPeriod] + p2[4000 + 20 * estPeriod]) / (h2[4000 + 19 * estPeriod] + p2[4000 + 19 * estPeriod])) /. mutDynamics5[[1]]


Max[1 - Exp[-2 * s5], 0]


Plot[(h2[t] /. mutDynamics5[[1]]) / (h[t] /. ecoEx4[[1]]), {t, 4000, 8000}, PlotRange -> Full]


Plot[(h2[t] + p2[t] /. mutDynamics5[[1]]), {t, 4000, 8000}, PlotRange -> Full]


(* ::Subsection::Closed:: *)
(*Check that a mutant with replication lower than the resident has s < 0*)


mutDynamics6 = rareMutantEcoDynamics[(w /. ecoEx4[[1]]),  (h /. ecoEx4[[1]]), (p /. ecoEx4[[1]]), 1, 1/10 - 1/100, 50, 
(*(h[tmin]/(h[tmin] + p[tmin]) /. ecoEx4[[1]])*) (h[tmin] /. ecoEx4[[1]]), (p[tmin] /. ecoEx4[[1]]),
tmin, 8000, conversionTradeOff -> True, attackTradeOff -> False];


s6 = (((h2[8000] + p2[8000]) / (h2[8000 - estPeriod] + p2[8000 - estPeriod])) /. mutDynamics6[[1]]) - 1


((h2[4000 + 20 * estPeriod] + p2[4000 + 20 * estPeriod]) / (h2[4000 + 19 * estPeriod] + p2[4000 + 19 * estPeriod])) /. mutDynamics6[[1]]


(* ::Text:: *)
(*Estimated fixation probability*)


Max[1 - Exp[-2 * s6], 0]


Plot[(h2[t] /. mutDynamics6[[1]]) / (h[t] /. ecoEx4[[1]]), {t, 4000, 8000}, PlotRange -> Full]


Plot[(h2[t] + p2[t] /. mutDynamics6[[1]]), {t, 4000, 8000}, PlotRange -> Full]


(* ::Subsection::Closed:: *)
(*Check that a mutant with a retention probability less than the resident has s < 0*)


mutDynamics7 = rareMutantEcoDynamics[(w /. ecoEx4[[1]]),  (h /. ecoEx4[[1]]), (p /. ecoEx4[[1]]), 1  - 1/100, 1/10, 50, 
(*(h[tmin]/(h[tmin] + p[tmin]) /. ecoEx4[[1]])*) (h[tmin] /. ecoEx4[[1]]), (p[tmin] /. ecoEx4[[1]]),
tmin, 8000, conversionTradeOff -> True, attackTradeOff -> False];


s7 = ((h2[8000] + p2[8000]) / (h2[8000 - estPeriod] + p2[8000 - estPeriod]) /. mutDynamics7[[1]]) - 1


ecoEx5 = residentEcoDynamics[7/50, 0, 50, 1, 1, 1, maxTime -> 8000, attackTradeOff -> False, conversionTradeOff -> True];


Plot[{(w[t] /. ecoEx5[[1]]), (h[t] /. ecoEx5[[1]]), (p[t] /. ecoEx5[[1]])}, {t, 0, 8000}]


mutDynamics8 = rareMutantEcoDynamics[(w /. ecoEx5[[1]]),  (h /. ecoEx5[[1]]), (p /. ecoEx5[[1]]), 13/100, 0, 50, 
(*(h[tmin]/(h[tmin] + p[tmin]) /. ecoEx4[[1]])*) (h[4000] /. ecoEx5[[1]]), (p[4000] /. ecoEx5[[1]]),
4000, 8000, conversionTradeOff -> True, attackTradeOff -> False];


mutDynamics8[[2]]


Plot[{(h2[t] /. mutDynamics8[[1]]), (p2[t] /. mutDynamics8[[1]])}, {t, 7950, 8000}]


s8 = ((h2[8000] + p2[8000]) / (h2[8000 - 50] + p2[8000 - 50]) /. mutDynamics8[[1]]) - 1


(* ::Section:: *)
(*Test evolution simulation*)


(* ::Subsection::Closed:: *)
(*Single run of the simulation (no attack trade-off)*)


testRun = Monitor[simulateEvolution[50, attackTradeOff -> False, conversionTradeOff -> True, maxMutations -> 500, initialRetention -> 0, initialReplication -> 0], i];


ListLinePlot[Transpose[testRun[[1;;2]]]]


Labeled[
	ListPlot[MapThread[Style[#1, If[#2, Red, Blue]]&, testRun[[3;;4]]], 
		PlotLabel -> "Computations take longer when resident shows cycling",
		Frame -> {True, True, False, False}, FrameLabel -> {"Mutational step number", "Computational time of step (sec)"},
		PlotLegends -> PointLegend[{Red, Blue}, {"Cycling", "No cycling"}],
		LabelStyle -> {Black, "Arial", 16}, ImageSize -> Large],
	Style["Total computational time: " <> ToString[UnitConvert[Quantity[Total[Select[testRun[[3]], NumericQ]], "Seconds"], MixedUnit[{"Minutes", "Seconds"}]]],
		{Black, "Arial", 14, Italic}]]


(* ::Subsubsection:: *)
(*Test that back mutations don't occur (okay if it happens rarely around the replication = 1 region where mutations are neutral)*)


backMuts = Pick[MapThread[{#1, #2, #3}&, {Transpose[testRun[[1;;2, 1;;499]]], Transpose[testRun[[1;;2, 2;;500]]], Range[1, 499]}],
 Table[testRun[[1;;2, i]] == testRun[[1;;2, i + 2]], {i, 1, 499}]];


(* ::Subsection::Closed:: *)
(*Simulate evolution at a variety of light levels, with multiple replicates (no attack trade-off)*)


testSet = Monitor[Table[{light, simulateEvolution[light, attackTradeOff -> False, conversionTradeOff -> True, maxMutations -> 500]}, 
	{light, {25, 35, 50}}, {j, 3}], {light, j, i}];


Table[
	Show[ListLinePlot[Transpose[testSet[[i, j, 2, 1;;2]]], AspectRatio -> 1, PlotStyle -> Black,
		Frame -> True, FrameLabel -> If[i == 3, If[j == 1, {"Retention", "Replication"}, {"Rentention", ""}], If[j == 1, {"", "Replication"}, {"", ""}]],
		PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.01, PlotRangeClipping -> False,
		PlotLabel -> "Sim " <> ToString[j + (i - 1)*3] <> ", light = " <> ToString[testSet[[i, j, 1]]] <> "\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s"],
	If[(i == j == 1),
		ListPlot[{{testSet[[i, j, 2, 1;;2, 1]]}, {testSet[[i, j, 2, 1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10},
		PlotLegends -> Placed[PointLegend[{"Simulation start", "Simulation end"}, 
			LegendLabel -> Placed[Style["No attack tradeoff", Italic], Below]], {0.25, 0.85}], LabelStyle -> 14],
	
		ListPlot[{{testSet[[i, j, 2, 1;;2, 1]]}, {testSet[[i, j, 2, 1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10}]],
		LabelStyle -> {14, "Arial", Black}, ImageSize -> 300],
{i, 3}, {j, 3}] // Grid

Export["sample_no_attack_tradeoff_evo_sims.pdf", %]


testSet2 = Monitor[Table[{light, simulateEvolution[light, attackTradeOff -> True, conversionTradeOff -> True, maxMutations -> 500]}, 
	{light, {25, 35, 50}}, {j, 3}], {light, j, i}];


Export["test_sims_CE.wdx", testSet]
Export["test_sims_CEAR.wdx", testSet2]
Export["test_sim_CEAR_extra.wdx", {{80, testRun3}, {50, testRun4}}]


Table[
	Show[ListLinePlot[Transpose[testSet2[[i, j, 2, 1;;2]]], AspectRatio -> 1, PlotStyle -> Black,
		Frame -> True, FrameLabel -> If[i == 3, If[j == 1, {"Retention", "Replication"}, {"Rentention", ""}], If[j == 1, {"", "Replication"}, {"", ""}]],
		PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.01, PlotRangeClipping -> False,
		PlotLabel -> "Sim " <> ToString[j + (i - 1)*3] <> ", light = " <> ToString[testSet2[[i, j, 1]]] <> "\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s"],
	If[(i == j == 1),
		ListPlot[{{testSet2[[i, j, 2, 1;;2, 1]]}, {testSet2[[i, j, 2, 1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10},
		PlotLegends -> Placed[PointLegend[{"Simulation start", "Simulation end"}, 
			LegendLabel -> Placed[Style["Attack tradeoff", Italic], Below]], {0.25, 0.85}], LabelStyle -> 14],
	
		ListPlot[{{testSet2[[i, j, 2, 1;;2, 1]]}, {testSet2[[i, j, 2, 1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10}]],
		LabelStyle -> {14, "Arial", Black}, ImageSize -> 300],
{i, 3}, {j, 3}] // Grid

Export["sample_attack_tradeoff_evo_sims.png", %]


(* ::Subsection::Closed:: *)
(*Simulate evolution at a variety of light levels, with multiple replicates (attack trade-off)*)


testRun2 = Monitor[simulateEvolution[50, attackTradeOff -> True, conversionTradeOff -> True, maxMutations -> 500, initialRetention -> 0, initialReplication -> 0], i];


ListLinePlot[Transpose[testRun2[[1;;2]]], PlotRange -> {{0, 1}, {0, 1}}]


(* ::Subsubsection:: *)
(*Simulate a higher light level (80 umol photons/m^2/s)*)


testRun3 = Monitor[simulateEvolution[80, attackTradeOff -> True, conversionTradeOff -> True, maxMutations -> 500, initialRetention -> 0, initialReplication -> 0], i];


ListLinePlot[Transpose[testRun3[[1;;2]]]]


(* ::Subsubsection:: *)
(*Simulate starting from a different initial traits combination (replication = retention = 0.1)*)


testRun4 = Monitor[simulateEvolution[50, attackTradeOff -> True, conversionTradeOff -> True, maxMutations -> 500, initialRetention -> 1/10, initialReplication -> 1/10], i];


(* ::Subsubsection:: *)
(*Plot all results together*)


Join[Table[
	Show[ListLinePlot[Transpose[testSet2[[i, j, 2, 1;;2]]], AspectRatio -> 1, PlotStyle -> Black,
		Frame -> True, FrameLabel -> If[i == 3, If[j == 1, {"", "Replication"}, If[j == 3, {"Rentention", ""}, {"",""}]], If[j == 1, {"", "Replication"}, {"", ""}]],
		PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.01, PlotRangeClipping -> False,
		PlotLabel -> "Sim " <> ToString[j + (i - 1)*3] <> ", light = " <> ToString[testSet2[[i, j, 1]]] <> "\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s"],
	If[(i == j == 1),
		ListPlot[{{testSet2[[i, j, 2, 1;;2, 1]]}, {testSet2[[i, j, 2, 1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10},
		PlotLegends -> Placed[PointLegend[{"Simulation start", "Simulation end"}, 
			LegendLabel -> Placed[Style["Attack tradeoff", Italic], Below]], {0.25, 0.85}], LabelStyle -> 14],
	
		ListPlot[{{testSet2[[i, j, 2, 1;;2, 1]]}, {testSet2[[i, j, 2, 1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10}]],
		LabelStyle -> {14, "Arial", Black}, ImageSize -> 300],
{i, 3}, {j, 3}],

	{{Show[ListLinePlot[Transpose[testRun4[[1;;2]]], AspectRatio -> 1, PlotStyle -> Black,
		Frame -> True, FrameLabel -> {"Retention", "Replication"},
			PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.01, PlotRangeClipping -> False,
			PlotLabel -> "Sim " <> ToString[10] <> ", light = " <> ToString[50] <> "\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s"],
		
		ListPlot[{{testRun4[[1;;2, 1]]}, {testRun4[[1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10}],
			LabelStyle -> {14, "Arial", Black}, ImageSize -> 300],
			
		Show[ListLinePlot[Transpose[testRun3[[1;;2]]], AspectRatio -> 1, PlotStyle -> Black,
		Frame -> True, FrameLabel -> {"Retention", ""},
			PlotRange -> {{0, 1}, {0, 1}}, PlotRangePadding -> 0.01, PlotRangeClipping -> False,
			PlotLabel -> "Sim " <> ToString[11] <> ", light = " <> ToString[80] <> "\[Mu]mol photons/\!\(\*SuperscriptBox[\(m\), \(2\)]\)/s"],
		
		ListPlot[{{testRun3[[1;;2, 1]]}, {testRun3[[1;;2, -1]]}}, 
			PlotStyle -> {Darker[Green], Magenta}, PlotMarkers -> {Automatic, 10}],
			LabelStyle -> {14, "Arial", Black}, ImageSize -> 300]}}] // Grid

Export["sample_attack_tradeoff_evo_sims.pdf", %]
