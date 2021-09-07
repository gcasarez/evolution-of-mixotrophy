#!/usr/bin/env wolframscript
(* ::Package:: *)

wDot[w_,hW_,lW_,a_,pW_,kappa_,iIn_,gammaA_,eta_,phi_,eta2_,phi2_,sigmaA_,p_,h_, a2_,h2_,p2_]:= 
(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-(a(1-sigmaA (phi eta))) w h-(gammaA a) w p-(a2(1-sigmaA (phi2 eta2))) w h2 - (gammaA a2) w p2 
hDot[eta_,c_,gammaC_,kappa_,gammaP_,pW_,sigmaP_,p_,hP_,iIn_,phi_,a_,w_,h_,gammaA_,sigmaA_,lH_,dMax_,sigmaD_]:=
(1-eta) (c/(gammaC c))((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi (a(1-sigmaA (phi eta))) w h + (1-phi)c (a(1-sigmaA (phi eta))) w h+ (1-eta) c (gammaA a) w p -lH h + (dMax (1-sigmaD eta))p
pDot[eta_,gammaP_,sigmaP_,pW_,p_,kappa_,hP_,iIn_,phi_,a_,w_,h_,gammaA_,sigmaA_,gammaC_,sigmaC_,c_,dMax_,sigmaD_,lP_] :=
eta ((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi (a(1-sigmaA (phi eta))) w h +eta (gammaA a) (gammaC c) p w - (dMax (1-sigmaD eta)) p - lP p
h2Dot[eta2_,c_,gammaC_,gammaP_,sigmaP_,pW_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,gammaA_,sigmaA_,lH2_,dMax_,sigmaD_] :=
(1-eta2) (c/(gammaC c))((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 (a2(1-sigmaA (phi2 eta2))) w h2 + (1-phi2) c (a2(1-sigmaA (phi2 eta2))) w h2 +(1-eta2) c (gammaA a2) w p2 - lH2 h2+ (dMax(1-sigmaD eta2)) p2
p2Dot[eta2_,gammaP_,sigmaP_,pW_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,gammaA_,sigmaA_,c_,gammaC_,sigmaC_,dMax_,sigmaD_,lP2_] :=
eta2 ((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 (a2(1-sigmaA (phi2 eta2))) w h2+eta2 (gammaA a2) (gammaC c) p2 w-(dMax(1-sigmaD eta2)) p2-lP2 p2
kappa[kW_, kH_,kP_,kH2_,kP2_,w_,h_,p_,h2_,p2_]:=
kW w+kH h+ kP p+kH2 h2+kP2 p2 
kappaNew[kW_, kH_,kP_,w_,h_,p_]:=
kW w+kH h+ kP p


pW=3
gammaP=2.9/3
sigmaP=1
kW=0.1
kH=0.05
kP=0.15
kH2=0.015
kP2=0.15
hW=30
hP=30
hP2=30
lW=0.5
lH=0.1
lP=0.4
lH2=0.1
lP2=0.4
dMax=0.3
sigmaD=0
a=0.15
a2=0.15
gammaA=0.7
sigmaA=1
phi=0
phi2=0
c=0.1
gammaC=0.7
eta=0
eta2=0
iIn=9


myJacobian[w_,hW_,lW_,a_,pW_,kappa_,iIn_,gammaA_,eta_,phi_,eta2_,phi2_,sigmaA_,p_,h_, a2_,h2_,p2_,c_,gammaC_,gammaP_,sigmaP_,hP_,lH_,dMax_,sigmaD_,sigmaC_,lP_]:=
 {{D[wDot[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,h2,p2], w], 
D[wDot[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,h2,p2], h],
D[wDot[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,h2,p2], p]}, 
{D[hDot[eta,c,gammaC,kappa,gammaP,pW,sigmaP,p,hP,iIn,phi,a,w,h,gammaA,sigmaA,lH,dMax,sigmaD], w], 
D[hDot[eta,c,gammaC,kappa,gammaP,pW,sigmaP,p,hP,iIn,phi,a,w,h,gammaA,sigmaA,lH,dMax,sigmaD], h],
D[hDot[eta,c,gammaC,kappa,gammaP,pW,sigmaP,p,hP,iIn,phi,a,w,h,gammaA,sigmaA,lH,dMax,sigmaD],p]},
{D[pDot[eta,gammaP,sigmaP,pW,p,kappa,hP,iIn,phi,a,w,h,gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],w],
D[pDot[eta,gammaP,sigmaP,pW,p,kappa,hP,iIn,phi,a,w,h,gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],h],
D[pDot[eta,gammaP,sigmaP,pW,p,kappa,hP,iIn,phi,a,w,h,gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],p]}}


negativeEigQ[jacobian_] :=
AllTrue[Eigenvalues[jacobian], Negative]


stableEquil[eta2_,c_,gammaC_,gammaP_,sigmaP_,pW_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,gammaA_,sigmaA_,lH2_,dMax_,sigmaD_,sigmaC_,lP2_, precision_: MachinePrecision] := 
Module[{equilibria},
  (* Find all equilibria *)
  equilibria = NSolve[{h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2,kappa,hP2,iIn,phi2,a2,w,h2,gammaA,sigmaA,lH2,dMax,sigmaD]==0,
  p2Dot[eta2,gammaP,sigmaP,pW,p2,kappa,hP2,iIn,phi2,a2,w,h2,gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2]==0}, {h2, p2}, Reals, WorkingPrecision -> precision];
  
  (* Select the stable equilibria *)
  equilibria = Select[equilibria, negativeEigQ[myJacobian[w/.#,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p/.#,h/.#, a2,h2,p2,c,gammaC,gammaP,sigmaP,hP,lH,dMax,sigmaD,sigmaC,lP]]&];
  
  (* Have equilibria be returned as numbers instead of iP \[Rule] number *)
  equilibria = Map[{w/.#, h/.#,p/.#}&, equilibria];
  
  Return[equilibria]
]


stableEquil[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,0,0,c,gammaC,gammaP,sigmaP,hP,lH,dMax,sigmaD,0,lP,lH2,lP2,hP2,kW,kH,kP,kH2,kP2]


(* ::Text:: *)
(*Maybe another way to look for equilibrium??*)


initialValueFunction=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==5,h[0]==5,p[0]==5, h2[0]==0,p2[0]==0},{w,h,p,h2,p2},{t,0,3000}]


wInit=w[2000]/.initialValueFunction[[1]]


hInit=h[2000]/.initialValueFunction[[1]]


pInit=Abs[p[2000]/.initialValueFunction[[1]]]


stableEquilList[traitList_, d_, precision_:MachinePrecision] := 
Module[{equil}, 

(* Make a table of stable equilibria *)
equil=Table[equil[trait[[1]], trait[[2]], d, precision], {trait, traitList}];

(* Sort the equilibria for each pair of transmission rates *)
equil = Map[Sort, equil];

(* Add the trait values that produced the equilibrium/equilibria *)
equil = Table[{traitList[[i, 1]], traitList[[i, 2]], equil[[i]]}, {i, Length[traitList]}];

Return[equil]]


mutX[w_,h_,p_,h2_, p2_, phi2_, eta2_, phi_, eta_] :=
{{-phi2 (a2(1-sigmaA (phi2 eta2))) w h2 + (1-phi2) c (a2(1-sigmaA (phi2 eta2))) w h2 - lH2 h2,
 (1-eta2) (c/(gammaC c))((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappaNew[kW, kH,kP,w,h,p]) Log[(hP2+iIn)/(hP2+iIn Exp[-kappaNew[kW, kH,kP,w,h,p]])]  +(1-eta2) c (gammaA a2) w p2 + (dMax(1-sigmaD eta2)) p2}, 
{phi2 (a2(1-sigmaA (phi2 eta2))) w h2, 
 eta2 ((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappaNew[kW, kH,kP,w,h,p]) Log[(hP2+iIn)/(hP2+iIn Exp[-kappaNew[kW, kH,kP,w,h,p]])]+eta2 (gammaA a2) (gammaC c) p2 w-(dMax(1-sigmaD eta2)) p2-lP2 p2}}


wInit


Eigenvalues[mutX[wInit,hInit,pInit,hInit,pInit,0.2,0.2,0.1,0.1],1]


mutGrowth[w_,h_, p_, phimut_, etamut_, phires_, etares_, d_, timeOut_:5] := 
TimeConstrained[
	Eigenvalues[mutX[iP, iQ, hmut, vmut, hres, vres, d], 1], (* Eigenvalues[matrix, 1] gives the leading eigenvalue *)
	timeOut, (* try to find eigenvalues for timeOut number of seconds *)
	Eigenvalues[N[mutX[iP, iQ, hmut, vmut, hres, vres, d]], 1] (* Do this if the previous calculation timed out *)
][[1]] (* This [[1]] is because eigenvalues are always returned as a list, even if there's only one *)



dGrowth[iP_, iQ_, hres_, vres_, d_, deltaH_, deltaV_, timeOut_:5] :=
Module[{dGdh, dGdv}, (* derivatives of growth in terms of horizontal transmission rate and vertical transmission rate *)
dGdh = (mutGrowth[iP, iQ, hres + deltaH, vres, hres, vres, d, timeOut] - 1)/deltaH; (* Subtract 1 from mutant growth rate b/c the resident growth rate in a population of residents is 1 (discrete time) *)
dGdv = (mutGrowth[iP, iQ, hres, vres + deltaV, hres, vres, d, timeOut] - 1)/deltaV;
Return[{dGdh, dGdv}]]


(* ::Section:: *)
(*Introduction*)


(* ::Text:: *)
(*Hi Grace! This script has the skeletons of functions to find ecological equilibria, invasion fitness, and the derivative of invasion fitness.*)
(*It's based on a model of organisms evolving to transmit a symbiont horizontally (contagiously) or vertically (from parent to offspring).*)
(*Like our model, the organisms can be in different "states," in this case locations P and Q (vs. heterotrophic and phototrophic in our model)*)


(* ::Text:: *)
(*State variables: (equivalent to our W, H, P)*)


(* ::Item:: *)
(*iP: the fraction of individuals in location P who are infected*)


(* ::Item:: *)
(*iQ: the fraction of individuals in location Q who are infected *)


(* ::Text:: *)
(*Traits under selection: (equivalent to our \[Eta], \[Phi])*)


(* ::Item:: *)
(*h: the horizontal transmission rate*)


(* ::Item:: *)
(*v: the vertical transmission rate*)


(* ::Text:: *)
(*Parameters: (equivalent to our d, c, etc.)*)


(* ::Item:: *)
(*d: the probability of dispersal between patches*)


(* ::Section:: *)
(*Function descriptions*)


(* ::Text:: *)
(*This is totally optional. It just makes ?Function say something about your function. *)
(*(You don't have to put them all together at the top like this. That's just a trick for making package-writing easier.)*)


(* Helper function *)
myJacobian::usage = "myJacobian finds the Jacobian of diPdt and diQdt, evaluated at ecological state (iP, iQ).";
negativeEigQ::usage = "negativeEigQ returns true if all eigenvalues of a matrix are negative.";

(* Ecological equilibria *)
diPdt::usage = "Rate of change in fraction of infected hosts in location P";
diQdt::usage = "Rate of change in fraction of infected hosts in location Q";
stableEquil::usage = "stableEquil finds the stable  equilibrium fraction of infected hosts (iP and iQ), in locations P and Q.";
stableEquilList::usage = "stableEquilList find stable ecological equilibria for a list of resident trait values.";

(* Mutant growth rate *)
mutX::usage = "mutX returns a matrix of the expected number of mutants produced by a mutant host when it is at low frequency";
mutGrowth::usage == "mutGrowth returns the growth rate of a mutant when rare";

(* Derivative of mutant growth rate *)
dGrowth::usage = "dGrowth estimates the partial derivative of the mutant growth rate, evaluated at the resident's transmission rates. Returns {\[PartialD]growth/\[PartialD]h, \[PartialD]growth/\[PartialD]v}.";
dGrowthList::usage = "dGrowthTable produces a list of derivatives of mutant growth rates when rare.";


(* ::Section:: *)
(*Helper functions*)


(* ::Text:: *)
(*myJacobian finds the Jacobian of diPdt and diQdt, the functions describing the ecological dynamics, for use in stability analysis.*)


myJacobian[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,h2,p2,c,gammaC,gammaP,sigmaP,hP,lH,dMax,sigmaD,sigmaC,lP]:=
 {{D[wDot[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,h2,p2], wDot], 
D[wDot[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,h2,p2], hDot],
D[wDot[w,hW,lW,a,pW,kappa,iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p,h, a2,h2,p2], pDot]}, 
{D[hDot[eta,c,gammaC,kappa,gammaP,pW,sigmaP,p,hP,iIn,phi,a,w,h,gammaA,sigmaA,lH,dMax,sigmaD], wDot], 
D[hDot[eta,c,gammaC,kappa,gammaP,pW,sigmaP,p,hP,iIn,phi,a,w,h,gammaA,sigmaA,lH,dMax,sigmaD], hDot],
D[hDot[eta,c,gammaC,kappa,gammaP,pW,sigmaP,p,hP,iIn,phi,a,w,h,gammaA,sigmaA,lH,dMax,sigmaD],pDot]},
{D[pDot[eta,gammaP,sigmaP,pW,p,kappa,hP,iIn,phi,a,w,h,gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],wDot],
D[pDot[eta,gammaP,sigmaP,pW,p,kappa,hP,iIn,phi,a,w,h,gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],hDot],
D[pDot[eta,gammaP,sigmaP,pW,p,kappa,hP,iIn,phi,a,w,h,gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],pDot]}}


(* ::Text:: *)
(*negativeEigQ takes a Jacobian and returns True if all eigenvalues are negative.*)
(*I know there's sometimes cycling in the mixotroph model, so potentially might also want to check for that.*)


negativeEigQ[jacobian_] :=
AllTrue[Eigenvalues[jacobian], Negative]


(* ::Section:: *)
(*Ecological Equilibria*)


(* ::Text:: *)
(*The functions diPdt and diQdt give the change in the fraction of infected hosts in each location per time step (time in units of host births).*)


diPdt[iP_, iQ_, h_, v_, d_] := 1 - (iP/2) (* function goes here *)


diQdt[iP_, iQ_, h_, v_, d_] := 1 - (iQ/3) (* function goes here *)


(* ::Text:: *)
(*stableEquil finds the stable  equilibrium fraction of infected hosts (iP and iQ), in patches P and Q.*)
(*"precision" is for the input "WorkingPrecision" in NSolve, which can be helpful if NSolve is being frustrating.*)
(*I'm using a module because I want to store the equilibria temporarily so I can do something with them (namely, discard the unstable equilibria).*)
(*The last bit of this function has an anonymous function. Anonymous functions use # as the variable (or #1, #2, etc for multiple variables) and end with &*)


stableEquil[h_, v_, d_, precision_: MachinePrecision] := 
Module[{equilibria},
  (* Find all equilibria *)
  equilibria = NSolve[{diPdt[iP, iQ, h, v, d] == 0, diQdt[iP, iQ, h, v, d] == 0}, {iP, iQ}, Reals, WorkingPrecision -> precision];
  
  (* Select the stable equilibria *)
  equilibria = Select[equilibria, negativeEigQ[myJacobian[iP/.#, iQ/.#, h, v, d]]&];
  
  (* Have equilibria be returned as numbers instead of iP \[Rule] number *)
  equilibria = Map[{iP/.#, iQ/.#}&, equilibria];
  
  Return[equilibria]
]


(* ::Text:: *)
(*stableEquilList takes as input a list of resident trait values at which to find stable equilibria. *)
(*Each entry of the list should be of the form {h, v}, where h is the horizontal transmission rate and v is the vertical transmission rate; e.g. {{0, 0}, {0, 0.1}, {0, 0.2}}*)
(**)
(*Output is returned in the form of a list with entries of the form {h, v,  equilibria}; e.g. *)
(*{{0, 0, {{0, 0}}}, (* only one stable equilibrium here *)*)
(*{0.5, 0.1, {{0.1,  0.5}, {0.24, 0.73}}}} (* Multiple stable equilibria here *)*)


stableEquilList[traitList_, d_, precision_:MachinePrecision] := 
Module[{equil}, 

(* Make a table of stable equilibria *)
equil=Table[equil[trait[[1]], trait[[2]], d, precision], {trait, traitList}];

(* Sort the equilibria for each pair of transmission rates *)
equil = Map[Sort, equil];

(* Add the trait values that produced the equilibrium/equilibria *)
equil = Table[{traitList[[i, 1]], traitList[[i, 2]], equil[[i]]}, {i, Length[traitList]}];

Return[equil]]


(* ::Section:: *)
(*Mutant growth rate (invasion fitness)*)


(* ::Text:: *)
(*To do adaptive dynamics, we need a formula for a mutant's growth rate when rare. The mutant's growth rate will probably depend on  its own traits, the resident's traits, and the environment.*)
(*Finding this growth rate falls in the gray area between coming up with a model and doing math on something you already have.*)


(* ::Text:: *)
(*In both the mixotroph model and the transmission model, one difficulty is that the mutant' s growth rate probably depends on what state it is in.*)
(*Is it a heterotroph or phototroph? Or, in the transmission model, is the mutant in location P or location Q?*)


(* ::Text:: *)
(*One way to find mutant fitness in this case to write down something like a Leslie matrix, that includes mutant birth, death, and transition between different states.*)
(*The matrix should look something like this:*)
(*X = ({*)
(* {H -> H, P -> H},*)
(* {H -> P, P -> P}*)
(*})*)
(*Where H -> H represents the expected number of heterotrophs in the next time step produced by 1 heterotroph in the current time step. That includes a heterotroph "producing" itself in the next time step by remaining alive and a heterotroph.*)
(*In the transmission model, I had to do some extra math at this point, because X depended on the population size, and I was interested in the infinite population case.*)
(**)
(*One you have X, its leading eigenvalue gives you the mutant growth rate.*)


mutX[iP_, iQ_, hmut_, vmut_, hres_, vres_, d_] :={{1, 0}, {0, 1}}


(* ::Text:: *)
(*Because finding eigenvalues exactly does not always go smoothly, I put a time constraint on it. If the exact eigenvalues can't be found quickly, it gives up and tries to find the eigenvalues of a numerical approximation.*)


mutGrowth[iP_, iQ_, hmut_, vmut_, hres_, vres_, d_, timeOut_:5] := 
TimeConstrained[
	Eigenvalues[mutX[iP, iQ, hmut, vmut, hres, vres, d], 1], (* Eigenvalues[matrix, 1] gives the leading eigenvalue *)
	timeOut, (* try to find eigenvalues for timeOut number of seconds *)
	Eigenvalues[N[mutX[iP, iQ, hmut, vmut, hres, vres, d]], 1] (* Do this if the previous calculation timed out *)
][[1]] (* This [[1]] is because eigenvalues are always returned as a list, even if there's only one *)


(* ::Section:: *)
(*Derivative of the mutant growth rate*)


(* ::Text:: *)
(*dGrowth estimates the derivative of a mutant's growth rate when rare. Because we are interested in two traits (h and v), it estimates the partial derivative for a change in each trait.*)
(*iP and iQ here are the ecological equilibrium fractions of infected hosts.*)
(*deltaH and deltaV are the amounts we'll add to the resident trait values when estimating the derivatives. I recommend make these really small. I used 1/1000 for the transmission evolution model.*)


dGrowth[iP_, iQ_, hres_, vres_, d_, deltaH_, deltaV_, timeOut_:5] :=
Module[{dGdh, dGdv}, (* derivatives of growth in terms of horizontal transmission rate and vertical transmission rate *)
dGdh = (mutGrowth[iP, iQ, hres + deltaH, vres, hres, vres, d, timeOut] - 1)/deltaH; (* Subtract 1 from mutant growth rate b/c the resident growth rate in a population of residents is 1 (discrete time) *)
dGdv = (mutGrowth[iP, iQ, hres, vres + deltaV, hres, vres, d, timeOut] - 1)/deltaV;
Return[{dGdh, dGdv}]]


(* ::Text:: *)
(*dGrowthList takes a list of resident traits and equilibria and finds the derivative of the mutant growth rate (derivative of invasion fitness) for each.*)
(*Input should be of the form produced by stableEquilList: a list with entries of the form {hres, vres, equilibria}. Equilibria is should be a list of stable equilibria, e.g. {{equil1P, equil1Q}, {equil2P,equil2Q}}.*)
(**)
(*Output is a list with entries of the form {hres, vres, derivatives}, where derivatives is a list of the form {{dGdh_equil1, dGdv_equil1}, {dGdh_equil2, dGdv_equil2}}.*)
(*You don't have to have 2 stable equilibria for this to work. Any number >= 1 will work.*)
(**)
(*dGrowthList takes into account that the traits are bounded between 0 and 1. It will subtract delta to get the mutant trait value when the resident trait value is 1.*)
(*delta_?Positive means the function only works for positive delta. You can get rid of that, it just means you'll need to deal with the case where delta is negative and the resident trait value is 0.*)


dGrowthList[equilList_, d_, delta_?Positive, timeOut_:5] :=
(* We'll go through all the entries of equilList in this table, calculating derivatives as we go *)
Table[ 
Module[{hres, vres, equil, deltaH, deltaV, derivatives}, 

(* Get the info from the entry of equilList we're currently on *)
hres=entry[[1]];
vres=entry[[2]]; 
equil = entry[[3]];

(* Set deltaH and deltaV to delta, unless that would make hmut or vmut > 1 *)
deltaH = If[hres + delta <= 1, delta, -delta];
deltaV = If[vres + delta <= 1, delta, -delta];

(* Use Map, which applies a function over a list, to get the derivatives of the growth rate for each equilibrium in equil *)
(* The #[[1]] and #[[2]] here refer to the first and second entries (iP and iQ) in each equilibrium listed in equil *)
derivatives = Map[dGrowth[#[[1]], #[[2]], hres, vres, d, deltaH, deltaV, timeOut]&, equil];

(* Join the resident trait values and derivatives into a list *)
{hres, vres, derivatives}], 

(* Last part of Table -- we're looping over the entry of equilList *)
{entry, equilList}]


(* ::Section:: *)
(*Plotting *)


(* ::Text:: *)
(*To plot stable equilibria, I use ListDensityPlot. The code below will make 2 separate plots:*)
(*One of the iP values for the first equilibrium for each entry in myStableEquil*)
(*One of the iQ values for the first equilibrium for each entry in myStableEquil*)


myStableEquil = {}; (* results of running stableEquilList *)

(* Get the values of iP and iQ at the first stable equilibrium listed for each trait combination *)
toPlotiP = Map[{#[[1]], #[[2]], #[[3, 1, 1]]}&, myStableEquil];
toPlotiQ = Map[{#[[1]], #[[2]], #[[3, 1, 2]]}&, myStableEquil];

ListDensityPlot[toPlotiP, PlotLabel -> "Ecological equilibrium values of iP", AxesLabel -> {"h", "v"}]
ListDensityPlot[toPlotiQ, PlotLabel -> "Ecological equilibrium values of iQ", AxesLabel -> {"h", "v"}]


(* Get any points with multiple stable equilibria *)
Select[myStableEquil, Length[#[[3]]] > 1&]

(* You can plot the second (or 3rd, 4th, etc.) equilibria of a point with this code. *)
myStableEquil2 = Select[myStableEquil, Length[#[[3]]] > 1&];
toPlotiP2 = Map[{#[[1]], #[[2]], #[[3, 2, 1]]}&, myStableEquil2];
toPlotiQ2 = Map[{#[[1]], #[[2]], #[[3, 2, 2]]}&, myStableEquil2];

ListDensityPlot[toPlotiP2, PlotLabel -> "2nd stable ecological equilibrium values of iP", AxesLabel -> {"h", "v"}]
ListDensityPlot[toPlotiQ2, PlotLabel -> "2nd stable ecological equilibrium values of iQ", AxesLabel -> {"h", "v"}]


(* ::Text:: *)
(*To plot evolutionary trajectories, I use ListStreamPlot*)


mydGrowth = {}; (* results of running dGrowthList *)

(* Get the values of dGrowthdh and dGrowthdV at the first stable equilibrium listed for each trait combination *)
toPlotGrowth = Map[{{#[[1]], #[[2]]}, {#[[3, 1, 1]], #[[3, 1, 2]]}}&, mydGrowth];

ListStreamPlot[toPlotGrowth, PlotLabel -> "Derivative of invasion fitness", AxesLabel -> {"h", "v"}]

(* Plot on top of the ecological equilibria for a fancy look *)
Show[ListDensityPlot[toPlotiP, PlotLabel -> "Ecol. equil. iP, arrows show derivative of invasion fitness", AxesLabel -> {"h", "v"}], 
ListStreamPlot[toPlotGrowth]]

