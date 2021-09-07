#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Initial Set-up*)


(* ::Subsection:: *)
(*Equations- Trade-offs included for d,d2, a, aP, cP,pP,pP2*)


(* ::Text:: *)
(*Here are the ODE equations that govern the system.*)


wDot[w_,hW_,lW_,a_,pW_,kappa_,iIn_,gammaA_,eta_,phi_,eta2_,phi2_,sigmaA_,p_,h_, a2_,h2_,p2_]:= 
(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-(a(1-sigmaA (phi eta))) w h-(gammaA a(1-sigmaA (phi eta))) w p-(a2(1-sigmaA (phi2 eta2))) w h2 - (gammaA a2(1-sigmaA (phi2 eta2))) w p2 
hDot[eta_,c_,gammaC_,kappa_,gammaP_,pW_,sigmaP_,p_,hP_,iIn_,phi_,a_,w_,h_,gammaA_,sigmaA_,lH_,dMax_,sigmaD_]:=
(1-eta) (c/(gammaC c))((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi (a(1-sigmaA (phi eta))) w h + (1-phi)c (a(1-sigmaA (phi eta))) w h+ (1-eta) c (gammaA a(1-sigmaA (phi eta))) w p -lH h + (dMax (1-sigmaD eta))p
pDot[eta_,gammaP_,sigmaP_,pW_,p_,kappa_,hP_,iIn_,phi_,a_,w_,h_,gammaA_,sigmaA_,gammaC_,sigmaC_,c_,dMax_,sigmaD_,lP_] :=
eta ((gammaP pW + sigmaP eta(pW-gammaP pW)) p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi (a(1-sigmaA (phi eta))) w h +eta (gammaA a(1-sigmaA (phi eta))) (gammaC c) p w - (dMax (1-sigmaD eta)) p - lP p
h2Dot[eta2_,c_,gammaC_,gammaP_,sigmaP_,pW_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,gammaA_,sigmaA_,lH2_,dMax_,sigmaD_] :=
(1-eta2) (c/(gammaC c))((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 (a2(1-sigmaA (phi2 eta2))) w h2 + (1-phi2) c (a2(1-sigmaA (phi2 eta2))) w h2 +(1-eta2) c (gammaA a2(1-sigmaA (phi2 eta2))) w p2 - lH2 h2+ (dMax(1-sigmaD eta2)) p2
p2Dot[eta2_,gammaP_,sigmaP_,pW_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,gammaA_,sigmaA_,c_,gammaC_,sigmaC_,dMax_,sigmaD_,lP2_] :=
eta2 ((gammaP pW + sigmaP eta2(pW-gammaP pW)) p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 (a2(1-sigmaA (phi2 eta2))) w h2+eta2 (gammaA a2(1-sigmaA (phi2 eta2))) (gammaC c) p2 w-(dMax(1-sigmaD eta2)) p2-lP2 p2
kappa[kW_, kH_,kP_,kH2_,kP2_,w_,h_,p_,h2_,p2_]:=
kW w+kH h+ kP p+kH2 h2+kP2 p2 


(* ::Subsection::Closed:: *)
(*Trade-Off Equations*)


(* ::Text:: *)
(*The trade-offs equations are as follows:*)
(*d = dMax (1 - sigmaD eta)*)
(*d2 = dMax (1 - sigmaD eta2)*)
(*aH = aMax[1 - sigmaA phi eta] **Edited** Where a=aMax in the code*)
(*aH2=a2Max[1 - sigmaA phi2 eta2] **Edited** Where a2=a2Max in the code*)
(*aP=gammaA aH **Edited***)
(*aP2 = gammaA aH2 **Edited***)
(*cP=gammaC c *)
(*pP = gammaP pW + sigmaP eta (pW - gammaP pW)*)
(*pP2 = gammaP pW + sigmaP eta2 (pW - gammaP pW)*)
(**)
(*Previously, we had d=d2=0.5, aP=aP2=0,cP=cP2=1, pP=pP2=2.9*)
(**)
(*To turn off the trade-offs,  set *)
(*sigmaA=sigmaD=sigmaC=sigmaP=0, *)
(*gammaP=2.9/3*)
(*gammaA=0*)
(*gammaC=1*)
(*dMax=0.5*)


(* ::Subsection::Closed:: *)
(*Parameters*)


(* ::Text:: *)
(*Set parameter values for the simulation. We also usually set phi=eta=phi2=eta2=0 to initialize.*)


pW=3
gammaP=2.9/3
sigmaP=0
kW=0.1
kH=0.05
kP=0.15
kH2=0.015
kP2=0.15
hW=50
hP=50
hP2=50
lW=0.5
lH=0.1
lP=0.4
lH2=0.1
lP2=0.4
dMax=0.3
sigmaD=0
a=0.15
a2=0.15
gammaA=0
sigmaA=0
phi=0
phi2=0
c=0.1
gammaC=0.7
eta=0
eta2=0
iIn=25


(* ::Subsection::Closed:: *)
(*Initial Values*)


(* ::Text:: *)
(*Here, we set the initial values that are dependent on a set of parameters. These values will be used later whenever we are computing equilibria.*)


h2Init=0
p2Init=0


initialValueFunction=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==5,h[0]==5,p[0]==5, h2[0]==0,p2[0]==0},{w,h,p,h2,p2},{t,0,3000}]


wInit=w[2000]/.initialValueFunction[[1]]


hInit=Abs[h[2000]/.initialValueFunction[[1]]]


pInit=Abs[p[2000]/.initialValueFunction[[1]]]


(* ::Section:: *)
(*Initialize: before each run*)


(* ::Text:: *)
(*Before an adaptive dynamics simulation, initialize the lists needed. Also, set the starting resident equilibria.*)


Do[{
listOfwStarList = {},
listOfhStarList = {},
listOfpStarList = {},
listOfEtaList = {},
listOfPhiList = {},
mutPhiList = {},
mutEtaList = {}
},1]


phi=0


eta=0


phi2=phi


eta2=eta


functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]]


wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat = Abs[pStarHat]


hStarHat = Abs[hStarHat]


(* ::Subsection::Closed:: *)
(*If saving start and end values of H,P,W,H2,P2 competition*)


(* ::Text:: *)
(*If we are interested in looking at the individual mutant-resident competitions during the adaptive dynamics simulation, we initialize those lists here.*)


Do[{
listOfMeanStartP2List = {},
listOfMeanStartH2List = {},
listOfMeanEndH2List = {},
listOfMeanEndP2List = {},
listOfMeanStartHList = {},
listOfMeanStartPList = {},
listOfMeanEndHList = {},
listOfMeanEndPList = {},
listOfMeanEndWList = {},
listOfMeanStartWList = {},
listOfEventTimeList={},
listOfEventTimeList2 = {}
},1]


(* ::Section:: *)
(*Adaptive Dynamics Algorithm*)


(* ::Text:: *)
(*Outer loop: Number of simulations*)
(*Inner loop: Number of mutations (usually fixed at 500)*)
(**)
(*Method:*)
(*1. Choose mutation and assign traits to mutant*)
(*2. Run the beginning of a competition*)
(*3. If mutant is doing well, its phi and eta values become the resident trait (Save phi and eta values)*)
(*4. Calculate equilibrium values for W,H,P (Save W,H,P values)*)
(*5. Restart simulation*)
(**)
(*Use Version 11 onwards*)


(* ::Subsection::Closed:: *)
(*Version 1*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
phi = 0,
eta = 0,
phi2 = 0,
eta2 = 0,
mutValList = {},
mutTraitList = {},
mutTypeList = {},
mutChoiceList = {},

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
AppendTo[mutChoiceList,mutTrait],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
AppendTo[mutValList,mutValue],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi]},500]}, {j,1,2}]


(* ::Text:: *)
(*Extra code that I took out: for calculating w/p/h min and max*)


wStarHatMin=FindMinValue[{w[t]/.functionEquilibrium,t>7000},t],
AppendTo[wStarHatMinList,wStarHatMin],
wStarHatMax=FindMaxValue[{w[t]/.functionEquilibrium,t>7000},t],
AppendTo[wStarHatMaxList,wStarHatMax],
hStarHatMin=FindMinValue[{h[t]/.functionEquilibrium,t>7000},t],
AppendTo[hStarHatMinList,hStarHatMin],
hStarHatMax=FindMaxValue[{h[t]/.functionEquilibrium,t>7000},t],
AppendTo[hStarHatMaxList,hStarHatMax],
pStarHatMin=FindMinValue[{p[t]/.functionEquilibrium,t>7000},t],
AppendTo[pStarHatMinList,pStarHatMin],
pStarHatMax=FindMaxValue[{p[t]/.functionEquilibrium,t>7000},t],
AppendTo[pStarHatMaxList,pStarHatMax],












(* ::Subsection::Closed:: *)
(*Version 2 (Single iteration)*)


Do[{
phi = 1,
eta = 0.8,
phi2 = 1,
eta2 = 0.8,
mutValList = {},
mutTraitList = {},
mutTypeList = {},
mutChoiceList = {},

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
AppendTo[mutChoiceList,mutTrait],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
AppendTo[mutValList,mutValue],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList,wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList,hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList,pStarHat],

AppendTo[listOfEtaList,eta],

AppendTo[listOfPhiList,phi]},500]},1]


(* ::Subsection::Closed:: *)
(*Version 3 (4/7/2021)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
phi = 0,
eta = 0,
phi2 = 0,
eta2 = 0,
mutValList = {},
mutTraitList = {},
mutTypeList = {},
mutChoiceList = {},
pStarHat = Abs[pStarHat],
hStarHat = Abs[hStarHat],
wStarHat = Abs[wStarHat],

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
AppendTo[mutChoiceList,mutTrait],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi]},500]}, {j,1,2}]


(* ::Subsection::Closed:: *)
(*Version 4 (4/23/2021)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
phi = 0,
eta = 0,
phi2 = 0,
eta2 = 0,
mutValList = {},
mutTraitList = {},
mutTypeList = {},
mutChoiceList = {},

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
AppendTo[mutChoiceList,mutTrait],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
AppendTo[mutValList,mutValue],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[(h2Star + p2Star) > (h2Pop+p2Pop)*(1.00),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi]},500]}, {j,3,5}]


(* ::Subsection::Closed:: *)
(*Version 5 (4/27/2021)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
phi = 0,
eta = 0,
phi2 = 0,
eta2 = 0,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[j]],phi2],
AppendTo[mutEtaList[[j]],eta2],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,8000}][[1]],

meanStartP2 = Integrate[p2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2[t]/.functionEquilibrium,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2[t]/.functionEquilibrium,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > (meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi]},600]}, {j,1,3}]


(* ::Subsection::Closed:: *)
(*Version 6 (5/7/2021)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
phi = 0,
eta = 0,
phi2 = 0,
eta2 = 0,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[j]],phi2],
AppendTo[mutEtaList[[j]],eta2],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,8000}][[1]],

meanStartP2 = Integrate[p2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2[t]/.functionEquilibrium,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2[t]/.functionEquilibrium,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi]},600], 

AppendTo[scatterPlotPhis[[5]], {16,phi}],
AppendTo[scatterPlotEtas[[5]], {16,eta}]},

{j,1,5}]


scatterPlotPhis = {}


scatterPlotEtas = {}


AppendTo[scatterPlotEtas,{}]


Save["/Users/Grace/Desktop/AOnCOn-etaScatter-5-10-2021.wdx",scatterPlotEtas]


scatterPlotEtas


RandomReal[1,{10}]


Flatten[scatterPlotEtas,1]


scatterPlotEtas


Count[scatterPlotEtas[[4]],{14,1}]


scatterPlotEtas[[1]]


Do[
{Do[
{AppendTo[scatterPlotPhis[[j]][[i]],Count[scatterPlotPhis[[j]],scatterPlotPhis[[j]][[i]]]]},
{i,1,5}]},
{j,1,5}]


scatterPlotEtas


BubbleChart[Flatten[scatterPlotEtas,1],ChartStyle->"Pastel"]


BubbleChart[Flatten[scatterPlotPhis,1],ChartStyle->"Pastel"]


scatterPlotEtas


dummy1 = {{{1,1},{1,2}},{{2,2},{2,3}},{{3,3},{3,4}}}


Flatten[dummy1,1]


ListPlot[Flatten[dummy1,1]]


(* ::Subsection::Closed:: *)
(*Version 7 (6/2/2021, With Equilibrium)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
AppendTo[listOfEventTimeList,{}],
AppendTo[listOfEventTimeList2,{}],
phi = 0,
eta = 0,
phi2 = phi,
eta2 = eta,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[j]],phi2],
AppendTo[mutEtaList[[j]],eta2],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

eventTime2 = 8000,
functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}],
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0]
},{w,h,p},{t,0,eventTime2}][[1]],

wFuncEq = Piecewise[{{w[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
hFuncEq = Piecewise[{{h[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
pFuncEq = Piecewise[{{p[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],

wStarHat=Abs[Integrate[wFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[hFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[pFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi],

AppendTo[listOfEventTimeList[[j]],eventTime],
AppendTo[listOfEventTimeList2[[j]],eventTime2]},400], 

AppendTo[scatterPlotPhis, phi],
AppendTo[scatterPlotEtas, eta]},

{j,1,2}]


(* ::Subsection::Closed:: *)
(*Version 11 (8/26/2021 With Equilibrium- Fixed 0 in denominator problem)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
AppendTo[listOfEventTimeList,{}],
AppendTo[listOfEventTimeList2,{}],
phi = 0,
eta = 0,
phi2 = phi,
eta2 = eta,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[j]],phi2],
AppendTo[mutEtaList[[j]],eta2],
If[hStarHat+pStarHat==0, {Break[],Print["Broke"]},
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}]],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

eventTime2 = 8000,
functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}],
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0]
},{w,h,p},{t,0,eventTime2}][[1]],

wFuncEq = Piecewise[{{w[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
hFuncEq = Piecewise[{{h[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
pFuncEq = Piecewise[{{p[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],

wStarHat=Abs[Integrate[wFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[hFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[pFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi]},600], 

AppendTo[barPlotPhis[[3]], phi],
AppendTo[barPlotEtas[[3]], eta]},

{j,1}]


phi


eta


Do[{wStarMax = First[Maximize[{wFuncEq,7500<t<8000},t]],
wStarMin = First[Minimize[{wFuncEq,7500<t<8000},t]],
AppendTo[listOfwMaxList,wStarMax],
AppendTo[listOfwMinList,wStarMin],

hStarMax = First[Maximize[{hFuncEq,7500<t<8000},t]],
hStarMin = First[Minimize[{hFuncEq,7500<t<8000},t]],
AppendTo[listOfhMaxList,hStarMax],
AppendTo[listOfhMinList,hStarMin],

pStarMax = First[Maximize[{pFuncEq,7500<t<8000},t]],
pStarMin = First[Minimize[{wFuncEq,7500<t<8000},t]],
AppendTo[listOfhMaxList,pStarMax],
AppendTo[listOfhMinList,pStarMin],

mixMax = pStarMax + hStarMax,
mixMin = pStarMin + hStarMin,
AppendTo[listOfmixMaxList,mixMax],
AppendTo[listOfmixMinList,mixMin]},1]


(* ::Subsection::Closed:: *)
(*Version 12 (8/31/2021 With Equilibrium- Iterating through light)*)


(* ::Subsubsection::Closed:: *)
(*Initials*)


listOfwMaxList = {}


listOfwMinList = {}


listOfhMaxList = {}


listOfhMinList = {}


listOfpMaxList = {}


listOfpMinList = {}


listOfmixMaxList = {}


listOfmixMinList = {}


barPlotEtas = {}


barPlotPhis = {}


(* ::Subsubsection:: *)
(*Algorithm*)


Do[{iIn=15+(q-1),
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
Do[
{
phi = 0,
eta = 0,
phi2 = phi,
eta2 = eta,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[q]],phi2],
AppendTo[mutEtaList[[q]],eta2],
If[hStarHat+pStarHat==0, {Break[],Print["Broke"]},
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}]],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

eventTime2 = 8000,
functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}],
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0]
},{w,h,p},{t,0,eventTime2}][[1]],

wFuncEq = Piecewise[{{w[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
hFuncEq = Piecewise[{{h[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
pFuncEq = Piecewise[{{p[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],

wStarHat=Abs[Integrate[wFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[q]],wStarHat],

hStarHat=Abs[Integrate[hFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[q]],hStarHat],

pStarHat=Abs[Integrate[pFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[q]],pStarHat],

AppendTo[listOfEtaList[[q]],eta],

AppendTo[listOfPhiList[[q]],phi]},600], 

wStarMax = First[Maximize[{wFuncEq,7500<t<8000},t]],
wStarMin = First[Minimize[{wFuncEq,7500<t<8000},t]],
AppendTo[listOfwMaxList,wStarMax],
AppendTo[listOfwMinList,wStarMin],

hStarMax = First[Maximize[{hFuncEq,7500<t<8000},t]],
hStarMin = First[Minimize[{hFuncEq,7500<t<8000},t]],
AppendTo[listOfhMaxList,hStarMax],
AppendTo[listOfhMinList,hStarMin],

pStarMax = First[Maximize[{pFuncEq,7500<t<8000},t]],
pStarMin = First[Minimize[{wFuncEq,7500<t<8000},t]],
AppendTo[listOfhMaxList,pStarMax],
AppendTo[listOfhMinList,pStarMin],

mixMax = pStarMax + hStarMax,
mixMin = pStarMin + hStarMin,
AppendTo[listOfmixMaxList,mixMax],
AppendTo[listOfmixMinList,mixMin],

AppendTo[barPlotPhis, phi],
AppendTo[barPlotEtas, eta]},

1]},{q,1,10}]


listOfmixMaxList


dummy = {1,2,3,4}


dummy = {dummy[[1]],dummy[[2]],dummy[[3]]}


listOfmixMinList = {listOfmixMinList[[1]],listOfmixMinList[[2]],listOfmixMinList[[3]],listOfmixMinList[[4]]}


listO


{dummy[[1]],dummy[[2]],dummy[[3]]}


(* ::Subsection::Closed:: *)
(*Version 9 (6/2/2021, With Equilibrium, Saving Start and End Vals of H,P,H2,P2)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
AppendTo[listOfEventTimeList,{}],
AppendTo[listOfEventTimeList2,{}],
AppendTo[listOfMeanStartH2List,{}],
AppendTo[listOfMeanStartP2List,{}],
AppendTo[listOfMeanEndH2List,{}],
AppendTo[listOfMeanEndP2List,{}],
AppendTo[listOfMeanStartHList,{}],
AppendTo[listOfMeanStartPList,{}],
AppendTo[listOfMeanEndHList,{}],
AppendTo[listOfMeanEndPList,{}],
AppendTo[listOfMeanStartWList,{}]
AppendTo[listOfMeanEndWList,{}],
phi = 0,
eta = 0,
phi2 = 0,
eta2 = 0,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[j]],phi2],
AppendTo[mutEtaList[[j]],eta2],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],


meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanStartW = Integrate[wFunc,{t,0,200}]/(200),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),
meanStartP = Integrate[pFunc,{t,0,200}]/(200-0),
meanStartH = Integrate[hFunc,{t,0,200}]/(200-0),
meanEndP = Integrate[pFunc,{t,7800,8000}]/(200),
meanEndH = Integrate[hFunc,{t,7800,8000}]/(200),
meanEndW = Integrate[wFunc,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

eventTime2 = 8000,
functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}],
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0]
},{w,h,p},{t,0,eventTime2}][[1]],

wFuncEq = Piecewise[{{w[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
hFuncEq = Piecewise[{{h[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],
pFuncEq = Piecewise[{{p[t]/.functionEquilibriumWHP,t<= eventTime},{0,t>eventTime}}],

wStarHat=Abs[Integrate[wFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[hFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[pFuncEq,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi],

AppendTo[listOfMeanStartH2List[[j]],meanStartH2],
AppendTo[listOfMeanStartP2List[[j]],meanStartP2],
AppendTo[listOfMeanEndH2List[[j]],meanEndH2],
AppendTo[listOfMeanEndP2List[[j]],meanEndP2],
AppendTo[listOfMeanStartHList[[j]],meanStartH],
AppendTo[listOfMeanStartPList[[j]],meanStartP],
AppendTo[listOfMeanEndHList[[j]],meanEndH],
AppendTo[listOfMeanEndPList[[j]],meanEndP],
AppendTo[listOfMeanStartWList[[j]],meanStartW],
AppendTo[listOfMeanEndWList[[j]],meanEndW],

AppendTo[listOfEventTimeList[[j]],eventTime],
AppendTo[listOfEventTimeList2[[j]],eventTime2]},400]},

{j,1,2}]


(* ::Subsection::Closed:: *)
(*Version 8 (6/2/2021, Without Equilibrium for Problem Diagnosis)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
AppendTo[listOfMeanStartH2List,{}],
AppendTo[listOfMeanStartP2List,{}],
AppendTo[listOfMeanEndH2List,{}],
AppendTo[listOfMeanEndP2List,{}],
AppendTo[listOfMeanStartHList,{}],
AppendTo[listOfMeanStartPList,{}],
AppendTo[listOfMeanEndHList,{}],
AppendTo[listOfMeanEndPList,{}],
AppendTo[listOfMeanStartWList,{}],
AppendTo[listOfMeanEndWList,{}],
AppendTo[listOfEventTimeList,{}],
phi = 0,
eta = 0,
phi2 = phi,
eta2 = eta,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[j]],phi2],
AppendTo[mutEtaList[[j]],eta2],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanStartW = Integrate[wFunc,{t,0,200}]/(200),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),
meanStartP = Integrate[pFunc,{t,0,200}]/(200-0),
meanStartH = Integrate[hFunc,{t,0,200}]/(200-0),
meanEndP = Integrate[pFunc,{t,7800,8000}]/(200),
meanEndH = Integrate[hFunc,{t,7800,8000}]/(200),
meanEndW = Integrate[wFunc,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],


AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi],
AppendTo[listOfMeanStartH2List[[j]],meanStartH2],
AppendTo[listOfMeanStartP2List[[j]],meanStartP2],
AppendTo[listOfMeanEndH2List[[j]],meanEndH2],
AppendTo[listOfMeanEndP2List[[j]],meanEndP2],
AppendTo[listOfMeanStartHList[[j]],meanStartH],
AppendTo[listOfMeanStartPList[[j]],meanStartP],
AppendTo[listOfMeanEndHList[[j]],meanEndH],
AppendTo[listOfMeanEndPList[[j]],meanEndP],
AppendTo[listOfMeanStartWList[[j]],meanStartW],
AppendTo[listOfMeanEndWList[[j]],meanEndW],
AppendTo[listOfEventTimeList[[j]],eventTime]
},700], 

AppendTo[scatterPlotPhis, phi],
AppendTo[scatterPlotEtas, eta]},

{j,1,2}]


eventTime


(* ::Subsection::Closed:: *)
(*Version 10 (7/27/2021 Smaller Mutation Size to Avoid Crash)*)


Do[
{AppendTo[listOfwStarList,{}],
AppendTo[listOfhStarList,{}],
AppendTo[listOfpStarList,{}],
AppendTo[listOfEtaList,{}],
AppendTo[listOfPhiList,{}],
AppendTo[mutPhiList,{}],
AppendTo[mutEtaList,{}],
AppendTo[listOfMeanStartH2List,{}],
AppendTo[listOfMeanStartP2List,{}],
AppendTo[listOfMeanEndH2List,{}],
AppendTo[listOfMeanEndP2List,{}],
AppendTo[listOfMeanStartHList,{}],
AppendTo[listOfMeanStartPList,{}],
AppendTo[listOfMeanEndHList,{}],
AppendTo[listOfMeanEndPList,{}],
AppendTo[listOfMeanStartWList,{}]
AppendTo[listOfMeanEndWList,{}],
AppendTo[listOfEventTimeList,{}],
phi = 0,
eta = 0,
phi2 = phi,
eta2 = eta,

Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.002]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.002]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
AppendTo[mutPhiList[[j]],phi2],
AppendTo[mutEtaList[[j]],eta2],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanStartW = Integrate[wFunc,{t,0,200}]/(200),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),
meanStartP = Integrate[pFunc,{t,0,200}]/(200-0),
meanStartH = Integrate[hFunc,{t,0,200}]/(200-0),
meanEndP = Integrate[pFunc,{t,7800,8000}]/(200),
meanEndH = Integrate[hFunc,{t,7800,8000}]/(200),
meanEndW = Integrate[wFunc,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],


AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi],
AppendTo[listOfMeanStartH2List[[j]],meanStartH2],
AppendTo[listOfMeanStartP2List[[j]],meanStartP2],
AppendTo[listOfMeanEndH2List[[j]],meanEndH2],
AppendTo[listOfMeanEndP2List[[j]],meanEndP2],
AppendTo[listOfMeanStartHList[[j]],meanStartH],
AppendTo[listOfMeanStartPList[[j]],meanStartP],
AppendTo[listOfMeanEndHList[[j]],meanEndH],
AppendTo[listOfMeanEndPList[[j]],meanEndP],
AppendTo[listOfMeanStartWList[[j]],meanStartW],
AppendTo[listOfMeanEndWList[[j]],meanEndW],
AppendTo[listOfEventTimeList[[j]],eventTime]
},700], 

AppendTo[scatterPlotPhis, phi],
AppendTo[scatterPlotEtas, eta]},

{j,1,2}]


(* ::Subsection::Closed:: *)
(*Time-Series Without Event Time*)


mutSuccess=0


Do[{
phi = 0.7,
eta = 0,
eta2 = eta,
mutValue = phi + 0.05,
phi2 = mutValue,
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], 
{hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, 
{hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],
eventTime=8000,

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop
},{w,h,p,h2,p2},{t,0,8000}][[1]],


meanStartP2 = Integrate[p2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2[t]/.functionEquilibrium,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2[t]/.functionEquilibrium,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02 * (meanStartH2 + meanStartP2),{mutSuccess = 1},{mutSuccess=0}]

}, 1]


mutSuccess


meanEndH2


(* ::Subsection::Closed:: *)
(*Time-Series (Edit 6/2/21)*)


hStarHat


mutSuccessList = {}


mutSuccess = 0


wStarHat


hStarHat


pStarHat


Do[{
phi = 1,
eta = 0,
eta2 = eta,
mutValue = phi - 0.01,
phi2 = mutValue,
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], 
{hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, 
{hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],
eventTime=8000,

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02 * (meanStartH2 + meanStartP2),{mutSuccess = 1},{mutSuccess=0}]

}, 1]


eventTime


mutSuccess


(* ::Subsection::Closed:: *)
(*Time-Series- ESS Competition*)


Do[{
phi = 1,
eta = 0,
phi2 = 0.5,
eta2 = 1,
hPop=0.5 * hStarHat,
pPop=0.5 * pStarHat,
h2Pop=0.5 * hStarHat,
p2Pop=0.5 * pStarHat, 

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,8000}][[1]],

meanStartP2 = Integrate[p2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2[t]/.functionEquilibrium,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2[t]/.functionEquilibrium,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02 * (meanStartH2 + meanStartP2),{mutSuccess = 1},{mutSuccess=0}]

}, 1]


(* ::Subsection::Closed:: *)
(*Time-Series- ESS Competition (Edit 6/4/21)*)


Do[{
phi = 0.66,
eta = 0,
phi2 = 0.7,
eta2 = 0,
hPop=0.5 * hStarHat,
pPop=0.5 * pStarHat,
h2Pop=0.5 * hStarHat,
p2Pop=0.5 * pStarHat, 

eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2Func,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02 * (meanStartH2 + meanStartP2),{mutSuccess = 1},{mutSuccess=0}]

}, 1]


phi = 0.7


eta = 0


functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}]
},{w,h,p},{t,0,eventTime2}]


Plot[{w[t]/.functionEquilibriumWHP,h[t]/.functionEquilibriumWHP,p[t]/.functionEquilibriumWHP},{t,0,8000}]


(* ::Text:: *)
(**)


mutSuccess


(* ::Subsection::Closed:: *)
(*Time-Series Plain (no mutant)*)


Do[{
phi = 0.5333860973900172,
eta = 0.6271899570747382,
eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}],
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0]
},{w,h,p},{t,0,eventTime2}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}]},1]


(* ::Subsection::Closed:: *)
(*Time-Series (no mutant) Iterating Through Light*)


Do[
{listOfwMaxList={},listOfwMinList={},
listOfpMaxList={},listOfpMinList={},
listOfhMaxList={},listOfhMinList={},
listOfmixMaxList={},listOfmixMinList={}}
,1]


Do[{iIn=15+(q-1),
Do[{
phi = 0.5,
eta = 0.5,
eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

wStarMax = First[Maximize[{wFunc,7500<t<8000},t]],
wStarMin = First[Minimize[{wFunc,7500<t<8000},t]],
AppendTo[listOfwMaxList,wStarMax],
AppendTo[listOfwMinList,wStarMin],
hStarMax = First[Maximize[{hFunc,7500<t<8000},t]],
hStarMin = First[Minimize[{hFunc,7500<t<8000},t]],
AppendTo[listOfhMaxList,hStarMax],
AppendTo[listOfhMinList,hStarMin],
pStarMax = First[Maximize[{pFunc,7500<t<8000},t]],
pStarMin = First[Minimize[{wFunc,7500<t<8000},t]],
AppendTo[listOfhMaxList,pStarMax],
AppendTo[listOfhMinList,pStarMin],
mixMax = pStarMax + hStarMax,
mixMin = pStarMin + hStarMin,
AppendTo[listOfmixMaxList,mixMax],
AppendTo[listOfmixMinList,mixMin]},1]},{q,1,11}]


(* ::Subsection::Closed:: *)
(*Time-Series Grid 6/29/21*)


phiEtaGrid


phiEtaGrid[[11]][[11]]


(* ::Subsubsection:: *)
(*Initialize*)


plotVals = {}


listOfhFuncMax = {}


listOfhFuncMin = {}


listOfhFuncMean = {}


listOfpFuncMax = {}


listOfpFuncMean = {}


listOfwFuncMax = {}


listOfwFuncMean = {}


phiEtaGridHMax = {}


phiEtaGridHMean ={}


phiEtaGridHMin = {}


phiEtaGridPMean = {}


phiEtaGridWMean = {}


(* ::Subsubsection:: *)
(*Algorithm*)


Do[{
Do[{
phi = phiEtaGrid[[i]][[j]][[1]],
eta = phiEtaGrid[[i]][[j]][[2]],
eventTime=8000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

hFuncMax = First[Maximize[{hFunc,7500<t<8000},t]],
AppendTo[listOfhFuncMax, hFuncMax],

hFuncMin = First[Minimize[{hFunc,7500<t<8000},t]],
AppendTo[listOfhFuncMin,hFuncMin],

hFuncMean= Integrate[hFunc,{t,7500,8000}]/500,
AppendTo[listOfhFuncMean,hFuncMean],

pFuncMax = First[Maximize[{pFunc,7500<t<8000},t]],
AppendTo[listOfpFuncMax,pFuncMax],

pFuncMean = Integrate[pFunc,{t,7500,8000}]/500,
AppendTo[listOfpFuncMean,pFuncMean],

wFuncMax = First[Maximize[{wFunc,7500<t<8000},t]],
AppendTo[listOfwFuncMax,wFuncMax],

wFuncMean = Integrate[wFunc,{t,7500,8000}]/500,
AppendTo[listOfwFuncMean,wFuncMean],

If[(hFuncMax-hFuncMin) >1,popStatus = 1, 
If[(hFuncMax + pFuncMax + wFuncMax) < 10^(-10),popStatus = 0, popStatus = 2]],


AppendTo[plotVals,{phiEtaGrid[[i]][[j]][[1]],phiEtaGrid[[i]][[j]][[2]],popStatus}],
AppendTo[phiEtaGridWMean,{phiEtaGrid[[i]][[j]][[1]],phiEtaGrid[[i]][[j]][[2]],wFuncMean}],
AppendTo[phiEtaGridPMean,{phiEtaGrid[[i]][[j]][[1]],phiEtaGrid[[i]][[j]][[2]],pFuncMean}],
AppendTo[phiEtaGridHMean,{phiEtaGrid[[i]][[j]][[1]],phiEtaGrid[[i]][[j]][[2]],hFuncMean}],
AppendTo[phiEtaGridHMax,{phiEtaGrid[[i]][[j]][[1]],phiEtaGrid[[i]][[j]][[2]],hFuncMax}],
AppendTo[phiEtaGridHMin,{phiEtaGrid[[i]][[j]][[1]],phiEtaGrid[[i]][[j]][[2]],hFuncMin}]},

{i,1,11}]},{j,1,11}]


(* ::Subsection::Closed:: *)
(*Type Frequencies (edit 8/30/21)*)


(* ::Text:: *)
(*Make room for plotting frequencies*)


Do[{plotVals = {},
AppendTo[plotVals,{}],
AppendTo[plotVals,{}]},1]


(* ::Text:: *)
(*Reset type counters*)


Do[
{deadCount = 0,
kpCount = 0,
okpfCount = 0,
okpvCount = 0,
skCount = 0}
,1]


(* ::Text:: *)
(*Count frequencies for length of scatterplotEtas&phis*)


Do[
{If[scatterPlotEtas[[i]]> 0.99 && scatterPlotPhis[[i]]>0.99,kpCount += 1,
If[scatterPlotEtas[[i]]<0.005, skCount +=1,
If[scatterPlotPhis[[i]]>0.47 && scatterPlotPhis[[i]]<0.48,okpfCount += 1,
okpvCount +=1]]]},
{i,1,5}
]


type2BCount


type2BCount = 1


type3Count = 19


(* ::Text:: *)
(*Append to Plotting values*)


Do[
{AppendTo[plotVals[[1]],deadCount],
AppendTo[plotVals[[1]],kpCount],
AppendTo[plotVals[[1]],okpfCount],
AppendTo[plotVals[[1]],okpvCount],
AppendTo[plotVals[[1]],skCount]
},1]


plotVals


(* ::Text:: *)
(*Save phi and eta vals*)


Export["/Users/Grace/Desktop/Iin16AOnCOnPoint3C-Etas-6-7-21.wdx",scatterPlotEtas]


Export["/Users/Grace/Desktop/Iin16AOnCOnPoint3C-Phis-6-7-21.wdx",scatterPlotPhis]


Export["/Users/Grace/Desktop/Iin16AOnCOnPoint3C-EtaEx-6-7-21.wdx",listOfEtaList[[1]]]


Export["/Users/Grace/Desktop/Iin16AOnCOnPoint3C-PhiEx-6-7-21.wdx",listOfPhiList[[1]]]


(* ::Text:: *)
(*Plot*)


BarChart[plotVals, ChartStyle->"CandyColors", ChartLayout-> "Stacked",
ChartLegends->{"Dead zone","KP","OKP (Fixed)","OKP (Var)","SK"}, 
Frame-> {True, True, False, False},
FrameTicks->{{{1,15},{2,20},{3,25}},Automatic},
PlotLabel->"No Trade-Offs"]


(* ::Subsection:: *)
(*Plotting*)


(* ::Subsubsection::Closed:: *)
(*Time Series*)


Plot[{wFuncEq,hFuncEq,pFuncEq},{t,0,8000},PlotRange->{0,15},
PlotLegends->{"W","H","P"},PlotStyle->{Darker[Green],Blue,Black}]


Plot[{wFunc,hFunc,pFunc},{t,0,8000},PlotRange->{0,20},
PlotLegends->{"W","H","P"},PlotStyle->{Darker[Green],Blue,Black}]


ListLinePlot[{listOfwMaxList,listOfwMinList,listOfmixMaxList,listOfmixMinList},PlotStyle->{{Darker[Green]},{Darker[Green],Dashed},{Blue},{Blue,Dashed}},
PlotRange-> {{0,15}},Ticks-> {{{1,15},{2,16},{3,17},{4,18},{5,19},{6,20},{7,21},{8,22},
{9,23},{10,24},{11,25}},Automatic},PlotLegends->{"Prey Max", "Prey Min", "Mixotroph Max", "Mixotroph Min"},
PlotLabel->"Attack & Conversion Efficiency Trade-Offs"]


Length[listOfmixMaxList]


listOfmixMaxList


meanStartH2


Export["/Users/Grace/Desktop/timeSeriesiIn10-h2p2-phi10eta00vsphi099eta00-7-14-21.pdf",
Plot[{h2Func,p2Func},{t,0,8000},
PlotLegends->{"H2","P2"},PlotStyle->{{Thick,Purple},{Thick,Cyan}}]]


Export["/Users/Grace/Desktop/timeSeriesiIn10-hph2p2-phi10eta00vsphi099eta00-7-14-21.pdf",
Plot[{hFunc,pFunc,h2Func,p2Func},{t,0,8000},PlotRange->{0,15},
PlotLegends->{"H","P","H2","P2"},
PlotStyle->{{Thick,Orange},{Thick,Blue},{Thick,Purple,Dashed},{Thick,Cyan}}]]


(* ::Subsubsection:: *)
(*General*)


OutofBoundsColors[c_] := (If[c ==0, Blue, If[c > 15, Orange,ColorData["Candy"][Rescale[c,{0,6}]]]])


ListDensityPlot[plotVals, InterpolationOrder->1,
ColorFunction->(ColorData["CandyColors"][Rescale[#,{0,2}]]&),
ColorFunctionScaling->False,
PlotLabel-> "Resident Population Status", FrameLabel->{"Phi","Eta","",""},FrameTicks->True,
Mesh-> 9,
MeshStyle->Black,PlotLegends->Automatic]


ListLinePlot[{listOfwMaxList,listOfwMinList,listOfhMaxList,listOfhMinList},PlotStyle->{{Green},{Green,Dashed},{Blue},{Blue,Dashed}},
PlotRange-> {0,10},Ticks-> {{{1,15},{2,16},{3,17},{4,18},{5,19},{6,20},{7,21},{8,22},
{9,23},{10,24},{11,25}},Automatic}]


Maximize[{hFunc,7500<t<7999},t]


Minimize[{hFunc,7500<t<7999},t]


ListLinePlot[{listOfwStarList[[1]],listOfhStarList[[1]],listOfpStarList[[1]]},PlotLegends->{"W","H","P"}]


Length[listOfEtaList[[1]]]


ListLinePlot[{listOfEtaList[[1]],listOfPhiList[[1]]},PlotRange->{0,1},
PlotLegends->{"Eta","Phi"},PlotLabel->"C On, iIn=25"]


(* ::Subsubsection::Closed:: *)
(*Bar Chart*)


Rasterize[Export["/Users/Grace/Desktop/COn-barchart-8-30-21.pdf",
BarChart[{{0,1,0,0,0},{0,0,0,1,0},{1,0,0,0,0}},
ChartLayout-> "Stacked", ChartStyle->"CandyColors", ChartLegends->{"Dead zone","KP","OKP (Fixed)","OKP (Var)","SK"},
Frame-> {True, True, False, False},
FrameTicks->{{{1,15},{2,20},{3,25}},Automatic},
PlotLabel->"Conversion Efficiency Trade-Off"
]]]


(* ::Subsubsection::Closed:: *)
(*Dot Diagrams of Start and End Pop Values*)


ListPlot[{eventTimeList},
PlotStyle->{PointSize[0.015]},PlotRange->{0,2.2},
PlotLegends->{"Start H2","End H2"}]


Export["/Users/Grace/Desktop/meanStartEndW-COnGamma5iIn16-regMut-7-28-21.pdf",
ListPlot[{listOfMeanStartWList[[1]],listOfMeanEndWList[[1]]},
PlotStyle->{PointSize[0.015]},PlotRange->{0,6},
PlotLegends->{"Start W","End W"}]]


listOfhStarList[[1]]


ListLinePlot[{listOfwStarList[[1]],listOfhStarList[[1]],listOfpStarList[[1]]},
PlotStyle->{PointSize[0.015]},
PlotLegends->{"W","H","P"}]


Export["/Users/Grace/Desktop/phiEtawEq-COnGamma3iIn10-7-14-21.pdf",
ListLinePlot[{listOfEtaList[[1]],listOfPhiList[[1]]},
PlotStyle->{PointSize[0.015]},
PlotLegends->{"Eta","Phi"}]]


Export["/Users/Grace/Desktop/hStarValswEq-COnGamma3iIn10-7-14-21.wdx",
listOfhStarList[[1]]]


listOfPhiList[[1]]


Export["/Users/Grace/Desktop/phiValswEq-COnGamma3iIn10-7-14-21.wdx",
listOfPhiList[[1]]]


Export["/Users/Grace/Desktop/meanStartWVals-COnGamma5iIn16-regMut-7-28-21.wdx",
listOfMeanStartWList[[1]]]


Export["/Users/Grace/Desktop/meanEndWVals-COnGamma5iIn16-regMut-7-28-21.wdx",
listOfMeanEndWList[[1]]]


(* ::Subsubsection::Closed:: *)
(*Heat Maps of Phi Eta Grid with Wmax,Hmax,Hmin,Pmax*)


Export["/Users/Grace/Desktop/resStatusVals-ACOnGamma7iIn24-8-2-21.wdx"
,plotVals]


Rasterize[Export["/Users/Grace/Desktop/resStatus-ACOnGamma7iIn24-8-2-21.pdf",
ListDensityPlot[plotVals,PlotLegends->Automatic,InterpolationOrder->1,ColorFunction->(ColorData["CandyColors"][Rescale[#,{0,3}]]&),
ColorFunctionScaling->False,Mesh-> 10]]]


ListDensityPlot[phiEtaGridHMax,PlotLegends->Automatic,InterpolationOrder->0,ColorFunction->(ColorData["CandyColors"][Rescale[#,{0,15}]]&),
ColorFunctionScaling->False,Mesh-> 10,PlotRange->All,
PlotLabel->"H max"]


ListDensityPlot[phiEtaGridHMin,PlotLegends->Automatic,InterpolationOrder->0,ColorFunction->(ColorData["CandyColors"][Rescale[#,{0,15}]]&),
ColorFunctionScaling->False,Mesh-> 10,PlotRange->All,
PlotLabel->"H min"]


Rasterize[Export["/Users/Grace/Desktop/wMeanHeatMap-ACOnGamma7iIn24-8-2-21.pdf",
ListDensityPlot[phiEtaGridWMean,PlotLegends->Automatic,InterpolationOrder->0,ColorFunction->(ColorData["CandyColors"][Rescale[#,{0,10}]]&),
ColorFunctionScaling->False,Mesh-> 9,PlotRange->All,
PlotLabel->"W Mean"]]]


Export["/Users/Grace/Desktop/wMeanHeatMapVals-ACOnGamma7iIn24-8-2-21.wdx",phiEtaGridWMean]


(* ::Section:: *)
(*Stream Plots and Heat Maps*)


(* ::Subsection:: *)
(*Initialize (Do this first)*)


(* ::Subsubsection::Closed:: *)
(*Preliminaries: Create Phi and Eta Grid*)


(* ::Text:: *)
(*Create an array of phi and eta values that we will use to iterate through. We only need to do this once.*)


phiEtaGrid=Array[1,{11,11}]


Do[{
Do[{
phiEtaGrid[[k]][[i]]={(i-1)*0.1,(k-1)*0.1}},{k,1,11}]},{i,1,11}]


(* ::Subsubsection:: *)
(*Initialize lists*)


(* ::Text:: *)
(*Initialize lists needed during the stream plot/heat map calculations.*)


Do[{
vector1List={},
pairOfVec={},
displacementVecList={},
meanPFracList = {},
meanMixotrophList = {},
meanWList = {},
pairOfVecPFrac = {},
pairOfVecMix = {},
pairOfVecW = {},
pairOfVecP = {},
pairOfVecH = {},
meanHList = {},
meanPList = {},
},1]


(* ::Subsection::Closed:: *)
(*Mean displacement: 25 random mutations*)


Do[
{Do[
{pStarListTest={},hStarListTest={},wStarListTest={},etaListTest={},phiListTest={},mutTraitList={},mutTypeList={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{mutTrait=RandomChoice[{phiChoice,etaChoice}],
pStarHat = Abs[pStarHat],
hStarHat = Abs[hStarHat],
wStarHat = Abs[wStarHat],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==Abs[p2Pop]},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==Abs[wInit],h[0]==Abs[hInit],p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]],
wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[wStarListTest,wStarHat],

hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[hStarListTest,hStarHat],

pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[pStarListTest,pStarHat],

AppendTo[etaListTest,eta],

AppendTo[phiListTest,phi]},25],

vector1 = {phiListTest[[25]],etaListTest[[25]]},

AppendTo[vector1List,vector1]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],displacementVec}]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*Mean displacement: 4 fixed mutations (7/19/21 Version), no equilibria calculated*)


Do[
{Do[
{mutTraitList={},mutTypeList={},hStarListTest={},pStarListTest={},wStarListTest={},phiListTest={},etaListTest={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=2000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,1800,2000}]/(200),
meanEndH2 = Integrate[h2Func,{t,1800,2000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],


AppendTo[etaListTest,eta],

AppendTo[phiListTest,phi],

vector1 = {phi,eta},

AppendTo[vector1List,vector1]},{q,1,4}]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],displacementVec}]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*Mean displacement: 4 fixed mutations (8/24/21 Version), no equilibria calculated, resident population heat map generated*)


Do[
{Do[
{mutTraitList={},mutTypeList={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=2000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,1800,2000}]/(200),
meanEndH2 = Integrate[h2Func,{t,1800,2000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

vector1 = {phi,eta},

AppendTo[vector1List,vector1],

phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
eventTime2=8000,
functionEquilibrium2 = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}]
},{w,h,p},{t,0,eventTime2}][[1]],

wFunc2 = Piecewise[{{w[t]/.functionEquilibrium2,t<= eventTime2},{0,t>eventTime2}}],
hFunc2 = Piecewise[{{h[t]/.functionEquilibrium2,t<= eventTime2},{0,t>eventTime2}}],
pFunc2 = Piecewise[{{p[t]/.functionEquilibrium2,t<= eventTime2},{0,t>eventTime2}}],

hFuncMean = Integrate[hFunc2,{t,7500,8000}]/500,
AppendTo[meanHList,hFuncMean],

pFuncMean = Integrate[pFunc2,{t,7500,8000}]/500,
AppendTo[meanPList,pFuncMean],

wFuncMean = Integrate[wFunc2,{t,7500,8000}]/500,
AppendTo[meanWList,wFuncMean]


},{q,1,4}]},

1],

meanVec = Mean[vector1List],

meanH = Mean[meanHList],

meanP = Mean[meanPList],

meanW = Mean[meanWList],

vector1List = {},

meanHList = {},
meanPList = {},
meanWList = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],displacementVec}],
AppendTo[pairOfVecH,{phiEtaGrid[[i]][[j]],{displacementVec,meanH}}],
AppendTo[pairOfVecP,{phiEtaGrid[[i]][[j]],{displacementVec,meanP}}],
AppendTo[pairOfVecW,{phiEtaGrid[[i]][[j]],{displacementVec,meanW}}]},

{j,1,11}]},
{i,1,11}]


(* ::Subsection:: *)
(*Mean displacement: 4 fixed mutations (8/26/21 Version), no equilibria calculated; W,H+P, P/(H+P) heat map generated*)


Do[
{Do[
{mutTraitList={},mutTypeList={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

eventTime=2000,
functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop,
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0], WhenEvent[h2[t]<0, h2[t]-> 0],
WhenEvent[p2[t]<0,p2[t]-> 0],
WhenEvent[w[t]+h[t]+p[t]+h2[t]+p2[t]<10^(-10),{eventTime=t, "StopIntegration"}]
},{w,h,p,h2,p2},{t,0,eventTime}][[1]],

wFunc = Piecewise[{{w[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
hFunc = Piecewise[{{h[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
pFunc = Piecewise[{{p[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
h2Func = Piecewise[{{h2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],
p2Func = Piecewise[{{p2[t]/.functionEquilibrium,t<= eventTime},{0,t>eventTime}}],

meanStartP2 = Integrate[p2Func,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2Func,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2Func,{t,1800,2000}]/(200),
meanEndH2 = Integrate[h2Func,{t,1800,2000}]/(200),

If[(meanEndH2+meanEndP2) > 1.02*(meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

vector1 = {phi,eta},

AppendTo[vector1List,vector1],

phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
eventTime2=8000,
functionEquilibrium2 = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit,
WhenEvent[w[t]+h[t]+p[t]<10^(-10),{eventTime2=t, "StopIntegration"}],
WhenEvent[w[t]<0, w[t]-> 0],WhenEvent[h[t]<0,h[t]-> 0], WhenEvent[p[t]<0,p[t]-> 0]
},{w,h,p},{t,0,eventTime2}][[1]],

wFunc2 = Piecewise[{{w[t]/.functionEquilibrium2,t<= eventTime2},{0,t>eventTime2}}],
hFunc2 = Piecewise[{{h[t]/.functionEquilibrium2,t<= eventTime2},{0,t>eventTime2}}],
pFunc2 = Piecewise[{{p[t]/.functionEquilibrium2,t<= eventTime2},{0,t>eventTime2}}],

hFuncMean = Integrate[hFunc2,{t,7500,8000}]/500,

pFuncMean = Integrate[pFunc2,{t,7500,8000}]/500,
AppendTo[meanMixotrophList,pFuncMean+hFuncMean],
If[pFuncMean+hFuncMean==0, {meanPFrac = 0}, {meanPFrac = Abs[pFuncMean] / (Abs[pFuncMean]+Abs[hFuncMean])}],
AppendTo[meanPFracList,meanPFrac],

wFuncMean = Integrate[wFunc2,{t,7500,8000}]/500,
AppendTo[meanWList,wFuncMean]


},{q,1,4}]},

1],

meanVec = Mean[vector1List],

meanMix = Mean[meanMixotrophList],

meanPFrac = Mean[meanPFracList],

meanW = Mean[meanWList],

vector1List = {},

meanMixotrophList = {},
meanPFracList = {},
meanWList = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],displacementVec}],
AppendTo[pairOfVecW,{phiEtaGrid[[i]][[j]],{displacementVec,meanW}}],
AppendTo[pairOfVecPFrac,{phiEtaGrid[[i]][[j]],{displacementVec,meanPFrac}}],
AppendTo[pairOfVecMix,{phiEtaGrid[[i]][[j]],{displacementVec,meanMix}}]},

{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*Mean displacement: 4 fixed mutations*)


Do[
{Do[
{mutTraitList={},mutTypeList={},hStarListTest={},pStarListTest={},wStarListTest={},phiListTest={},etaListTest={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]],
wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[wStarListTest,wStarHat],

hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[hStarListTest,hStarHat],

pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[pStarListTest,pStarHat],

AppendTo[etaListTest,eta],

AppendTo[phiListTest,phi],

vector1 = {phi,eta},

AppendTo[vector1List,vector1]},{q,1,4}]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],displacementVec}],

meanW = Mean[wStarListTest],

AppendTo[meanWList,meanW]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*4/28/21 Version: 4 fixed mutations*)


Do[
{Do[
{mutTraitList={},mutTypeList={},hStarListTest={},pStarListTest={},wStarListTest={},phiListTest={},etaListTest={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,8000}][[1]],

meanStartP2 = Integrate[p2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2[t]/.functionEquilibrium,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2[t]/.functionEquilibrium,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > (meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[wStarListTest,wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[hStarListTest,hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[pStarListTest,pStarHat],

AppendTo[etaListTest,eta],
AppendTo[phiListTest,phi],

vector1 = {phi,eta},

AppendTo[vector1List,vector1]},{q,1,4}]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],displacementVec}]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*Mean displacement: fixed mutations, prey heat-map*)


Do[
{Do[
{mutTraitList={},mutTypeList={},hStarListTest={},pStarListTest={},wStarListTest={},phiListTest={},etaListTest={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]],
wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[wStarListTest,wStarHat],

hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[hStarListTest,hStarHat],

pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[pStarListTest,pStarHat],

AppendTo[etaListTest,eta],

AppendTo[phiListTest,phi],

vector1 = {phi,eta},

AppendTo[vector1List,vector1]},{q,1,4}]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

meanW = Mean[wStarListTest],

AppendTo[meanWList,meanW],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],{displacementVec,meanW}}]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*4/29/2021 Version: 4 fixed mutations with prey equilibria*)


Do[
{Do[
{mutTraitList={},mutTypeList={},hStarListTest={},pStarListTest={},wStarListTest={},phiListTest={},etaListTest={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,8000}][[1]],

meanStartP2 = Integrate[p2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2[t]/.functionEquilibrium,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2[t]/.functionEquilibrium,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > (meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[wStarListTest,wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[hStarListTest,hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[pStarListTest,pStarHat],

AppendTo[etaListTest,eta],
AppendTo[phiListTest,phi],

vector1 = {phi,eta},

AppendTo[vector1List,vector1]},{q,1,4}]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

meanW = Mean[wStarListTest],

AppendTo[meanWList,meanW],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],{displacementVec,meanW}}]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*Mean displacement: fixed mutations, phototroph heat-map*)


Do[
{Do[
{mutTraitList={},mutTypeList={},hStarListTest={},pStarListTest={},wStarListTest={},phiListTest={},etaListTest={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]],
wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[wStarListTest,wStarHat],

hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[hStarListTest,hStarHat],

pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),

AppendTo[pStarListTest,pStarHat],

AppendTo[etaListTest,eta],

AppendTo[phiListTest,phi],

vector1 = {phi,eta},

AppendTo[vector1List,vector1]},{q,1,4}]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

meanP = Mean[pStarListTest],

AppendTo[meanWList,meanP],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],{displacementVec,meanP}}]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*4/29/2021 Version: 4 fixed mutations with phototroph equilibria*)


Do[
{Do[
{mutTraitList={},mutTypeList={},hStarListTest={},pStarListTest={},wStarListTest={},phiListTest={},etaListTest={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{If[q<=2,{mutTrait=phiChoice}, {mutTrait=etaChoice}],
If[q==1,{mutValueLB=phi+0.02, mutValue=Min[mutValueLB,1]}],
If[q==2,{mutValueLB=phi-0.02, mutValue=Max[mutValueLB,0]}],
If[q==3,{mutValueLB=eta+0.02, mutValue=Min[mutValueLB,1]}],
If[q==4,{mutValueLB=eta-0.02, mutValue=Max[mutValueLB,0]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

functionEquilibrium = NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,8000}][[1]],

meanStartP2 = Integrate[p2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanStartH2 = Integrate[h2[t]/.functionEquilibrium,{t,0,200}]/(200-0),
meanEndP2 = Integrate[p2[t]/.functionEquilibrium,{t,7800,8000}]/(200),
meanEndH2 = Integrate[h2[t]/.functionEquilibrium,{t,7800,8000}]/(200),

If[(meanEndH2+meanEndP2) > (meanStartH2 + meanStartP2),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

functionEquilibriumWHP=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[wStarListTest,wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[hStarListTest,hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibriumWHP,{t,7000,8000}]/(8000-7000)],
AppendTo[pStarListTest,pStarHat],

AppendTo[etaListTest,eta],
AppendTo[phiListTest,phi],

vector1 = {phi,eta},

AppendTo[vector1List,vector1]},{q,1,4}]},

1],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

meanP = Mean[pStarListTest],

AppendTo[meanPList,meanP],

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],{displacementVec,meanP}}]},
{j,1,11}]},
{i,1,11}]


(* ::Subsection::Closed:: *)
(*Mean displacement: full, random mutations, ESS heat-map*)


(* ::Text:: *)
(*Note: ESS 1 refers to when phi-> 1 and eta-> 0*)


Do[
{Do[
{pStarListTest={},hStarListTest={},wStarListTest={},etaListTest={},phiListTest={},mutTraitList={},mutTypeList={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{mutTrait=RandomChoice[{phiChoice,etaChoice}],
pStarHat = Abs[pStarHat],
hStarHat = Abs[hStarHat],
wStarHat = Abs[wStarHat],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==Abs[wInit],h[0]==Abs[hInit],p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

hStarHat=Abs[Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

pStarHat=Abs[Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[etaListTest,eta],

AppendTo[phiListTest,phi]},500],

vector1 = {phiListTest[[25]],etaListTest[[25]]},

AppendTo[vector1List,vector1],

If[etaListTest[[500]] < 0.6,{ess1Frequency += 1}],

etaListTest = {},

phiListTest = {}},

10],

meanVec = Mean[vector1List],

vector1List = {},

displacementVec = meanVec-phiEtaGrid[[i]][[j]],

ess1PercentFrequency = ess1Frequency / 10,

AppendTo[ess1PercentFrequencyList,ess1PercentFrequency],

ess1Frequency = 0,

AppendTo[pairOfVec,{phiEtaGrid[[i]][[j]],{displacementVec,ess1PercentFrequency}}]},
{j,1,11}]},
{i,1,11}]





(* ::Subsection::Closed:: *)
(*C=0.5 Eta Investigations*)


Do[
{mutChoiceList={},mutValList={},listOfEtaList={},listOfPhiList={},listOfhStarList={},listOfpStarList={},listOfwStarList={},
phi=phiEtaGrid[[11]][[11]][[1]],eta=phiEtaGrid[[i]][[1]][[2]],phi2=phi,eta2=eta,
Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
AppendTo[mutChoiceList,mutTrait],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
AppendTo[mutValList,mutValue],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList,wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList,hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList,pStarHat],

AppendTo[listOfEtaList,eta],

AppendTo[listOfPhiList,phi]},500],

AppendTo[plotEtaVals,{phiEtaGrid[[i]][[1]][[2]],listOfEtaList[[500]]}],
AppendTo[plotPhiVals,{phiEtaGrid[[i]][[1]][[2]],listOfPhiList[[500]]}],
AppendTo[plotHVals,{phiEtaGrid[[i]][[1]][[2]],listOfhStarList[[500]]}],
AppendTo[plotPVals,{phiEtaGrid[[i]][[1]][[2]],listOfpStarList[[500]]}],
AppendTo[plotWVals,{phiEtaGrid[[i]][[1]][[2]],listOfwStarList[[500]]}]},

{i,1,11}]


(* ::Subsection:: *)
(*Plotting*)


Save["/Users/Grace/Desktop/NoSigmaNoGammavectorListOnlyALongRunPartial.wdx",pairOfVec]


ListPlot[{plotEtaVals,plotPhiVals,plotWVals,plotHVals,plotPVals},PlotLegends->{"End Eta Value","End Phi Value","End W Value","End H Value","End P Value"},
PlotStyle->PointSize[Large],AxesLabel->{"Beginning Eta Value"}]


Export["/Users/Grace/Desktop/Point5CInvestigationLine-4-20-2021.pdf",ListLinePlot[{plotEtaVals,plotPhiVals,plotWVals,plotHVals,plotPVals},PlotLegends->{"End Eta Value","End Phi Value","End W Value","End H Value","End P Value"},
PlotMarkers->Automatic,AxesLabel->{"Beginning Eta Value"}]]


Save["/Users/Grace/Desktop/14IinAOffCOn-VectorList-5-5-2021.wdx",pairOfVec]


Export["/Users/Grace/Desktop/vectorTrajectoryVals-ACOnGamma7iIn24-8-2-21.wdx",pairOfVec]


ListDensityPlot[plotVals, InterpolationOrder->1,
ColorFunction->(ColorData["CandyColors"][Rescale[#,{0,2}]]&),
ColorFunctionScaling->False,
PlotLabel-> "Resident Population Status", FrameLabel->{"Phi","Eta","",""},FrameTicks->True,
Mesh-> 9,
MeshStyle->Black,PlotLegends->Automatic]


Rasterize[Export["/Users/Grace/Desktop/NOTOiIn10-VectorPlotW-8-24-21.pdf",



restrictedRangeColors[c_] := ColorData["Candy"][(4/20) + (16 c/20)]


OutofBoundsColors[c_] := If[c == 0, Blue, ColorData["Candy"][(4/20) + (16 c/20)]]


OutofBoundsColors[c_] := (If[c ==0, Blue, ColorData["Candy"][Rescale[c,{0,6}]]];
If[c > 10, Red, ColorData["Candy"][Rescale[c,{0,6}]]])


OutofBoundsColors[c_] := (If[c ==0, Black, If[c > 15, Orange,ColorData["CandyColors"][Rescale[c,{0,0.25}]]]])


streamPlot = ListVectorDensityPlot[pairOfVecPFrac,PlotLabel->"Trajectories with Resident W",FrameLabel->{"Phi","Eta","",""},
ColorFunction-> OutofBoundsColors, 
PlotRange->Full,InterpolationOrder->1,
VectorStyle-> Black,
PlotLegends->Automatic,VectorScale-> {Automatic, Automatic,Automatic}]


streamPlot = ListVectorDensityPlot[pairOfVecW,PlotLabel->"Trajectories with Resident W",FrameLabel->{"Phi","Eta","",""},
ColorFunction-> OutofBoundsColors,
ColorFunctionScaling->False,
VectorStyle-> Black,
PlotLegends->Automatic,VectorScale-> {Automatic, Automatic,Automatic}]


starborderESS = ListPlot[{{0.43,1},{1,0}}, PlotStyle -> Black, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 38]]];
starESS = ListPlot[{{0.43,0.995},{1,0}}, PlotStyle -> Yellow, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 30]]]


starESS = ListPlot[{{0.43,0.995},{1,0}}, PlotStyle -> Yellow, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 30]]]


thickerlineEss = Plot[1,{t,0.1,0.8},PlotStyle -> {Directive[Thickness[0.02], Black], Black}]


lineEss = Plot[1,{t,0.1,0.8},PlotStyle -> {Directive[Thickness[0.01], Yellow], Yellow}]


Show[thickerlineEss,lineEss,starborderESS,starESS]


Rasterize[Export["/Users/Grace/Desktop/ACOniIn25-VectorPlotW-blue-8-26-21.pdf",
Show[streamPlot,thickerlineEss,lineEss,starborderESS,starESS]]]


Clear[p]
