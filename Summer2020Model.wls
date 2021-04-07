#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Initial Set-up*)


(* ::Subsection:: *)
(*Equations- Trade-offs included for d,d2, a, aP, cP,pP,pP2*)


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


pW=3
gammaP=2.9/3
sigmaP=0
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
gammaC=1
eta=0
eta2=0
iIn=9


(* ::Subsection::Closed:: *)
(*Initial Values*)


h2Init=0
p2Init=0


initialValueFunction=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==5,h[0]==5,p[0]==5, h2[0]==0,p2[0]==0},{w,h,p,h2,p2},{t,0,3000}]


wInit=w[2000]/.initialValueFunction[[1]]


hInit=h[2000]/.initialValueFunction[[1]]


pInit=p[2000]/.initialValueFunction[[1]]


(* ::Section::Closed:: *)
(*Initialize: before each run*)


listOfwStarList = {}


listOfhStarList = {}


listOfpStarList = {}


listOfEtaList = {}


listOfPhiList = {}


wStarHatMinList={}


wStarHatMaxList={}


hStarHatMinList={}


hStarHatMaxList={}


pStarHatMinList={}


pStarHatMaxList={}


mutValList={}


mutTraitList={}


mutTypeList={}


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


plotListEta= {}


plotListPhi = {}


(* ::Section::Closed:: *)
(*Adaptive Dynamics Algorithm*)


(* ::Text:: *)
(*Outer loop: Number of simulations*)
(*Inner loop: Number of mutations (usually fixed at 500)*)
(**)
(*Method:*)
(*1. Choose mutation and assign traits to mutant*)
(*2. Run the beginning of a competition*)
(*3. If mutant is doing well, its phi and eta values become the resident trait*)
(*4. Calculate equilibrium values for W,H,P*)
(*5. Restart simulation*)
(**)


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
wStarHat = Abs[Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],
hStarHat = Abs[Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],
pStarHat = Abs[Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

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
w[0]==Abs[wInit],h[0]==Abs[hInit],p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]],

wStarHat=Abs[Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfwStarList[[j]],wStarHat],

hStarHat=Abs[Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfhStarList[[j]],hStarHat],

pStarHat=Abs[Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)],

AppendTo[listOfpStarList[[j]],pStarHat],

AppendTo[listOfEtaList[[j]],eta],

AppendTo[listOfPhiList[[j]],phi]},500]}, {j,1,3}]


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








(* ::Section::Closed:: *)
(*Stability Analysis*)


(* ::Subsection::Closed:: *)
(*Preliminaries- Creating grid of phi and eta test values*)


(* ::Text:: *)
(*Create array of phi and eta values*)


phiEtaGrid=Array[1,{11,11}]


Do[{
Do[{
phiEtaGrid[[k]][[i]]={(i-1)*0.1,(k-1)*0.1}},{k,1,11}]},{i,1,11}]


phiEtaGrid


phiEtaGrid[[1]][[3]][[1]]


(* ::Subsubsection::Closed:: *)
(*Resetting- Before each stability analysis*)


vector1List={}


pairOfVec={}


displacementVecList={}


etaListTest={}


phiListTest={}


wStarListTest={}


pStarListTest={}


hStarListTest={}


meanWList={}


ess1Frequency = 0


ess1PercentFrequencyList = {}


(* ::Subsection::Closed:: *)
(*Mean displacement: 25 random mutations*)


Do[
{Do[
{pStarListTest={},hStarListTest={},wStarListTest={},etaListTest={},phiListTest={},mutTraitList={},mutTypeList={},
phi=phiEtaGrid[[i]][[j]][[1]],eta=phiEtaGrid[[i]][[j]][[2]],phi2=phi,eta2=eta,
Do[
{Do[
{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000},Method->{"StiffnessSwitching"}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000},Method->{"StiffnessSwitching"}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==Abs[pInit]},{w,h,p},{t,0,8000},Method->{"StiffnessSwitching"}][[1]],
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
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{Abs[hStarHat/(hStarHat+pStarHat)],Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],

h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000},Method->{"StiffnessSwitching"}][[1]],

p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,h2[t],p2[t]],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
h2'[t]==h2Dot[eta2,c,gammaC,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,lH2,dMax,sigmaD],
p2'[t]==p2Dot[eta2,gammaP,sigmaP,pW,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],gammaA,sigmaA,c,gammaC,sigmaC,dMax,sigmaD,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000},Method->{"StiffnessSwitching"}][[1]],

If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],

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
(*Plotting*)


Save["/Users/Grace/Desktop/vectorListOnlyALongRun4.wdx",pairOfVec]


Export["/Users/Grace/Desktop/OnlyALongRunVectorPlot4.pdf",ListVectorDensityPlot[pairOfVec,PlotLabel->"Mutant Trajectories",FrameLabel->{"Phi","Eta","",""},
ColorFunction->ColorData["Pastel"],PlotLegends->Automatic,VectorScale-> {Automatic, Automatic,Automatic}]]
