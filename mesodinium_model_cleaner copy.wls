#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Preliminaries*)


(* ::Subsection:: *)
(*Equations*)


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


(* ::Subsection:: *)
(*Trade-Off Equations*)


(* ::Text:: *)
(*The trade-offs equations are as follows:*)
(*d = dMax (1 - sigmaD eta)*)
(*d2 = dMax (1 - sigmaD eta2)*)
(*aH = aMax[1 - sigmaA phi eta] Where a=aMax in the code*)
(*aH2=a2Max[1 - sigmaA phi2 eta2] Where a2=a2Max in the code*)
(*aP=gammaA aH *)
(*aP2 = gammaA aH2 *)
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
(**)
(*To turn on:*)
(*Attack: Set gammaA=0.7, sigmaA=1,*)
(*Conversion: Set gammaC=0.7*)


(* ::Subsection:: *)
(*Parameters*)


pW=3;
gammaP=2.9/3;
sigmaP=0;
kW=0.1;
kH=0.05;
kP=0.15;
kH2=0.05; (* Changing this to match kH2 *)
kP2=0.15;
hW=50;
hP=50;
hP2=50;
lW=0.5;
lH=0.1;
lP=0.4;
lH2=0.1;
lP2=0.4;
dMax=0.3;
sigmaD=0;
a=0.15;
a2=0.15;
gammaA=1;
sigmaA=1;
c=0.1;
gammaC=1;
sigmaC = 0; (* Added sigmaC because it looks like it's not defined here *)
iIn=25;


(* ::Subsection:: *)
(*Initial Values*)


(* ::Text:: *)
(*Here, we set the initial values that are dependent on a set of parameters. These values will be used later whenever we are computing equilibria. *)
(*We also find the first equilibria values of W,H, and P.*)


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


functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,gammaA,eta,phi,eta2,phi2,sigmaA,p[t],h[t], a2,0,0],
h'[t]==hDot[eta,c,gammaC,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],gammaP,pW,sigmaP,p[t],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,lH,dMax,sigmaD],
p'[t]==pDot[eta,gammaP,sigmaP,pW,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],gammaA,sigmaA,gammaC,sigmaC,c,dMax,sigmaD,lP],
w[0]==wInit,h[0]==hInit,p[0]==Abs[pInit]},{w,h,p},{t,0,8000}][[1]]


wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat = Abs[pStarHat]


hStarHat = Abs[hStarHat]


(* ::Section:: *)
(*Adaptive Dynamics*)


(* ::Text:: *)
(*Outer loop: Number of simulations*)
(*Inner loop: Number of mutations (usually ~600)*)
(**)
(*Method:*)
(*1. Choose mutation and assign traits to mutant*)
(*2. Run the beginning of a competition*)
(*3. If mutant is doing well, its phi and eta values become the resident trait (Save phi and eta values)*)
(*4. Calculate equilibrium values for W,H,P (Save W,H,P values) If the population crashes, exit the simulation.*)
(*5. Restart simulation*)


(* ::Subsection:: *)
(*Initialize, Before Each *)


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


(* ::Subsection:: *)
(*Version 11 (8/26/2021)*)


(* I've added these because I think they also need to be initialized for the code block below to run *)
listOfEventTimeList = {};
listOfEventTimeList2 = {2};


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

AppendTo[listOfPhiList[[j]],phi]},600]},

{j,1}]


(* ::Subsection:: *)
(*Plotting*)


(* ::Text:: *)
(*Plot the equilibrium values of W,H,P for each simulation step.*)


ListLinePlot[{listOfwStarList[[1]],listOfhStarList[[1]],listOfpStarList[[1]]},PlotLegends->{"W","H","P"}]


(* ::Text:: *)
(*Plot the resulting phi and eta for each simulation step.*)


ListLinePlot[{listOfEtaList[[1]],listOfPhiList[[1]]},PlotRange->{0,1},
PlotLegends->{"Eta","Phi"}]


(* ::Section:: *)
(*Time-Series*)


(* ::Text:: *)
(*If we are interested in looking at the population dynamics on an ecological time-scale (either resident-only, or resident vs mutant), we can generate a time-series. *)
(*Choose the desired phi and eta and plug in those values.*)


(* ::Subsection:: *)
(*Resident Only: Choose phi and eta*)


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


(* ::Subsection:: *)
(*Resident Only: Iterate Through Light (Bifurcation Diagram)*)


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


(* ::Subsection:: *)
(*Resident vs Mutant: Chose phi, eta, phi2, eta2*)


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


(* ::Subsection:: *)
(*Plotting*)


Plot[{wFuncEq,hFuncEq,pFuncEq},{t,0,8000},PlotRange->{0,15},
PlotLegends->{"W","H","P"},PlotStyle->{Darker[Green],Blue,Black}]


(* ::Section:: *)
(*Stream Plots and Heat Maps*)


(* ::Text:: *)
(*In general, to create a stream plot, we first must create an 11x11 grid of phi and eta values.*)
(**)
(*Then,*)
(*For each value in the Phi/Eta Grid, *)
(*	Run a competition between resident and mutant (with a fixed mutation size) in each of the cardinal directions. *)
(*	Then save the average displacement from the original coordinate.*)
(*	If we are interested in generating a heat map, we also save the average values of whatever we are interested in.*)
(*Now, we have a list of coordinates and their respective displacements.*)


(* ::Subsection:: *)
(*Preliminaries*)


(* ::Subsubsection:: *)
(*Create Phi and Eta Grid*)


(* ::Text:: *)
(*Create an array of phi and eta values that we will use to iterate through. We only need to do this once.*)


phiEtaGrid=Array[1,{11,11}]


Do[{
Do[{
phiEtaGrid[[k]][[i]]={(i-1)*0.1,(k-1)*0.1}},{k,1,11}]},{i,1,11}]


(* ::Subsection:: *)
(*Stream Plot with Resident Pop (W,H, or P) Heat Map*)


(* ::Text:: *)
(*Initialize lists needed to create the plot.*)


Do[{
vector1List={},
pairOfVec={},
displacementVecList={},
meanPFracList = {},
meanWList = {},
pairOfVecW = {},
pairOfVecP = {},
pairOfVecH = {},
meanHList = {},
meanPList = {},
},1]


(* ::Text:: *)
(*Create the plot. *)
(*For each value in the Phi/Eta Grid, *)
(*	Run a competition between resident and mutant (with a fixed mutation size) in each of the cardinal directions. *)
(*	Then save the average displacement from the original coordinate.*)
(*	Additionally, save the average W, H, P values.*)
(*Now, we have a list of coordinates and their respective displacements and pop values.*)


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
(*Stream Plot with W Heat Map, H+P Heat Map, and P/(P+H) Heat Map*)


(* ::Text:: *)
(*Initialize lists needed to create the plot.*)


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
pairOfVecMix = {},
meanHList = {},
meanPList = {},
},1]


(* ::Text:: *)
(*Create the plot. *)
(*For each value in the Phi/Eta Grid, *)
(*	Run a competition between resident and mutant (with a fixed mutation size) in each of the cardinal directions. *)
(*	Then save the average displacement from the original coordinate.*)
(*	Additionally, save the average W, H+P, and P/(H+P) values.*)
(*Now, we have a list of coordinates and their respective displacements and pop values.*)


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
(*Plotting*)


(* ::Text:: *)
(*Define a color function to color 0 values Black and to rescale the color scheme to appropriate values. If there are outliers, color them orange.*)


OutofBoundsColors[c_] := (If[c ==0, Black, If[c > 15, Orange,ColorData["Candy"][Rescale[c,{0,6}]]]])


(* ::Text:: *)
(*Create the stream plot, overlaid with the desired values (such as the W Population heat map)*)


streamPlot = ListVectorDensityPlot[pairOfVecW,PlotLabel->"Trajectories with Resident Prey",FrameLabel->{"Phi","Eta","",""},
ColorFunction-> OutofBoundsColors, 
PlotRange->Full,InterpolationOrder->1,
VectorStyle-> Black,
PlotLegends->Automatic,VectorScale-> {Automatic, Automatic,Automatic}]


(* ::Text:: *)
(*Example of marking ESS locations with yellow stars or a thick yellow line*)


starborderESS = ListPlot[{{0.43,1},{1,0}}, PlotStyle -> Black, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 38]]]


starESS = ListPlot[{{0.43,0.995},{1,0}}, PlotStyle -> Yellow, PlotMarkers -> Style["\[FivePointedStar]", Rule[FontSize, 30]]]


thickerlineEss = Plot[1,{t,0.1,0.8},PlotStyle -> {Directive[Thickness[0.02], Black], Black}]


lineEss = Plot[1,{t,0.1,0.8},PlotStyle -> {Directive[Thickness[0.01], Yellow], Yellow}]


(* ::Text:: *)
(*Combine the desired plots*)


Show[streamPlot,thickerlineEss,lineEss,starborderESS,starESS]]]
