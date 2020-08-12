#!/usr/bin/env wolframscript
(* ::Package:: *)

wDot[w_,hW_,lW_,a_,pW_,kappa_,iIn_,aP_,p_,h_, a2_,h2_,aP2_,p2_]:= 
(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-a w h-aP w p-a2 w h2 - aP2 w p2 
hDot[eta_,c_,cP_,kappa_,pP_,p_,hP_,iIn_,phi_,a_,w_,h_,aP_,lH_,d_]:=
(1-eta) (pP p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi a w h + (1-phi)c a w h+ (1-eta) c aP w p -lH h + d p
pDot[eta_,pP_,p_,kappa_,hP_,iIn_,phi_,a_,w_,h_,aP_,cP_,d_,lP_] :=
eta (pP p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi a w h +eta aP cP p w - d p - lP p
h2Dot[eta2_,c_,cP_,pP2_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,aP2_,lH2_,d2_] :=
(1-eta2) (pP2 p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 a2 w h2 + (1-phi2) c a2 w h2 +(1-eta2) c aP2 w p2 - lH2 h2+d2 p2
p2Dot[eta2_,pP2_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,aP2_,cP_,d2_,lP2_] :=
eta2 (pP2 p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 a2 w h2+eta2 aP2 cP p2 w-d2 p2-lP2 p2
kappa[kW_, kH_,kP_,kH2_,kP2_,w_,h_,p_,h2_,p2_]:=
kW w+kH h+ kP p+kH2 h2+kP2 p2 





pW=3
pP=2.9
pP2=2.9
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
d=0.5
d2=0.5
a=0.15
a2=0.15
aP=0
aP2=0
phi=0
phi2=0
c=0.1
cP=0.1
eta=0
eta2=0
iIn=9


h2Init=0
p2Init=0


NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,aP,p[t],h[t],a2,h2[t],aP2,p2[t]],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
h2'[t]==h2Dot[eta2,c,cP,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,lH2,d2],
p2'[t]==p2Dot[eta2,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,cP,d2,lP2],
w[0]==5,h[0]==5,p[0]==5, h2[0]==0,p2[0]==0},{w,h,p,h2,p2},{t,0,3000}]


wInit=w[2000]/.%7952[[1]]


hInit=h[2000]/.%7952[[1]]


pInit=p[2000]/.%7989[[1]]


(* ::Section:: *)
(*Initialize*)


wStarHatList={}


wStarHatMinList={}


wStarHatMaxList={}


hStarHatList={}


hStarHatMinList={}


hStarHatMaxList={}


pStarHatList={}


pStarHatMinList={}


pStarHatMaxList={}


mutValList={}


mutTraitList={}


mutTypeList={}


phiMutationList={}


etaMutationList={}


phi=0


eta=0


phi2=phi


eta2=eta


functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,aP,p[t],h[t],a2,0,aP2,0],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]]


wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat=-pStarHat


(* ::Section:: *)
(*The Loop (Individual Components Below)*)


Do[{mutTrait=RandomChoice[{phiChoice,etaChoice}],
If[String[mutTrait]===String[phiChoice],{mutValueLB=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Min[mutValueLB,1]},{mutValueLB=Max[RandomVariate[NormalDistribution[eta,0.02]],0],mutValue=Min[mutValueLB,1]}],
If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue],
mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),Abs[pStarHat/(hStarHat+pStarHat)]}->{hChoice,pChoice}],
If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}],
h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,aP,p[t],h[t],a2,h2[t],aP2,p2[t]],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
h2'[t]==h2Dot[eta2,c,cP,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,lH2,d2],
p2'[t]==p2Dot[eta2,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,cP,d2,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],
p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,aP,p[t],h[t],a2,h2[t],aP2,p2[t]],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
h2'[t]==h2Dot[eta2,c,cP,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,lH2,d2],
p2'[t]==p2Dot[eta2,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,cP,d2,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]],
If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2},{phi2=phi,eta2=eta}],
functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,aP,p[t],h[t],a2,0,aP2,0],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]],
wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),
AppendTo[wStarHatList,wStarHat],
hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),
AppendTo[hStarHatList,hStarHat],
pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000),
AppendTo[pStarHatList,pStarHat],
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
AppendTo[etaMutationList,eta],
AppendTo[phiMutationList,phi]},10]


ListLinePlot[{phiMutationList,etaMutationList},PlotStyle->{Red,Purple},
PlotLegends->{"Phi","Eta"},PlotLabel->"Phi and Eta Mutations"]


ListLinePlot[{wStarHatMaxList,wStarHatMinList,wStarHatList},PlotStyle->{{Dashing[Medium],Darker[Green]},{Dashing[Tiny],Darker[Green]},Darker[Green]},
PlotLegends->{"W* Max","W* Min","W* Average"},PlotLabel->"Equilibrium Values of W"]


ListLinePlot[{hStarHatMaxList,hStarHatMinList,hStarHatList},PlotStyle->{{Dashing[Medium],Orange},{Dashing[Tiny],Orange},Orange},PlotRange->{0,2},
PlotLegends->{"H* Max","H* Min","H* Average"},PlotLabel->"Equilibrium Values of H"]


ListLinePlot[{pStarHatMaxList,pStarHatMinList,pStarHatList},PlotStyle->{{Dashing[Medium],Darker[Cyan]},{Dashing[Tiny],Darker[Cyan]},Darker[Cyan]},PlotRange->{0,1.5},
PlotLegends->{"P* Max","P* Min","P* Average"},PlotLabel->"Equilibrium Values of P"]


ListLinePlot[{wStarHatList,hStarHatList,pStarHatList},PlotLegends->{"W*","H*","P*"},PlotStyle->{Darker[Green],Orange,Darker[Cyan]},
PlotLabel->"Equilibrium Values of W,H,P"]


Plot[{w[t]/.functionEquilibrium,h[t]/.functionEquilibrium,p[t]/.functionEquilibrium},{t,0,8000}]


(* ::Section:: *)
(*Adaptive Dynamics Cycle*)


mutTrait=RandomChoice[{phiChoice,etaChoice}]


If[String[mutTrait]===String[phiChoice],mutValue=Max[RandomVariate[NormalDistribution[phi,0.02]],0],mutValue=Max[RandomVariate[NormalDistribution[eta,0.02]],0]]


If[String[mutTrait]===String[phiChoice],phi2=mutValue,eta2=mutValue]


mutType=RandomChoice[{hStarHat/(hStarHat+pStarHat),pStarHat/(hStarHat+pStarHat)}->{hChoice,pChoice}]


If[String[mutType]===String[hChoice], {hPop=0.99*(hStarHat+pStarHat)-pStarHat,pPop=pStarHat,h2Pop=0.01*(hStarHat+pStarHat),p2Pop=p2Init}, {hPop=hStarHat,pPop=0.99*(hStarHat+pStarHat)-hStarHat,h2Pop=h2Init,p2Pop=0.01*(pStarHat+hStarHat)}]


NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,aP,p[t],h[t],a2,h2[t],aP2,p2[t]],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
h2'[t]==h2Dot[eta2,c,cP,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,lH2,d2],
p2'[t]==p2Dot[eta2,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,cP,d2,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}]


h2Star=h2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,aP,p[t],h[t],a2,h2[t],aP2,p2[t]],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
h2'[t]==h2Dot[eta2,c,cP,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,lH2,d2],
p2'[t]==p2Dot[eta2,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,cP,d2,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]]


p2Star=p2[2000]/.NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],iIn,aP,p[t],h[t],a2,h2[t],aP2,p2[t]],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
h2'[t]==h2Dot[eta2,c,cP,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,lH2,d2],
p2'[t]==p2Dot[eta2,pP2,p2[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],h2[t],p2[t]],hP2,iIn,phi2,a2,w[t],h2[t],aP2,cP,d2,lP2],
w[0]==wStarHat,h[0]==hPop,p[0]==pPop, h2[0]==h2Pop,p2[0]==p2Pop},{w,h,p,h2,p2},{t,0,2000}][[1]]


h2Star> h2Pop(1.02)


p2Star> p2Pop*1.02


(* ::Section:: *)
(*Determine winner*)


If[h2Star > h2Pop(1.02)&&p2Star>p2Pop(1.02),{phi=phi2,eta=eta2, AppendTo[etaMutationList,eta],AppendTo[phiMutationList,phi]},{phi2=phi,eta2=eta}]


phiMutationList


etaMutationList


(* ::Section:: *)
(*Find Equilibrium Values*)


functionEquilibrium=NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,aP,p[t],h[t],a2,0,aP2,0],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,d],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],aP,cP,d,lP],
w[0]==wInit,h[0]==hInit,p[0]==pInit},{w,h,p},{t,0,8000}][[1]]


wStarHat=Integrate[w[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


hStarHat=Integrate[h[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat=Integrate[p[t]/.functionEquilibrium,{t,7000,8000}]/(8000-7000)


pStarHat=-pStarHat
