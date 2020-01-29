#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*Equations for prey, predator 1, predator 2*)
(*List of variables:*)
(*pW,pP,pP2: max photosynthetic rate*)
(*kW,kP,kH,kP2,kH2: absorption*)
(*hW,hP,hP2: half-saturation coefficient *)
(*lW,lP,lH,lP2,lH2: mortality rate*)
(*d1,d2: plastid decay rate*)
(*a1,a2,aP,aP2: attack rate*)
(*phi,phi2: fraction of predation events that lead to plastid acquisition*)
(*eta, eta2: fraction of division events leading to photosynthetic daughter cell*)
(*c,cP: conversion efficiency of prey to predator*)
(*w,p,h,p2,h2: population of prey(w), predator 1 (h,p), and predator 2 (h2,p2)*)
(**)


(* ::Input:: *)
(*wDot[w_,hW_,lW_,a_,pW_,kappa_,iIn_,aP_,p_,h_, a2_,h2_,aP2_,p2_]:= *)
(*(pW w)/kappa Log[(hW+iIn)/(hW+iIn Exp[-kappa])]-lW w-a w h-aP w p-a2 w h2 - aP2 w p2 *)
(*hDot[eta_,c_,cP_,kappa_,pP_,p_,hP_,iIn_,phi_,a_,w_,h_,aP_,lH_,d_]:=*)
(*(1-eta) c/cP (pP p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])] - phi a w h + (1-phi)c a w h+ (1-eta) c aP w p -lH h + d p*)
(*pDot[eta_,pP_,p_,kappa_,hP_,iIn_,phi_,a_,w_,h_,aP_,cP_,d_,lP_] :=*)
(*eta (pP p)/(kappa) Log[(hP+iIn)/(hP+iIn Exp[-kappa])]+phi a w h +eta aP cP p w - d p - lP p*)
(*h2Dot[eta2_,c_,cP_,pP2_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,aP2_,lH2_,d2_] :=*)
(*(1-eta2) c/cP (pP2 p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]-phi2 a2 w h2 + (1-phi2) c a2 w h2 +(1-eta2) c aP2 w p2 - lH2 h2+d2 p2*)
(*p2Dot[eta2_,pP2_,p2_,kappa_,hP2_,iIn_,phi2_,a2_,w_,h2_,aP2_,cP_,d2_,lP2_] :=*)
(*eta2 (pP2 p2)/(kappa) Log[(hP2+iIn)/(hP2+iIn Exp[-kappa])]+phi2 a2 w h2+eta2 aP2 cP p2 w-d2 p2-lP2 p2*)
(*kappa[kW_, kH_,kP_,kH2_,kP2_,w_,h_,p_,h2_,p2_]:=*)
(*kW w+kH h+ kP p+kH2 h2+kP2 p2 *)
(*?DSolve*)
(*DSolve[{y'[t]== y[t],x'[t]==1, y[0]==2,x[0]==5},{x,y},t]*)


(* ::InheritFromParent:: *)
(*iIn=9*)


(* ::Input:: *)
(*iIn ==9*)


(* ::Input:: *)
(*pW=3*)


(* ::Input:: *)
(*pP=2.9*)


(* ::Input:: *)
(*pP2=2.9*)


(* ::Input:: *)
(*kW=0.01*)


(* ::Input:: *)
(*kH=0.015*)


(* ::Input:: *)
(*kP=0.15*)


(* ::Input:: *)
(*kH2=0.015*)


(* ::Input:: *)
(*kP2=0.15*)


(* ::Input:: *)
(*hW=42*)


(* ::Input:: *)
(*hP=42*)


(* ::Input:: *)
(*hP2=42*)


(* ::Input:: *)
(*lW=0.5*)


(* ::Input:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::Input:: *)
(*lH=0.5*)


(* ::Input:: *)
(*lP=0.1*)


(* ::Input:: *)
(*lH2=0.5*)


(* ::Input:: *)
(*lP2=0.1*)


(* ::Input:: *)
(*d1=0.8*)


(* ::Input:: *)
(*d2=0.8*)


(* ::Input:: *)
(*a1=0.15*)


(* ::Input:: *)
(*a2=0.15*)


(* ::Input:: *)
(*aP=0*)


(* ::Input:: *)
(*aP2=0*)


(* ::Input:: *)
(*phi=1*)


(* ::Input:: *)
(*phi2=1*)


(* ::Input:: *)
(*c=0.1*)


(* ::Input:: *)
(*cP=0.1*)


(* ::Input:: *)
(*eta=0*)


(* ::Input:: *)
(*a = 0.15*)


(* ::Input:: *)
(*eta2=0*)


(* ::Input:: *)
(*?DSolve*)


(* ::Input:: *)
(**)
(*NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,aP,p[t],h[t],a2,0,aP2,0],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,0.8],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],aP,cP,0.8,lP],w[0]==1,h[0]==1,p[0]==1},{w,h,p},{t,0,100}]*)
(*?wDot[w_,hW_,lW_,a_,pW_,kappa_,iIn_,aP_,p_,h_,a2_,h2_,aP2_,p2_]*)


(* ::Input:: *)
(*Plot[ %51[[1,1]],{t,0,100}]*)


(* ::Input:: *)
(*%51*)


(* ::Input:: *)
(*w/.%51*)


(* ::Input:: *)
(*Plot[w[t]/.%51,{t,0,100}]*)


(* ::Input:: *)
(*Plot[{w[t]/.%51,p[t]/.%51,h[t]/.%51},{t,0,100},PlotStyle->{Blue, Black, Green}]*)


(* ::Input:: *)
(*w[20]*)


(* ::Input:: *)
(*w[20]/.%51*)


(* ::Input:: *)
(*Table[NDSolve[{w'[t]==wDot[w[t],hW,lW,a,pW,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],iIn,aP,p[t],h[t],a2,0,aP2,0],h'[t]==hDot[eta,c,cP,kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],pP,p[t],hP,iIn,phi,a,w[t],h[t],aP,lH,0.8],p'[t]==pDot[eta,pP,p[t],kappa[kW,kH,kP,kH2,kP2,w[t],h[t],p[t],0,0],hP,iIn,phi,a,w[t],h[t],aP,cP,0.8,lP],w[0]==1,h[0]==1,p[0]==1},{w,h,p},{t,0,100}],{pP,2.5,2.99,0.1}]*)


(* ::Input:: *)
(*?DSolve*)


(* ::Input:: *)
(*?NDSolve*)


(* ::InheritFromParent:: *)
(*{{y->Function[{t},2 E^t],x->Function[{t},5+t]}}*)


(* ::Input:: *)
(*wDot[1,40,0.5,0.15,0.3,1,10,0.05,1,1]*)


(* ::Input:: *)
(*?hDot*)


(* ::Input:: *)
(*?pDot*)


(* ::Input:: *)
(*?kappa*)


(* ::Input:: *)
(*%13*)


(* ::Input:: *)
(*%13[[1]]*)


(* ::Input:: *)
(*Plot[y[t]/.%13[[1,1]],{t,0,10}]*)


(* ::Input:: *)
(*ClearAll[wDot]*)
(**)
(**)
