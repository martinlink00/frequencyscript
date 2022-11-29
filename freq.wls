#!/usr/bin/env wolframscript
(* ::Package:: *)

(*

Notes for usage:

- Effective potential terms will be calculated up to Subscript[V, 2*highestcoeff] and the boolean incond can be switched on if initial conditions are to be calculated.
- Before running the script, the kernel needs to be reset in order to delete all global parameters from memory.
- Note, that results are dumpsaved into the subfolder "exports", in order to save RAM and in order to not lose results if Mathematica crashes which it tends to do.

*)




(**************************************************************************)




total=
Timing[


(*INIT*)

ClearSystemCache[];
Off[FrontEndObject::notavail];

highestcoeffresp = ReadList["params.txt"][[1]];

If[Not[DirectoryQ["exports"]],
  Print["Creating exports directory."];CreateDirectory["exports"];SetDirectory["exports"];CreateDirectory["logs"];,
  SetDirectory["exports"];];
		
(*MISCDEF*)

logfilename = "freq_"<>DateString["ISODateTime"]<>".log";

LogPrint[x_,type_] := Module[{},
	Print[x];
	outp = "["<>DateString["ISODateTime"]<>"]\t["<>type<>"]\t[Maximum RAM: "<>ToString[Round[N[MaxMemoryUsed[]/(10^9) *1.0],0.01]]<>" GB ("<>ToString[Round[N[MaxMemoryUsed[]/MemoryAvailable[] *100.0],0.01]]<>"%)]\t"<>x;
	WriteLine[FileNameJoin[{"logs",logfilename}],outp];];

myByteCount[symbolName_String]:=
	Replace[
		ToHeldExpression[symbolName],
		Hold[x__]:>If[MemberQ[Attributes[x],
		Protected|ReadProtected],
		Sequence@@{},
		{ByteCount[Through[{OwnValues,DownValues,UpValues,SubValues,DefaultValues,FormatValues,NValues}[Unevaluated@x,Sort->False]]],symbolName}]
	];
CheckRAM[] :=
	With[
		{listing=myByteCount/@Names[]},
		Labeled[Grid[Reverse@Take[Sort[listing],-100],
		Frame->True,Alignment->Left],
		Column[{Style["ByteCount for symbols without attributes Protected and ReadProtected in all contexts",16,FontFamily->"Times"],
		Style[Row@{"Total: ",Total[listing[[All,1]]]," bytes for ",Length[listing]," symbols"},Bold]},Center,1.5],Top]
	];
PD[]:=Print["---------------------------------------------------------------------"];

AppendTo[$Messages,OpenWrite[FileNameJoin[{"logs",logfilename}]]];

(*INIT PRINTOUT*)

PD[];

If[$VersionNumber!=12.1,LogPrint["Version error: Script was developed on and for version 12.1.","ERROR"];Quit[],
LogPrint["### KAPITZA FREQUENCY SCRIPT ###","RUNNING"];LogPrint["Acquiring parameters from param.txt","RUN:PRM"]];

PD[];

highestcoeff = If[IntegerQ[highestcoeffresp],highestcoeffresp,LogPrint["Wrong input in params.txt file!","ERROR"];Quit[]];

userinfo = "Script has been executed on "<>ToString[$OperatingSystem]<>" by the user " <>ToString[$UserName]<>".";
versioninfo = "Running on Mathematica "<>ToString[$VersionNumber]<>".";
processorinfo = ToString[$ProcessorCount]<>" processors of type "<>ToString[$ProcessorType]<>" available.";
memoryinfo = ToString[N[DecimalForm[MemoryAvailable[]/10^6]]]<>" MB of memory available for allocation.";
paraminfo = "Calculating up to non-trivial correction "<>ToString[highestcoeff]<>"."(*<>" with params (rho,om_0) = ("<>ToString[\[Rho]]<>","<>ToString[\[Omega]]<>")."*);
esttimeinfo = "Estimate resource usage from logfiles by running paramest.wls.";
exportinfo = "Results are dumpsaved into exports.";
loginfo = "Log is saved under exports/logs/"<>logfilename<>".";
LogPrint["### SYSTEM INFORMATION ###","RUNNING"];
LogPrint[userinfo,"RUN:SYS"];LogPrint[versioninfo,"RUN:SYS"];LogPrint[processorinfo,"RUN:SYS"];LogPrint[memoryinfo,"RUN:SYS"];PD[];
LogPrint["### RUN INFORMATION ###","RUNNING"];LogPrint[paraminfo,"RUN:PRM"];LogPrint[esttimeinfo,"INFORM"];LogPrint[loginfo,"INFORM"];LogPrint[exportinfo,"INFORM"];PD[];
LogPrint["### NOW RUNNING SCRIPT ###","RUNNING"];PD[];

(*SCRIPT*)

LogPrint["Initiating.","PROGRESS"];
\[Epsilon]max = 2 * highestcoeff + 2; (*Max collected \[Epsilon]*)
highestpowercorr=highestcoeff; (*Highest power correction*)
e[n_,t]:=Exp[I \[Omega] n t];
\[CapitalPhi][n_] := \[Epsilon]^(Abs[n])*Piecewise[{{Sum[\[Phi][n,m][t]\[Epsilon]^m,{m,0,\[Epsilon]max-n}],n>=0}},Sum[(-1)^(n+m)\[Phi][-n,m][t] \[Epsilon]^m,{m,0,\[Epsilon]max-n}]];
\[Delta]\[Theta]= Sum[If[n==0,0,e[n,t]\[CapitalPhi][n]],{n,-\[Epsilon]max,\[Epsilon]max}]; (*\[Delta]\[Theta] Ansatz*)
\[Theta]= \[CapitalPhi][0] +\[Delta]\[Theta];
oddphi0repl =Flatten[Table[D[\[Phi][0,2i+1][t],{t,k}]->0,{i,0,highestcoeff+2},{k,0,\[Epsilon]max}]];
newnot = Flatten[Table[D[\[Phi][0,2j][t],{t,k}]->D[\[CurlyPhi][2j][t],{t,k}],{j,0,highestcoeff+2},{k,0,\[Epsilon]max}],2];
enrepl = Table[Exp[I \[Omega] n t]->x^n,{n,-\[Epsilon]max,\[Epsilon]max}];
PD[];

algtime =
Timing[

\[Theta]sercoeff = Series[\[Theta]/.oddphi0repl/.newnot,{\[Epsilon],0,\[Epsilon]max}][[3]];
\[Omega]sercoeff = D[\[Theta]sercoeff,t]/.enrepl;
\[Theta]sercoeff = \[Theta]sercoeff/.enrepl;
dercoeff = Series[\[Epsilon]^2 D[\[CapitalPhi][0],{t,2}] + Sum[If[n==0,0,e[n,t](\[Epsilon]^2 D[\[CapitalPhi][n],{t,2}] + 2 I \[Omega] n \[Epsilon] D[\[CapitalPhi][n],{t,1}]-\[Omega]^2 n^2 \[CapitalPhi][n])],{n,-\[Epsilon]max,\[Epsilon]max}]/.oddphi0repl /.newnot/.enrepl,{\[Epsilon],0,\[Epsilon]max}][[3]];
f[0][t_] := 1;
f[n_]'[t_] := n I \[CurlyPhi][0]'[t]f[n][t];
additiontheorem = {f[n_][t]f[m_][t]:>f[n+m][t],(f[n_][t])^m_:>f[n*m][t]};
addi[exp_] := Expand[exp,f[n_][t_]]//.additiontheorem;
sinser = Series[Sin[\[CurlyPhi][0][t]+\[Delta]],{\[Delta],0,\[Epsilon]max-1}][[3]]/.{Sin[\[CurlyPhi][0][t]]->1/(2 I) (f[1][t]-f[-1][t]),Cos[\[CurlyPhi][0][t]]->1/2 (f[1][t]+f[-1][t])};
bellargument = Table[Table[Factorial[i]\[Theta]sercoeff[[i+1]],{i,1,n}],{n,1,\[Epsilon]max}];
sinexpfunc[epsord_] := If[epsord<0,0,If[epsord==0,1/(2 I) (f[1][t]-f[-1][t]),1/Factorial[epsord] Sum[Factorial[n]sinser[[n+1]]BellY[epsord,n,bellargument[[epsord-n+1]]],{n,1,epsord}]]];

DumpSave["ini.mx",{(*\[Rho],\[Omega],*)LogPrint,logfilename,\[Theta],\[Delta]\[Theta],\[CapitalPhi],\[Epsilon]max,highestcoeff,myByteCount,CheckRAM,PD,newnot,oddphi0repl,additiontheorem,addi}]; (*Dumping basic definitions for later usage.*)


(*Sin treatment*)

LogPrint["Acquiring sin(theta) expansion.","PRG:INIT"];
sinexplist = Table[sinexpfunc[ee],{ee,0,\[Epsilon]max-1}]; (*Need only up to \[Epsilon]max-1, since one \[Epsilon] already on lhs.*)
sinexp[\[Epsilon]ord_] := If[\[Epsilon]ord<0,0,sinexplist[[\[Epsilon]ord+1]]];
dgl[\[Epsilon]ord_,nord_] := (Coefficient[Expand[dercoeff[[\[Epsilon]ord]]],x,nord]==Coefficient[Expand[\[Omega]^2 (sinexp[\[Epsilon]ord-2]-\[Rho]  * 1/(2 I) * (x-1/x)sinexp[\[Epsilon]ord-1])],x,nord]);
Share[]; (*In hope of saving memory.*)

(*Declare \[Phi][n\[NotEqual]0,m,t] solutions*)

PD[];LogPrint["Acquiring algebraic solutions up to eps^"<>ToString[\[Epsilon]max]<>".","PROGRESS"];


For[i=1,i<\[Epsilon]max+1,i++,
	PD[];LogPrint["Finding relevant terms at order eps^"<>ToString[i]<>".","PRG:ALG"];
	ex[i] = Expand[(\[Omega]^2 (sinexp[i-2]-\[Rho]  * 1/(2 I) * (x-x^-1)sinexp[i-1]))//.additiontheorem,Trig->False];
	Monitor[
		For[j=0,j<i,j++,
			If[i==1&&j==0,Clear[\[Phi]],<<phi.mx];
			cef = Coefficient[ex[i],x,i-j];
			If[j==0,
				solex=-1/((i-j)^2 \[Omega]^2) (cef),
				If[j==1,
				solex=-1/((i-j)^2 \[Omega]^2) (cef - 2 I \[Omega] (i-j)D[\[Phi][i-j,j-1][t],t]),
				solex=-1/((i-j)^2 \[Omega]^2) (cef- 2 I \[Omega] (i-j)D[\[Phi][i-j,j-1][t],t]-D[\[Phi][i-j,j-2][t],{t,2}])
				];
			];
			LogPrint["Declaring algebraic solution of phi[+-"<>ToString[i-j]<>","<>ToString[j]<>"][t].","PRG:ALG"];
			\[Phi][i-j,j][t_]= addi[solex];
			DumpSave["phi.mx",\[Phi]];
			Clear[\[Phi],solex];
		];
	,ProgressIndicator[j,{0,i-1}]];
];

ClearAll[ex,cef,solex];
<<phi.mx;

][[1]];

PD[];LogPrint["All phi[n!=0,m,t] declared. It took "<>ToString[Round[algtime,0.1]]<>" s.","INFORM"];PD[];

(*Effective potential treatment*)

vefftime = Timing[

LogPrint["Acquiring acceleration expressions.","PROGRESS"];
acc0[ord_] :=  Coefficient[Expand[\[Omega]^2 (sinexp[ord-2]-\[Rho]  * 1/(2 I) * (x-1/x)sinexp[ord-1]),x||f[n_][t_]],x,0]//.additiontheorem;
acc0list = Table[acc0[2 ee+2],{ee,0,highestpowercorr}];

For[m=0,m<highestcoeff+1,m++,
	For[n=0,n<\[Epsilon]max-1,n++,
		\!\(\*SuperscriptBox[\(\[CurlyPhi][2  m]\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", "+", "n"}], ")"}],
Derivative],
MultilineFunction->None]\)[t]=D[a[2m][t],{t,n}];
	];
];

LogPrint["Declaring varphi_2m''[t] and derivatives.","PRG:ACC"];

Monitor[
	For[m=0,m<highestcoeff,m++,
		a[2m][t_]=acc0list[[m+1]];
	],
	ProgressIndicator[m,{0,highestcoeff-1}]
];

DumpSave["a.mx",a];
ClearAll[a];
<<a.mx;
acc0list = addi[acc0list];

LogPrint["Declaring varphi_2m'[t].","PRG:ACC"];

Monitor[
	For[m=0,m<highestcoeff+1,m++,
		\[CurlyPhi][2m][t] = 0;
		\[CurlyPhi][2m]'[t] = 0;
		V[2m][t] = 0;
		V[2m]'[t] = 0;
	],
	ProgressIndicator[m,{0,highestcoeff+1}]
];


PD[];
LogPrint["Expanding acceleration expressions.","PRG:POT"];
acc0list = addi[acc0list];
DumpSave["acclist.mx",acc0list];
ClearAll[acc0list];
<<acclist.mx;

LogPrint["Acquiring frequency expressions..","PRG:POT"];

For[n=0,n<highestcoeff+1,n++,
	\[Omega][2n] = Collect[addi[-acc0list[[n+1]]]/.{f[n_][_]:>I n}/.\[Omega]->1,\[Rho],Expand];
];
DumpSave["om.mx",\[Omega]];

][[1]];

PD[];LogPrint["Frequency acquired up to correction "<>ToString[highestpowercorr]<>". It took "<>ToString[Round[vefftime,0.1]]<>" s.","INFORM"];

][[1]];

(*CLEANUP*)

PD[];LogPrint["Frequency calculation completed. It took "<>ToString[Round[total,0.1]]<>" s.","INFORM"];
LogPrint["Access solutions using Results.nb.","INFORM"];
PD[];
LogPrint["### FREQUENCY CALCULATION DONE ###","RUNNING"];
PD[];
