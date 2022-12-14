(*

Note, that results are dumpsaved into the subfolder "exports", in order to save RAM and in order to not lose results if Mathematica crashes which it tends to do.

*)

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
PD[]:=Print["---------------------------------------------------------------------"];

SmarterExpand[exp_] := Module[{list,expanded,head},
	head = ToString[Head[exp]];
	If[head==ToString[Plus],
		list = List@@exp;
		expanded = Map[Expand[#]&,list];
		Plus@@expanded,
	Print["Head was "<>head];
	Expand[exp]
	]
];

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
paraminfo = "Calculating up to non-trivial correction "<>ToString[highestcoeff]<>".";
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
e[n_,t]:=Exp[I n t];
\[CapitalPhi][n_] := ep^(Abs[n])*Piecewise[{{Sum[\[Phi][n,m][t]ep^m,{m,0,\[Epsilon]max-n}],n>=0}},Sum[(-1)^(n+m)\[Phi][-n,m][t] ep^m,{m,0,\[Epsilon]max-n}]];
\[Delta]\[Theta]= Sum[If[n==0,0,e[n,t]\[CapitalPhi][n]],{n,-\[Epsilon]max,\[Epsilon]max}]; (*\[Delta]\[Theta] Ansatz*)
\[Theta]= \[CapitalPhi][0] +\[Delta]\[Theta];
oddphi0repl =Flatten[Table[D[\[Phi][0,2i+1][t],{t,k}]->0,{i,0,highestcoeff+2},{k,0,\[Epsilon]max}]];
newnot = Flatten[Table[D[\[Phi][0,2j][t],{t,k}]->D[\[CurlyPhi][2j][t],{t,k}],{j,0,highestcoeff+2},{k,0,\[Epsilon]max}],2];
enrepl = Table[Exp[I n t]->x^n,{n,-\[Epsilon]max,\[Epsilon]max}];
PD[];

algtime =
Timing[

additiontheorem = {f[n_][t_]f[m_][t_]:>f[n+m][t],(f[n_][t_])^m_:>f[n*m][t]};
\[Theta]sercoeff = Series[\[Theta]/.oddphi0repl/.newnot,{ep,0,\[Epsilon]max}][[3]];
\[Omega]sercoeff = D[\[Theta]sercoeff,t]/.enrepl;
\[Theta]sercoeff = \[Theta]sercoeff/.enrepl;
dercoeff = Series[ep^2 D[\[CapitalPhi][0],{t,2}] + Sum[If[n==0,0,e[n,t](ep^2 D[\[CapitalPhi][n],{t,2}] + 2 I n ep D[\[CapitalPhi][n],{t,1}]- n^2 \[CapitalPhi][n])],{n,-\[Epsilon]max,\[Epsilon]max}]/.oddphi0repl /.newnot/.enrepl,{ep,0,\[Epsilon]max}][[3]];
f[0][t_] := 1;
f[n_]'[t_] := n I \[CurlyPhi][0]'[t]f[n][t];
sinser = Series[Sin[\[CurlyPhi][0][t]+\[Delta]],{\[Delta],0,\[Epsilon]max-1}][[3]]/.{Sin[\[CurlyPhi][0][t]]->1/(2 I) f[1][t] - 1/(2 I) f[-1][t],Cos[\[CurlyPhi][0][t]]->1/2 f[1][t]+1/2 f[-1][t]};
bellargument = Table[Table[Factorial[i]\[Theta]sercoeff[[i+1]],{i,1,n}],{n,1,\[Epsilon]max}];
sinexpfunc[epsord_] := If[epsord<0,0,If[epsord==0,1/(2 I) f[1][t]-1/(2 I) f[-1][t],1/Factorial[epsord] Sum[Factorial[n]sinser[[n+1]]BellY[epsord,n,bellargument[[epsord-n+1]]],{n,1,epsord}]]];

DumpSave["ini.mx",{LogPrint,logfilename,\[Theta],\[Delta]\[Theta],\[CapitalPhi],\[Epsilon]max,highestcoeff,PD,newnot,oddphi0repl,SmarterExpand,additiontheorem}]; (*Dumping basic definitions for later usage.*)

(*Sin treatment*)

LogPrint["Acquiring sin(theta) expansion.","PRG:INIT"];
sinexplist = Table[sinexpfunc[ee],{ee,0,\[Epsilon]max-1}]; (*Need only up to \[Epsilon]max-1, since one \[Epsilon] already on lhs.*)
sinexp[\[Epsilon]ord_] := If[\[Epsilon]ord<0,0,sinexplist[[\[Epsilon]ord+1]]];

(*Declare \[Phi][n\[NotEqual]0,m,t] solutions*)

PD[];LogPrint["Acquiring algebraic solutions up to eps^"<>ToString[\[Epsilon]max]<>".","PROGRESS"];

LogPrint["Finding relevant coefficients.","PRG:ALG"];
For[i=1,i<\[Epsilon]max+1,i++,
	LogPrint["Finding relevant terms at order eps^"<>ToString[i]<>".","PRG:ALG"];
	exbef[i] = Expand[sinexp[i-2]- r * 1/(2 I) * (x-x^-1)sinexp[i-1],x];
	For[j=0,j<i,j++,
		cf = Coefficient[exbef[i],x,i-j];
		cefbef[i,j] = Expand[-1/((i-j)^2)*cf]; 
		DumpSave["cefbef.mx",cefbef];
		ClearAll[cefbef];
		<<cefbef.mx;
	];
	DumpSave["exbef.mx",exbef];
	ClearAll[exbef];
	<<exbef.mx;
];
ClearAll[exbef];
For[i=1,i<\[Epsilon]max+1,i++,
	PD[];LogPrint["Declaring algebraic solutions at order eps^"<>ToString[i]<>".","PRG:ALG"];
	For[j=0,j<i,j++,
		LogPrint["Declaring algebraic solution of phi["<>ToString[i-j]<>","<>ToString[j]<>"][t].","PRG:ALG"];
		If[i==1&&j==0,Clear[\[Phi]],<<phi.mx];
		cef = SmarterExpand[cefbef[i,j]]//.additiontheorem;
		DumpSave["cef.mx",cef];
		ClearAll[cef];
		<<cef.mx;
		If[j==0,
			\[Phi][i-j,j][t_] = cef,
			If[j==1,
			\[Phi][i-j,j][t_] = cef + Expand[2/(i-j) I D[\[Phi][i-j,j-1][t],t]],
			\[Phi][i-j,j][t_] = cef + Expand[2/(i-j) I D[\[Phi][i-j,j-1][t],t]] + Expand[1/(i-j)^2 D[\[Phi][i-j,j-2][t],{t,2}]]
			];
		];
		DumpSave["phi.mx",\[Phi]];
		Clear[\[Phi]];
	];
];
ClearAll[ex,cef];
<<phi.mx;
Share[];
][[1]];
PD[];LogPrint["All phi[n!=0,m,t] declared. This took "<>ToString[Round[algtime,0.1]]<>" s.","INFORM"];PD[];

(*Effective potential treatment*)

vefftime = Timing[

LogPrint["Acquiring acceleration expressions.","PROGRESS"];
LogPrint["Declaring varphi_2m''[t] and derivatives.","PRG:ACC"];
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

acc0[ord_] :=  SmarterExpand[Coefficient[Expand[sinexp[ord-2]- r * 1/(2 I) * (x-1/x)sinexp[ord-1],x],x,0]]//.additiontheorem;
acc0list = Table[acc0[2 ee+2],{ee,0,highestpowercorr}];

For[m=0,m<highestcoeff,m++,
	LogPrint["Inserting into phi[0,"<>ToString[2 m]<>"]''[t].","PRG:ACC"];
	a[2m][t_]=acc0list[[m+1]];
	DumpSave["a.mx",a];
	ClearAll[a];
	<<a.mx;
	acc0list = Map[SmarterExpand[#]&,acc0list]//.additiontheorem;
	DumpSave["acclist.mx",acc0list];
	ClearAll[acc0list];
	<<acclist.mx;
];

PD[];LogPrint["Setting varphis to zero.","PRG:ACC"];
For[m=0,m<highestcoeff+1,m++,
	\[CurlyPhi][2m][t] = 0;
	\[CurlyPhi][2m]'[t] = 0;
];

LogPrint["Expanding varphi_2m''[t].","PRG:ACC"];
acc0list = Map[SmarterExpand[#]&,acc0list]//.additiontheorem;
DumpSave["acclist.mx",{acc0list}];
ClearAll[acc0list];
<<acclist.mx;
][[1]];
PD[];LogPrint["varphi[2m]''[t] acquired up to correction "<>ToString[highestpowercorr]<>". This took "<>ToString[Round[vefftime,0.1]]<>" s.","INFORM"];

LogPrint["Acquiring frequency expressions.","PROGRESS"];
For[n=0,n<highestcoeff+1,n++,
	\[Omega][2n] = Collect[-acc0list[[n+1]]/.{f[n_][_]:>I n},r,Expand];
];
DumpSave["om.mx",\[Omega]];
][[1]];

(*CLEANUP*)

PD[];LogPrint["Script run completed. It took "<>ToString[Round[total,0.1]]<>" s.","INFORM"];
PD[];
LogPrint["### FREQUENCY CALCULATION DONE ###","RUNNING"];
PD[];
