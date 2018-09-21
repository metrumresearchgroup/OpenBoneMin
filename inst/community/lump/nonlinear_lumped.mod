$PROBLEM  Final lumped model
$INPUT  DV DOSE TIME ID AMT MDV EVID II ADDL KSSINI R00 KDEG
$DATA  nmdata.csv IGNORE=@
$SUBROUTINES  ADVAN13 TOL=9

$MODEL
COMP=(FOB)
COMP=(SOB)
COMP=(OC)
COMP=(ROB)
COMP=(RANKL)
COMP=(TGFB)
COMP=(COMPLEX)
COMP=(RANK)
COMP=(BMD)

COMP=(ABS,DEFDOSE)
COMP=(CENTRAL)
COMP=(PERI)
COMP=(RANKL)

$PK

IC1 = 0.0040
IC2 = 0.0010
IC3 = 0.0012
IC4 = 0.0010
IC5 = 1.4001
IC6 = 0.2281
IC7 = 0.00022283
IC8 = 30082

; 4 ROB
PIC0 = IC6                      ; 0.228142
RIN4 = THETA(1)*24*30/1000          ; *(IC1+IC2)/PIC0
K40 = THETA(2)*24*30
E0PICROB = 0.8838*PIC0
EMAXPICROB = 3.9745
GAMMA = 1.810
EC50PICROBPAREN = (EMAXPICROB*IC6 / (PIC0 - E0PICROB)) - IC6
EC50PICROB = EXP(LOG(EC50PICROBPAREN))
A_INITIAL(4) = IC4

; 1 fast OB
K41 = K40*(1-THETA(3))
K10 = THETA(4)*24*30       ; k17aD
K20 = THETA(5)*24*30/1000
A_INITIAL(1) = IC1

; 2 slow OB
K42 = K40*THETA(3)
A_INITIAL(2) = IC2

; 3 OC
RIN3 = THETA(6)*24*30/1000
K30 = THETA(7)*24*30
E0PICOC = 0.878215*PIC0
EMAXPICOC = 2
EC50PICOCPAREN = (EMAXPICOC*IC6/(PIC0 - E0PICOC)) - IC6
EC50PICOC = EXP(LOG(EC50PICOCPAREN))           

E0RANKL = 3.8034
EMAXL = 0.4698
PIL0 = 0.000022283
LSURVOCGAM = 3.0902
EC50SURVINPAR = (E0RANKL - EMAXL)*(PIL0**LSURVOCGAM/(E0RANKL - 1)) - PIL0**LSURVOCGAM
EC50SURV = EXP(LOG(EC50SURVINPAR)/LSURVOCGAM)  
A_INITIAL(3) = IC3

; 5 RANKL-RELATED LUMPED STATE
RIN5 = THETA(8)*24*30/1000
K50 = THETA(9)*24*30
K75 = THETA(10)*24*30
K15 = THETA(17)*24*30
K25 = K15
KSS = THETA(21)
A_INITIAL(5) = IC5

; 6 ACTIVE TGF_BETA
K60 = THETA(11)*24*30
K36 = K60*IC6/IC3
A_INITIAL(6) = IC6

; 7 RANK-RANKL COMPLEX
K57 = THETA(12)*24*30         ; as a function of RANK[0]?
K70 = K75
A_INITIAL(7) = IC7

; 8 RANK-RELATED LUMPED STATE
RIN8 = THETA(13)*24*30
K18 = THETA(14)*24*30
K28 = K18
K58 = THETA(15)*24*30
K80 = THETA(16)*24*30
K78 = K75
A_INITIAL(8) = IC8

; 9 BMD
KOUT_BMD = THETA(18)*24*30/1000
POWOB = THETA(19)
POWOC = THETA(20)
BMD0 = 100
A_INITIAL(9) = BMD0

; 10-13 PK
CL   = THETA(22)/1000*24*30                  ; L/month
V1   = THETA(23)/1000
Q    = THETA(24)/1000*24*30
V2   = THETA(25)/1000
TVKA = THETA(26)/24*24*30
KA   = TVKA*(66/71.5)**(-0.577) 
F10  = THETA(27)
KSYN = R00*KDEG
A_INITIAL(13) = R00

;-------------------------------------------------------------
$DES

KINT  = 0.00795*24*30
KOUTL = 0.00290*24*30

C = 0.5*((A(11)/V1 - A(13) - KSSINI) + SQRT((A(11)/V1 - A(13) - KSSINI)**2 + 4*KSSINI*A(11)/V1))
CLTOT = CL + KINT*V1*A(13)/(KSSINI+C)
MIC = (KINT-KDEG)*C/(KSSINI+C)

DADT(10) = -KA*A(10)
DADT(11) = KA*A(10) - CLTOT*C - Q*(C-A(12)/V2)
DADT(12) = Q*(C-A(12)/V2)
DADT(13) = KSYN - (KDEG+MIC)*A(13)

DRUG = (KINT-KOUTL)*C/(KSS+C)

; H+2016
PICROB = E0PICROB + EMAXPICROB*A(6)**GAMMA/(A(6)**GAMMA + EC50PICROB**GAMMA)

; H+2018D
PICOC = E0PICOC + EMAXPICOC*A(6)/(A(6) + EC50PICOC)

; H-2218D
PIL = A(7)/10
LSURVOC = E0RANKL - (E0RANKL - EMAXL)*(PIL**LSURVOCGAM/(PIL**LSURVOCGAM + EC50SURV**LSURVOCGAM))
   
DADT(1) =K41*A(4) - K10*A(1)                                         ; fast OB 1
DADT(2) =K42*A(4) - K20*A(2)                                         ; slow OB 2
DADT(3) =RIN3 - K30*PICOC*LSURVOC*A(3)                               ; OC      3
DADT(4) =RIN4*PICROB - K40*A(4)                                      ; ROB     4
DADT(5) =RIN5 + K15*A(1) + K25*A(2) + K75*A(7) - (K50+DRUG/3)*A(5)   ; RANKL-RELATED LUMPED STATE 5
DADT(6) =K36*A(3) - K60*A(6)                                         ; TGFBETA 6
DADT(7) =K57*A(5) - K70*A(7)                                         ; RANK-RANKL COMPLEX 7
DADT(8) =RIN8 + K18*A(1) + K28*A(2) + K58*A(5) + K78*A(7) - K80*A(8) ; RANK-RELATED LUMPED STATE 8

BSAP = A(1) + A(2)
SCTX = A(3)

RIN_BMD = KOUT_BMD*BMD0*(BSAP/(IC1+IC2))**POWOB;
K_BMD = KOUT_BMD*(SCTX/IC3)**POWOC;

DADT(9)=RIN_BMD - K_BMD*A(9)                                        ; BMD 9

;-------------------------------------------------------------
$ERROR

BSAP2 = A(1) + A(2)
SCTX2 = A(3)
ROB = A(4)
RANKL = A(5)
ATGF = A(6)
COMPLEX = A(7)
BMD = A(9)
PBMD = (BMD-BMD0)/BMD0 * 100

Y=PBMD+ETA(1)

$THETA 
.003 FIX          ; RIN4        1   (0 .60552)
.003 FIX          ; K40         2

.06 FIX           ; FRAC        3
.01 FIX           ; K10         4
.001 FIX          ; K20         5

.00298 FIX        ; RIN3        6
(0 .0111)         ; K30         7   .0292 FIX

(0 .16)           ; RIN5        8
.0011 FIX         ; K50         9
.1120 FIX         ; K75         10

.0298 FIX         ; k60         11

.000019 FIX       ; k57         12  as a function of RANK[0]?

75 FIX            ; RIN8        13
55.3 FIX          ; k18         14
160 FIX           ; k58         15
.97 FIX           ; k80         16

.234 FIX          ; K15 and K25 17

(0 .146)          ; kout_BMD    18
.0739 FIX         ; POWOB       19  (0 .0739)
.0779 FIX         ; POWOC       20  (0 .0779) 

(0 100)           ; KSS         21

3.06 FIX          ; CL 22
2490 FIX          ; V1 23
37.9 FIX          ; Q  24
1360 FIX          ; V2 25
0.212 FIX         ; KA 26
0.638 FIX         ; F1 27

$OMEGA
1

$ESTIMATION  NSIG=3 SIGL=9 MAXEVAL=9999  PRINT=5
$COV
$TABLE   DV DOSE TIME ID AMT MDV EVID II ADDL C EC50SURV PBMD BSAP BSAP2 SCTX SCTX2 RANKL COMPLEX ATGF ROB 
EC50PICOC EC50SURV EC50PICROB
NOPRINT ONEHEADER FILE=sdtab

;---------------------------------------------------
;---------------------------------------------------
; END OF THE CODE
;---------------------------------------------------
;---------------------------------------------------
