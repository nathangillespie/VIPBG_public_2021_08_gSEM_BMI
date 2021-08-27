library(GenomicSEM)
library(lavaan)

# Read in S and V matrices
# Lavaan syntax
# Autoregression models
# 1-factor EFA
# Fit indices
# Model checks



# Read in S and V matrices

load("UKB23704_bmi_6ageTranche_20200618_men.RData")

S 					<- LDSCoutput$S
S 					<- as.matrix(Matrix::nearPD(S)$mat)
rownames(S) 		<- colnames(S)
eigen(S,T,T)$values
eigen(LDSCoutput$V,T,T)$values
V 					<- LDSCoutput$V
dim(LDSCoutput$V)
rownames(V) 		<- as.character(1:21)
colnames(V) 		<- rownames(V)
W 					<- chol2inv(chol(LDSCoutput$V))
#dataDWLS 			<- mxData(S, numObs = 2, means = NA, type = "acov", acov=diag(diag(V)), fullWeight=W)

S
cov2cor(S)

selVars 	<- c("BMI1","BMI2","BMI3","BMI4","BMI5","BMI6")

# mxOption(NULL,"mvnRelEps",0.0045)
 nv     	= length(selVars)
 nVariables	= nv
 nFactors	= nv
 psi_lab_a	= c("ia11","ia22","ia33","ia44","ia55","ia66")
 res_lab_e	= "res_e"
 psi_lab_vals=c(2.6,0.01,0.01,0.01,0.01,0.01)

# # Full
#  summary( auto_BMI_fit 		<- mxTryHard( auto_BMI,  extraTries=35, greenOK=FALSE,checkHess=FALSE,intervals=F) )
#  summary( auto_BMI_fit 		<- mxRun( auto_BMI_fit ) ) 
#  summary( auto_BMI_fit )  
# 
# # Sub1 
#  sub1 						<- mxModel( auto_BMI_fit,name="sub1")
#  sub1 						<- omxSetParameters( sub1, label="ia66",free=F,values=0)
#  summary( sub1_fit			<- mxRun( sub1,intervals=F))
#  summary( sub1_fit ) 
# 
# # Sub2
#  sub2 						<- mxModel( auto_BMI_fit,name="sub2")
#  sub2 						<- omxSetParameters( sub2, label="ia55",free=F,values=0)
#  summary( sub2_fit			<- mxRun( 		 sub2,intervals=F))
#  summary( sub2_fit ) 
# 
# # Sub3
#  sub3 						<- mxModel( auto_BMI_fit,name="sub3")
#  sub3 						<- omxSetParameters( sub3, label="ia44",free=F,values=0)
#  summary( sub3_fit			<- mxRun( 		 sub3,intervals=F))
#  summary( sub3_fit ) 
# 
# # Sub4
#  sub4 						<- mxModel( auto_BMI_fit,name="sub4")
#  sub4 						<- omxSetParameters( sub4, label="ia33",free=F,values=0)
#  summary( sub4_fit			<- mxRun( 		 sub4,intervals=F))
#  summary( sub4_fit ) 
# 
# # Sub5
#  sub5 						<- mxModel( auto_BMI_fit,name="sub5")
#  sub5 						<- omxSetParameters( sub5, label="ia22",free=F,values=0)
#  summary( sub5_fit			<- mxRun( 		 sub5,intervals=F))
#  summary( sub5_fit ) 
# 
# # Sub6
#  sub6 						<- mxModel( auto_BMI_fit,name="sub6")
#  sub6 						<- omxSetParameters( sub6, label=c("ia22","ia33","ia44","ia55","ia66"),free=F,values=0)
#  summary( sub6_fit			<- mxRun( 		 sub6,intervals=F))
#  summary( sub6_fit ) 
# 
#  subs <- c(sub1_fit,sub2_fit,sub3_fit,sub4_fit,sub5_fit,sub6_fit)
#  #mxCompare(auto_BMI_fit,subs)
 
 
# Autoregression models 
# Lavaan syntax:
auto_sntx <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ ia22*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ ia33*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ ia44*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ ia55*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ ia66*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
ia11 > 0
ia22 > 0
ia33 > 0
ia44 > 0
ia55 > 0
ia66 > 0
res > 0
'
autolav <- lavaan(
	model=auto_sntx,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav)
autolav@optim$fx

#Drop ia66:
auto_sntx_sub1 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ ia22*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ ia33*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ ia44*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ ia55*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ 0*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
ia11 > 0
ia22 > 0
ia33 > 0
ia44 > 0
ia55 > 0
res > 0
'
autolav_sub1 <- lavaan(
	model=auto_sntx_sub1,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub1)
autolav_sub1@optim$fx

#drop ia55:
auto_sntx_sub2 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ ia22*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ ia33*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ ia44*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ 0*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ ia66*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
ia11 > 0
ia22 > 0
ia33 > 0
ia44 > 0
ia66 > 0
res > 0
'
autolav_sub2 <- lavaan(
	model=auto_sntx_sub2,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub2)
autolav_sub2@optim$fx

#Drop ia44:
auto_sntx_sub3 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ ia22*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ ia33*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ 0*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ ia55*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ ia66*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
ia11 > 0
ia22 > 0
ia33 > 0
ia55 > 0
ia66 > 0
res > 0
'
autolav_sub3 <- lavaan(
	model=auto_sntx_sub3,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub3)
autolav_sub3@optim$fx

#Drop ia33:
auto_sntx_sub4 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ ia22*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ 0*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ ia44*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ ia55*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ ia66*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
ia11 > 0
ia22 > 0
ia44 > 0
ia55 > 0
ia66 > 0
res > 0
'
autolav_sub4 <- lavaan(
	model=auto_sntx_sub4,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub4)
autolav_sub4@optim$fx

#Drop ia22:
auto_sntx_sub5 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ 0*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ ia33*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ ia44*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ ia55*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ ia66*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
ia11 > 0
ia33 > 0
ia44 > 0
ia55 > 0
ia66 > 0
res > 0
'
autolav_sub5 <- lavaan(
	model=auto_sntx_sub5,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub5)
autolav_sub5@optim$fx

#Drop ia22 thru ia66:
auto_sntx_sub6 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ 0*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ 0*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ 0*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ 0*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ 0*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
ia11 > 0
res > 0
'
autolav_sub6 <- lavaan(
	model=auto_sntx_sub6,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub6)
autolav_sub6@optim$fx

#Model comparison:
#lavTestLRT(autolav,autolav_sub1,autolav_sub2,autolav_sub3,autolav_sub4,autolav_sub5,autolav_sub6)
#^^^Hmm, lavTestLRT() doesn't work the way I thought it might.
 
 
 
 
# 1-factor EFA

#lavaan syntax:
sntx_cfa1 <- '
F =~ 1*BMI1 + BMI2 + BMI3 + BMI4 + BMI5 + BMI6
F ~~ F
BMI1 ~~ BMI1
BMI2 ~~ BMI2
BMI3 ~~ BMI3
BMI4 ~~ BMI4
BMI5 ~~ BMI5
BMI6 ~~ BMI6
'  
lav_cfa1 <- lavaan(
	model=sntx_cfa1,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(lav_cfa1)
lav_cfa1@optim$fx



# Fit indices

#Pseudo-AICs:
autolav@Fit@test$scaled.shifted$stat + 2*autolav@Fit@npar
autolav_sub1@Fit@test$scaled.shifted$stat + 2*autolav_sub1@Fit@npar
autolav_sub2@Fit@test$scaled.shifted$stat + 2*autolav_sub2@Fit@npar
autolav_sub3@Fit@test$scaled.shifted$stat + 2*autolav_sub3@Fit@npar
autolav_sub4@Fit@test$scaled.shifted$stat + 2*autolav_sub4@Fit@npar
autolav_sub5@Fit@test$scaled.shifted$stat + 2*autolav_sub5@Fit@npar
autolav_sub6@Fit@test$scaled.shifted$stat + 2*autolav_sub6@Fit@npar
lav_cfa1@Fit@test$scaled.shifted$stat + 2*lav_cfa1@Fit@npar

#SRMRs:
fitmeasures(autolav)["srmr"]
fitmeasures(autolav_sub1)["srmr"]
fitmeasures(autolav_sub2)["srmr"]
fitmeasures(autolav_sub3)["srmr"]
fitmeasures(autolav_sub4)["srmr"]
fitmeasures(autolav_sub5)["srmr"]
fitmeasures(autolav_sub6)["srmr"]
fitmeasures(lav_cfa1)["srmr"]

#CFIs:
fitmeasures(autolav)["cfi"]
fitmeasures(autolav_sub1)["cfi"]
fitmeasures(autolav_sub2)["cfi"]
fitmeasures(autolav_sub3)["cfi"]
fitmeasures(autolav_sub4)["cfi"]
fitmeasures(autolav_sub5)["cfi"]
fitmeasures(autolav_sub6)["cfi"]
fitmeasures(lav_cfa1)["cfi"]

#TFIs:
fitmeasures(autolav)["tli"]
fitmeasures(autolav_sub1)["tli"]
fitmeasures(autolav_sub2)["tli"]
fitmeasures(autolav_sub3)["tli"]
fitmeasures(autolav_sub4)["tli"]
fitmeasures(autolav_sub5)["tli"]
fitmeasures(autolav_sub6)["tli"]
fitmeasures(lav_cfa1)["tli"]


# Model checks

auto_sntx_sub6_2 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ ia11*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ 0*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ 0*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ 0*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ 0*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ 0*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
res > 0
'
#^^^No lower bound on `ia11`.
autolav_sub6_2 <- lavaan(
	model=auto_sntx_sub6_2,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub6_2)
autolav_sub6_2@optim$fx

auto_sntx_sub6_3 <- '
#BMI1:
A1 =~ 1*BMI1
BMI1 ~~ res*BMI1 + BMI1
A1 ~~ 0*A1 + A1
#BMI2:
A2 =~ 1*BMI2
BMI2 ~~ res*BMI2 + BMI2
A2 ~~ 0*A2 + A2
A2 ~ b1*A1 + A1
#BMI3:
A3 =~ 1*BMI3
BMI3 ~~ res*BMI3 + BMI3
A3 ~~ 0*A3 + A3
A3 ~ b1*A2 + A2
#BMI4:
A4 =~ 1*BMI4
BMI4 ~~ res*BMI4 + BMI4
A4 ~~ 0*A4 + A4
A4 ~ b1*A3 + A3
#BMI5:
A5 =~ 1*BMI5
BMI5 ~~ res*BMI5 + BMI5
A5 ~~ 0*A5 + A5
A5 ~ b1*A4 + A4
#BMI6:
A6 =~ 1*BMI6
BMI6 ~~ res*BMI6 + BMI6
A6 ~~ 0*A6 + A6
A6 ~ b1*A5 + A5
#Constraints:
res > 0
'
#^^^No `ia11` at all.
autolav_sub6_3 <- lavaan(
	model=auto_sntx_sub6_3,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(autolav_sub6_3)
autolav_sub6_3@optim$fx

# Removing the lower bound on `ia11` does not adequately improve this model's fit; fixing `ia11` to zero made it slightly worse:
autolav_sub6@optim$fx
autolav_sub6_2@optim$fx
autolav_sub6_3@optim$fx


# Can we resolve the Haywood case in `lav_cfa1`?:
sntx_cfa1_2 <- '
F =~ 1*BMI1 + BMI2 + BMI3 + BMI4 + BMI5 + BMI6
F ~~ F
BMI1 ~~ u1*BMI1
BMI2 ~~ BMI2
BMI3 ~~ BMI3
BMI4 ~~ BMI4
BMI5 ~~ BMI5
BMI6 ~~ BMI6
u1 > 0
'  
lav_cfa1_2 <- lavaan(
	model=sntx_cfa1_2,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	estimator="WLSMV"
)
summary(lav_cfa1_2)
lav_cfa1_2@optim$fx
lav_cfa1@optim$fx
#^^^Haywood case can be resolved by placing a lower bound, but makes negligible difference.
