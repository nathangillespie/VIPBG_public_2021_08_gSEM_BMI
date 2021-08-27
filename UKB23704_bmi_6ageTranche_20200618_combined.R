library(GenomicSEM)
library(lavaan)

# Read in S and V matrices
# Autoregression models
# 1-factor EFA
# Fit indices



# Read in S and V matrices

load("BMI_6ageTranche_20200605.RData")


S 					<- LDSCoutput_v3$S
S 					<- as.matrix(Matrix::nearPD(S)$mat)
rownames(S) 		<- colnames(S)
eigen(S,T,T)$values
eigen(LDSCoutput_v3$V,T,T)$values
V 					<- LDSCoutput_v3$V
dim(LDSCoutput_v3$V)
rownames(V) 		<- as.character(1:21)
colnames(V) 		<- rownames(V)
W 					<- chol2inv(chol(LDSCoutput_v3$V))
#dataDWLS 			<- mxData(S, numObs = 2, means = NA, type = "acov", acov=diag(diag(V)), fullWeight=W)

S
 #           BMI1      BMI2      BMI3      BMI4      BMI5      BMI6
 # BMI1 0.2422834 0.2537639 0.2527534 0.2437357 0.2321645 0.2228504
 # BMI2 0.2537639 0.2665734 0.2662147 0.2571411 0.2395505 0.2310082
 # BMI3 0.2527534 0.2662147 0.2673344 0.2633383 0.2374354 0.2302735
 # BMI4 0.2437357 0.2571411 0.2633383 0.2859213 0.2385714 0.2337707
 # BMI5 0.2321645 0.2395505 0.2374354 0.2385714 0.2441764 0.2303000
 # BMI6 0.2228504 0.2310082 0.2302735 0.2337707 0.2303000 0.2187359

 cov2cor(S)
 #           BMI1      BMI2      BMI3      BMI4      BMI5      BMI6
 # BMI1 1.0000000 0.9985263 0.9931336 0.9260489 0.9545137 0.9680361
 # BMI2 0.9985263 1.0000000 0.9972319 0.9314076 0.9389377 0.9566629
 # BMI3 0.9931336 0.9972319 1.0000000 0.9524964 0.9293220 0.9522619
 # BMI4 0.9260489 0.9314076 0.9524964 1.0000000 0.9029073 0.9347744
 # BMI5 0.9545137 0.9389377 0.9293220 0.9029073 1.0000000 0.9965107
 # BMI6 0.9680361 0.9566629 0.9522619 0.9347744 0.9965107 1.0000000





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
lavTestLRT(autolav,autolav_sub1,autolav_sub2,autolav_sub3,autolav_sub4,autolav_sub5,autolav_sub6)
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
lav_cfa1@Fit@test$scaled.shifted$stat + lav_cfa1@Fit@npar

#SRMRs:
fitmeasures(autolav)["srmr"]
fitmeasures(autolav_sub1)["srmr"]
fitmeasures(autolav_sub2)["srmr"]
fitmeasures(autolav_sub3)["srmr"]
fitmeasures(autolav_sub4)["srmr"]
fitmeasures(autolav_sub5)["srmr"]
fitmeasures(autolav_sub6)["srmr"]
fitmeasures(lav_cfa1)["srmr"]
