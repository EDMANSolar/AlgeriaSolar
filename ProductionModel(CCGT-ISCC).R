## load of required libraries
library(raster)
library(solaR)
library(financial)
## introduce the working directory
dir<-as.character('/Users/usuario/Dropbox/Desertec/DocumentosActivos/CCGTAlgeria')
setwd(dir)
## definition of the daily temporal series. In this case 2005 daily time series (dTs)
## is defined with fBTd (solaR function)
dTs<-fBTd(mode='serie',start='01-01-2005',end='31-12-2005',format='%d-%m-%Y')
## definition of the intradaily temporal series (iTs) with as.POSIXct
iTs<-seq(as.POSIXct('2005-01-01 00:30:00'),as.POSIXct('2005-12-31 23:30:00'),by='hour')
## definition of month, day,hour, julian day and julian hour
dom<-c(31,28,31,30,31,30,31,31,30,31,30,31)
month<-rep(seq(1:12),dom*24)
df<-lapply(dom,function(x){t<-as.numeric(seq(1:x))})
day<-rep(do.call(c,df),each=24)
hour<-rep(1:24,365)
julianday<-rep(seq(1:365),each=24)
julianhour<-seq(1:8760)
## Solar field specifications
## reflectividad
LS3reflec<-0.94
## transmisividad
LS3trans<-0.955
## interceptation factor
LS3intf<-0.997
## absortivity
LS3abs<-0.955
## peak efficiency
nopico<-LS3reflec*LS3trans*LS3intf*LS3abs
LS3width<-5.76
LS3focdis<-1.71
LS3length<-99
LS3oparea<-545
mirrorclean<-0.98
##fluid inlet temperature
HTFin<-293
HTFout<-393
HTFm<-343
##other data
PCI<-39.900 ## PCI of natural gas [kJ/kg]

## Supposing cummulative efficiency curves: GTrTaRHP (power)
GTrTaRHP<-function(Ta,RH,P,Pb){
  efTa<-(-0.5024*Ta+107.536)/100
  efRH<-((1.05/90)*RH+99.3)/100
  efP<-((27/25)*((P/Pb)*100-100)+100)/100
  ef<-efTa*efRH*efP
  ef}
## Supposing cummulative efficiency curves: STrTaRHP (power)
STrTaRHP<-function(Ta){
  #efTa<-(-0.053185*Ta+100.79)/100 ## Shi.Agnew.ea2010
  efTa<-(6e-4*Ta^2-0.1579*Ta+102.24)/100 ## Excel Javier
  #efRH<-(-8E-06*RH^2 + 0.0125*RH + 99.251)/100 
  #efP<-((27/25)*((P/Pb)*100-100)+100)/100 
  ef<-efTa#*efRH*efP
  ef}

## definition of stations
ccgt<-as.data.frame(read.csv2('CCGTv1.csv',sep=';',dec=',',header=TRUE),sep=';',header=TRUE,colClasses='numeric',dec='.')
est<-as.character(ccgt$Name[which(ccgt$Land==1)])
###
## definition of lapply function to be applied along est

estccgt<-est[which(ccgt$Cycle=='CCGT')]

for(n_loops in 1:50){

modProd<-lapply(estccgt,function(x){
## introduction of coordinates of the CCGT studied
###
Lat<-as.numeric(as.character(ccgt$Latitude))[which(ccgt$Land==1)][which(est==x)]
Lon<-as.numeric(as.character(ccgt$Longitude))[which(ccgt$Land==1)][which(est==x)]
#Pot<-tapply(ccgt$Gross_Power[which(ccgt$Room==1)],ccgt$Same[which(ccgt$Room==1)],sum)
#Pot<-Pot[!is.na(Pot)]
## introduction of meteo data
setwd('/Users/usuario/Dropbox/Desertec/DocumentosActivos/meteo')
meteo<-as.data.frame(read.csv2(paste0('est',which(est==x),'.csv'),sep=';',header=TRUE,colClasses='numeric',dec='.'))
names(meteo)<-c('y','m','d','h','ghi','TempMed','dni','HumMed','P')
## Apparent movement of the Sun from the Earth: incidence angle
sol<-calcSol(Lat,local2Solar(iTs,Lon),sample='hour',EoT=TRUE,method='michalsky')
incang<-r2d(acos(as.numeric(fTheta(sol,modeTrk='horiz')$cosTheta)[25:8784]))
incang[which(is.na(incang))]<-0
## Definition of modificator of incidence angle as per Gonz?lez.Zarza.ea2001
which(incang>80)
## No angle higher than 80?, which would mean modincang<-0
modincang<-1-2.23073e-4*(incang)-1.1e-4*(incang^2)+3.18596e-6*(incang^3)-4.85509e-8*(incang^4)
modincang[which(incang>80)]<-0
## Optical efficiency
optef<-modincang*nopico
## Loss area per collector
Lossarea<-LS3width*(LS3focdis+((LS3focdis*(LS3width^2))/(48*(LS3focdis^2))))*tan(d2r(incang))
## Specific heat inpu8t to absorber (kWth/m2)
#spheatin<-meteo$dni*optef*mirrorclean/1000
## Heat loss from collector to environment (kWth per collector)
Pcolenv<-(0.00154*(HTFm-meteo$TempMed)^2+0.2021*(HTFm-meteo$TempMed)-24.899+
  ((0.00036*(HTFm-meteo$TempMed)^2+0.2029*(HTFm-meteo$TempMed)+24.899)*(meteo$dni/900)*
     cos(d2r(incang))))*LS3length/1000
## Potential thermal Power (kWth per collector)
Psuncol<-LS3oparea*meteo$dni*cos(d2r(incang))/1000
## Thermal power from collector to fluid
Pcolfluid<-(LS3oparea-Lossarea)*meteo$dni*cos(d2r(incang))*optef*mirrorclean/1000-Pcolenv
Pcolfluid[which(Pcolfluid<0)]<-0
Psolar_th<-n_loops*4*Pcolfluid ## 4 Solar Collection Assemblies per loop
## we consider a 34.25% steam turbine efficiency from exergy graph book
ef_st<-0.3425
ef_hrsg<-0.95
Psolar_el<-Psolar_th*ef_st*ef_hrsg/1000 ## in MW
## remove NA
Psolar_el[which(is.na(Psolar_el))]<-0
## Annual thermal power from collector to fluid kWhth
Pcolannual<-sum(Pcolfluid,na.rm=TRUE)
## Losses of power of combined cycle, gas turbine and steam turbine
## Definition of inputs
Ta<-meteo$TempMed
RH<-meteo$HumMed
P<-meteo$P
DNI<-meteo$dni
## base pressure
Pb<-((-27/2400)*ccgt$Elevation[which(ccgt$Name==x)]+100)/100*1013
## gas turbine related to air temperature ef=-0.002*x^2-0.1237*x+102.28 (R2=0.99949)
eGTrTa<-function(x){ef<-(-0.002*x^2-0.1237*x+102.28)}
## Natural gas consumption rate: NECESITAMOS ENCONTRARLO!!!!!!
#NG_cons_rate<-(16.85/266)*ccgt$TG[which(ccgt$Name==x)]#kg/s (266 MW GT (Palos de la Frontera book CCGT))
NG_cons_rate<-(20/260)*as.numeric(as.character(ccgt$TG[which(ccgt$Name==x)]))
efGTrTa<-eGTrTa(Ta)

### SCENARIO 1: solar boosting mode
## a) Natural gas consumption scenario 1
NG_cons_sc1<-NG_cons_rate*efGTrTa/100
## b) GT Power output scenario 1
Pgas_sc1<-GTrTaRHP(Ta,RH,P,Pb)*as.numeric(as.character(ccgt$TG[which(ccgt$Name==x)]))
## c) Integrable solar power scenario 1
Psol_int_sc1<-(-STrTaRHP(Ta)+1)*as.numeric(as.character(ccgt$TV[which(ccgt$Name==x)]))
## remove negative values in Psol_int_sc1
Psol_int_sc1[which(Psol_int_sc1<0)]<-0
## from integrable to energy integrated
Psol_int_sc1i<-lapply(c(1:8760),function(t){
  if(Psolar_el[t]>Psol_int_sc1[t]){Psol_int_sc1[t]<-Psol_int_sc1[t]}else{Psol_int_sc1[t]<-Psolar_el[t]}
Total<-Psol_int_sc1[t]})
Psol_int_sc1<-do.call(c,Psol_int_sc1i)

## e) Dumping scenario 1
Dumping_sc1<-Psolar_el-Psol_int_sc1
Dumping_sc1[Dumping_sc1<0]<-0
## d) ST Power output scenario 1
Pst_sc1<-STrTaRHP(Ta)*as.numeric(as.character(ccgt$TV[which(ccgt$Name==x)]))+Psol_int_sc1
## d) Total CCGT power output scenario 1
Pccgt_sc1<-Pst_sc1+Pgas_sc1
## f) Overall efficiency
ef_ccgt_sc1<-Pccgt_sc1/(NG_cons_sc1*PCI)

## Electrical Power
Aux_con<-1.1
ef_sf<-0.52

Rad<-600
Pe_fixed<-n_loops*LS3oparea*4*Rad*ef_st*ef_sf/(Aux_con*1000) # in KW
# precio conversion CCGT en ISCC
p_conv<-3392 # USD/kWe
Cost_fixed<-p_conv*Pe_fixed # cost of hybridation in dollars 

## Financial Model
## discount rate
dr<-0.04
infl<-0.02
years<-25
## Positive cash flows

## electricity price for <5% solar share
ep<-0.027214 #USD/kWh
infla<-lapply(c(1:(years-1)),function(x){infl<-(1+infl)^x}); inflation<-do.call(c,infla);inflation<-c(1,inflation)
CF_sell<-sum(Psol_int_sc1*ep*1000)*inflation

## Negative cash flows
CF_fixed<--49.9*Pe_fixed*inflation
CF_var<- -0.7*sum(Psolar_el)*inflation

# Total cash flows
CF0<-Cost_fixed
CF<-CF_sell+CF_fixed+CF_var
residual_cost<-0.1
CF[25]<-CF[25]+residual_cost*CF0
# Residual cost
npv<-sum( do.call(c,lapply(c(1:years),function(t){a<-CF[t]/((1+dr)^t)})))-CF0

CF<-c(-CF0,CF)

irr<-cf(CF,dr)$irr[1]

npv<-rep(npv,8760)
irr<-rep(irr,8760)

# levelized cost of electricity
num<-lapply(c(1:years),function(t){(CF_fixed[t]+CF_var[t])/((1+dr)^t)});num<-do.call(c,num)
num[25]<-CF_fixed[25]+CF_var[25]+residual_cost*CF0
num<--CF0-sum(num)
den<-lapply(c(1:years),function(t){sum(Psol_int_sc1)*1000/((1+dr)^t)});den<-do.call(c,den);den<-sum(den)
lcoe<-num/den
lcoe<-rep(lcoe,8760)
  
  
### SCENARIO 0: CCGT
## a) Natural gas consumption
NG_cons_sc0<-NG_cons_rate*efGTrTa/100
## b) GT power output scenario 0
Pgas_sc0<-GTrTaRHP(Ta,RH,P,Pb)*as.numeric(as.character(ccgt$TG[which(ccgt$Name==x)]))
## c) ST power output sceneario 0
Pst_sc0<-STrTaRHP(Ta)*as.numeric(as.character(ccgt$TV[which(ccgt$Name==x)]))
## d) CCGT power output scenario 0
Pccgt_sc0<-Pgas_sc0+Pst_sc0
## e) Solar power integration
Psol_int_sc0<-rep(0,8760)
## g) Overall efficiency
ef_ccgt_sc0<-Pccgt_sc0/(NG_cons_sc0*PCI)

Total<-data.frame(iTs,month,day,julianday,hour,julianhour,DNI,Ta,RH,P,
                  incang,modincang,optef,Lossarea,Pcolenv,Psuncol,Pcolfluid,Psolar_th,Psolar_el,efGTrTa,
                  NG_cons_sc1,Pgas_sc1,Psol_int_sc1,Dumping_sc1,Pccgt_sc1,Pst_sc1,ef_ccgt_sc1,
                  NG_cons_sc0,Pgas_sc0,Pst_sc0,Pccgt_sc0,Psol_int_sc0,ef_ccgt_sc0,npv,irr,lcoe)

Total})
old<-getwd()
setwd('/Users/usuario/Dropbox/Desertec/DocumentosActivos/Resultados/CCGT-ISCC')
save(modProd,file=paste('Todo',n_loops,'.RData',sep=''))
setwd(old)
}

setwd('/Users/usuario/Dropbox/Desertec/DocumentosActivos/Resultados/CCGT-ISCC')
load('Todo1.RData')
horas_operacion<-seq(11,17,1)
horas<-lapply(horas_operacion,function(h){
  which(modProd[[1]]$hour==h)})
horas<-sort(do.call(c,horas))
setwd('/Users/usuario/Dropbox/Desertec/DocumentosActivos/Resultados/CCGT-ISCC')
Annual<-lapply(c(1:50),function(x){
  load(paste('Todo',x,'.RData',sep=''))
 t<- lapply(modProd,function(x){
    Pccgt_sc0<-sum(x$Pccgt_sc0[horas])
    Pccgt_sc1<-sum(x$Pgas_sc1[horas])+sum(x$Pst_sc1[horas])
    Psol_int_sc1<-sum(x$Psol_int_sc1[horas])
    Dumping_sc1<-sum(x$Dumping_sc1[horas])
    NG_cons_sc0<-sum(x$NG_cons_sc0[horas])
    Pgas_sc1<-sum(x$Pgas_sc1[horas])
    Pst_sc1<-sum(x$Pst_sc1[horas])
    Pst_sc0<-sum(x$Pst_sc0[horas])
    Pgt_sc0<-sum(x$Pgas_sc0[horas])
    solar_share<-(Psol_int_sc1/Pccgt_sc1)*100
    irr<-x$irr[1]
    lcoe<-x$lcoe[1]
    Total<-c(Pccgt_sc0,Pst_sc0,Pgt_sc0,Pccgt_sc1,Psol_int_sc1,Dumping_sc1,NG_cons_sc0,Pgas_sc1,Pst_sc1,solar_share,irr,lcoe)
  })
  Total<-data.frame(do.call(rbind,t))
  names(Total)<-c('Pccgt_sc0','Pst_sc0','Pgt_sc0','Pccgt_sc1','Psol_int_sc1',
                  'Dumping_sc1','NG_cons_sc0','Pgas_sc1','Pst_sc1','solar_share','irr','lcoe')
  return(Total)
})
save(Annual,file='Annual.RData')

## Selecciono el número de lazos para el cual el dumping es mayor al solar integrable
K<-lapply(c(1:3),function(st_ccgt){
  m<-lapply(c(1:50),function(nloop){t<-Annual[[nloop]][st_ccgt,]$Dumping_sc1>Annual[[nloop]][st_ccgt,]$Psol_int_sc1;t})
  m<-do.call(c,m)
  m<-which(m==TRUE)[1]-1
  m<-data.frame(Annual[[m]][st_ccgt,],m)})

## Unidades
KK<-lapply(K,function(x){
  m<-data.frame(x[,c(1,2,3,4,5,6)]/1000,x[,10],x[,12]*100,x[,11],x[,13])
  names(m)<-c('Pccgt','Pst','Pgt','Piscc','Psol','Dumping','sol_share','lcoe','irr','nloop')
  m
})
est<-c(1,2,3)
R_ccgt<-data.frame(est,do.call(rbind,KK))
print(xtable(R_ccgt,digits=c(1,0,1,1,1,1,1,1,2,2,2,0)),include.rownames=FALSE)


load('Annual.RData')

LCOE<-do.call(rbind,Annual)$lcoe
IRR<-do.call(rbind,Annual)$irr/100
Station<-rep(c('Station 1','Station 2','Station 3'),50)
Loop<-rep(c(1:50),each=3)
LCOE<-data.frame(LCOE,IRR,Station,Loop)
LCOE$Station<-as.factor(LCOE$Station)
trellis.device(pdf, file='LCOEccgt-iscc.pdf')
xyplot(LCOE~Loop,groups=Station,data=LCOE,auto.key=list(space='right', points=TRUE, lines=FALSE))
dev.off() 


#############################################################################################################
## Gráficos de resultados ###################################################################################
#############################################################################################################

## 
setwd('/Users/usuario/Dropbox/Desertec/DocumentosActivos/Resultados/CCGT-ISCC')
load('Todo10.RData')
Lazo1<-modProd[[3]]

Sys.setlocale("LC_TIME", 'C')
P0<-Lazo1$Pccgt_sc0[25:48]
P1<-Lazo1$Pccgt_sc1[25:48]
#P2<-Lazo1$Pccgt_sc2[4273:4296]
D1<-Lazo1$Dumping_sc1[25:48]
Ta<-Lazo1$Ta[25:48]
day<-seq(as.POSIXct('2005-01-02 00:00:00'),as.POSIXct('2005-01-02 23:00:00'),by='hour')

obj1<-xyplot(P1~day,ylim=c(1180,1220),type='l',col='mediumblue',ylab='P (MW)',
             auto.key=FALSE,panel=function(x, y, ...){
               panel.xyplot(x, y, ...)
               panel.abline(v=c(day[12],day[18]),col="black",alpha=0.5,lwd=1.2)})
#obj2<-xyplot(P2~day,ylim=c(405,437),type='l',col='green',lwd=1.2)
obj3<-xyplot(D1~day,type='l',ylim=c(0,20))
obj4<-xyplot(Ta~day,type='l',ylim=c(0,20),col='indianred')
obj5<-xyplot(P0~day,ylim=c(1180,1220),type='l',col='orange',lwd=1.2)
trellis.device(pdf, file='diaconcreto.pdf')
update(doubleYScale(obj1, obj3,text = c("Piscc_1",'Dumping','Ta','Piscc_0')),
       par.settings = simpleTheme(col = c('mediumblue','black','indianred','orange'), lty = c(1,1,1,1)))+
  update(doubleYScale(obj1, obj4,text = c("Piscc_1",'Ta')),
         par.settings = simpleTheme(col = c('black','black','indianred','black'), lty = c(1,1,1,1)))+
  obj5
dev.off()


M<-data.frame(do.call(rbind,Annual)[,c(5,6)],rep(c('St 1','St 2','St 3'),50),rep(c(1:50),each=3))
names(M)<-c('Psol','Dumping','St','Loops')

trellis.device(pdf, height=4, width=6,file='Dumping-SolarInt_CCGT.pdf')
xyplot(Dumping~Loops,groups=St,superpose=T,
       auto.key=list(c('a','b','c'),space='right', points=TRUE, lines=FALSE),
       data=M,ylab='MWh',cex=0.4, type='p')+
  xyplot(Psol~Loops,groups=St,superpose=T,
         auto.key=list(space='right', points=TRUE, lines=FALSE),
         data=M,ylab='MWh',cex=0.4, type='p')
dev.off()
