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
## partial load of steam turbine
STrLoad<-function(ratio){
  efratio<--3.3784*(ratio^2)+3.9324*ratio-0.1716
  efratio
}
## gas turbine related to air temperature ef=-0.002*x^2-0.1237*x+102.28 (R2=0.99949)
eGTrTa<-function(x){ef<-(-0.002*x^2-0.1237*x+102.28)}

## definition of stations
ocgt<-as.data.frame(read.csv2('CCGTv1.csv',sep=';',dec=',',header=TRUE),sep=';',header=TRUE,colClasses='numeric',dec='.')
est<-as.character(ocgt$Name[which(ocgt$Land==1)])
###
## definition of lapply function to be applied along est

estocgt<-est[which(ocgt$Cycle=='OCGT')]

n_loops_max<-seq(0,500,by=5)

## ST/GT from 53% to 80%
ratio_stgt_tot<-c(0.53,seq(0.55,0.8,by=0.05))

modProd<-lapply(estocgt,function(x){
  
## introduction of meteo data
setwd('/Users/usuario/Dropbox/Desertec/DocumentosActivos/meteo')
meteo<-as.data.frame(read.csv2(paste0('est',which(est==x),'.csv'),sep=';',header=TRUE,colClasses='numeric',dec='.'))
names(meteo)<-c('y','m','d','h','ghi','TempMed','dni','HumMed','P')
# wet or dry cooling
w<-ocgt$W[which(est==x)]
  
## introduction of coordinates of the ocgt studied
###
Lat<-as.numeric(as.character(ocgt$Latitude))[which(ocgt$Land==1)][which(est==x)]
Lon<-as.numeric(as.character(ocgt$Longitude))[which(ocgt$Land==1)][which(est==x)]

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
## Calculation for each loop

PerLoopTotal<-lapply(ratio_stgt_tot,function(ratio_stgt){

PerLoop<-lapply(n_loops_max,function(n_loops){
  
Psolar_th<-n_loops*4*Pcolfluid ## 4 Solar Collection Assemblies per loop
## we consider a 34.25% steam turbine efficiency from exergy graph book
ST<-ratio_stgt*as.numeric(as.character(ocgt$TG[which(ocgt$Name==x)]))
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
Pb<-((-27/2400)*ocgt$Elevation[which(ocgt$Name==x)]+100)/100*1013
## Natural gas consumption rate: NECESITAMOS ENCONTRARLO!!!!!!
#NG_cons_rate<-(16.85/266)*ocgt$TG[which(ocgt$Name==x)]#kg/s (266 MW GT (Palos de la Frontera book ocgt))
NG_cons_rate<-(20/260)*as.numeric(as.character(ocgt$TG[which(ocgt$Name==x)]))
efGTrTa<-eGTrTa(Ta)


potenciainicialst<-(ratio_stgt_tot[1]*as.numeric(as.character(ocgt$TG[which(ocgt$Name==x)])))
diferencia<-ST-potenciainicialst

### SCENARIO 1: solar boosting mode
## a) Natural gas consumption scenario 1
NG_cons_sc1<-NG_cons_rate*efGTrTa/100
## b) GT Power output scenario 1
Pgas_sc1<-GTrTaRHP(Ta,RH,P,Pb)*as.numeric(as.character(ocgt$TG[which(ocgt$Name==x)]))
## c) Integrable solar power scenario 1
Psol_int_sc1<-(-STrTaRHP(Ta)+1)*potenciainicialst
## remove negative values in Psol_int_sc1
Psol_int_sc1[which(Psol_int_sc1<0)]<-0
Psol_int_sc1<-Psol_int_sc1+diferencia
## from integrable to energy integrated
Psol_int_sc1i<-lapply(c(1:8760),function(t){
  if(Psolar_el[t]>Psol_int_sc1[t]){Psol_int_sc1[t]<-Psol_int_sc1[t]}else{Psol_int_sc1[t]<-Psolar_el[t]}
Total<-Psol_int_sc1[t]})
Psol_int_sc1<-do.call(c,Psol_int_sc1i)

## e) Dumping scenario 1
Dumping_sc1<-Psolar_el-Psol_int_sc1
Dumping_sc1[Dumping_sc1<0]<-0
## d) ST Power output scenario 1
Pst_sc1<-STrTaRHP(Ta)*potenciainicialst+Psol_int_sc1
## dd) Partial Load performance
ratio<-Pst_sc1/ST
efLoad<-c()
efLoad[which(ratio>0.6)]<-1
if(length(which(ratio<0.6))==0){a<-1}else{efLoad[which(ratio<0.6)]<-STrLoad(ratio[which(ratio<0.6)])}

Pst_sc1<-Pst_sc1*efLoad
## d) Total ccgt power output scenario 1
Pccgt_sc1<-Pst_sc1+Pgas_sc1
if(w==1){Pccgt_sc1<-Pccgt_sc1}else{Pccgt_sc1<-Pccgt_sc1*0.9}
## f) Overall efficiency
ef_ocgt_sc1<-Pccgt_sc1/(NG_cons_sc1*PCI)

## Electrical Power
Aux_con<-1.1
ef_sf<-0.52

Rad<-600
Pe_fixed<-n_loops*LS3oparea*4*Rad*ef_st*ef_sf/(Aux_con*1000) # in KW
# precio conversion OCGT en CCGT
p_conv<-605.4 # USD/kWe
pot<-(ST+as.numeric(as.character(ocgt$TG[which(ocgt$Name==x)])))
Cost_fixed<-p_conv*pot*1000+15*pot+0.85*sum(Pccgt_sc1)
# cost of hybridation in dollars 
p_conv_solar<-3956
Cost_fixed_solar<-p_conv_solar*Pe_fixed

## Financial Model
## discount rate
dr<-0.04
infl<-0.02
years<-25
## Positive cash flows

## electricity price for <5% solar share
ep<-0.027214 #USD/kWh

infla<-lapply(c(1:(years-1)),function(x){infl<-(1+infl)^x}); inflation<-do.call(c,infla);inflation<-c(1,inflation)
CF_sell<-sum(Pst_sc1*ep*1000)*inflation


## Negative cash flows
CF_fixed<--49.9*Pe_fixed*inflation
CF_var<- -0.7*sum(Psolar_el)*inflation
CF_st_fixed<--18*ST*100
CF_st_var<--1.5*sum(Pst_sc1)
# Total cash flows
CF0<-Cost_fixed+Cost_fixed_solar
CF<-CF_sell+CF_fixed+CF_var+CF_st_fixed+CF_st_var

residual_cost<-0.1
CF[25]<-CF[25]+residual_cost*CF0
npv<-sum(do.call(c,lapply(c(1:years),function(t){a<-CF[t]/((1+dr)^t)})))-CF0
CF<-c(-CF0,CF)
irr<-cf(CF,dr)$irr[1]

npv<-rep(npv,8760)
irr<-rep(irr,8760)

# levelized cost of electricity
CF<-CF[-1]
num<-lapply(c(1:years),function(t){(CF[t])/((1+dr)^t)});num<-do.call(c,num)
num<-sum(num)+CF0
den<-lapply(c(1:years),function(t){sum(Pst_sc1)*1000/((1+dr)^t)});den<-do.call(c,den);den<-sum(den)
lcoe<-num/den
lcoe<-rep(lcoe,8760)


Total<-data.frame(iTs,month,day,julianday,hour,julianhour,DNI,Ta,RH,P,
                  incang,modincang,optef,Lossarea,Pcolenv,Psuncol,Pcolfluid,Psolar_th,Psolar_el,efGTrTa,
                  NG_cons_sc1,Pgas_sc1,Psol_int_sc1,Dumping_sc1,Pccgt_sc1,Pst_sc1,ef_ocgt_sc1,
                  npv,irr,lcoe)

Total})

PerLoop

})

num<-ocgt$Number[which(est==x)]
old<-getwd()
setwd('/Users/usuario/Desktop/OCGT-ISCC')
save(PerLoopTotal,file=paste('Todo',num,'.RData',sep=''))
setwd(old)
})

num<-which(ocgt$Cycle=='OCGT')
setwd('/Users/usuario/Desktop/OCGT-ISCC')
load('Todo4.RData')
## Escenario de operación: horas centrales del día

horas_operacion<-seq(11,17,1)
horas<-lapply(horas_operacion,function(x){
  t<-which(PerLoopTotal[[1]][[1]]$hour==x)
})
horas<-sort(do.call(c,horas))

Annual<-lapply(num,function(est){
  load(paste('Todo',est,'.RData',sep=''))
tt<-lapply(c(1:7),function(conf_st){ 
ttt<-lapply(c(1:length(n_loops_max)),function(lazos){
    x<-PerLoopTotal[[conf_st]][[lazos]][horas,]
    Pocgt_sc0<-sum(x$Pgas_sc0)
    Pccgt_sc1<-sum(x$Pgas_sc1)+sum(x$Pst_sc1)
    Psol_int_sc1<-sum(x$Psol_int_sc1)
    Dumping_sc1<-sum(x$Dumping_sc1)
    NG_cons_sc0<-sum(x$NG_cons_sc0)
    Pgas_sc1<-sum(x$Pgas_sc1)
    Pst_sc1<-sum(x$Pst_sc1)
    Pgt_sc0<-sum(x$Pgas_sc0)
    lcoe<-x$lcoe[1]
    irr<-x$irr[1]
    npv<-x$npv[1]
    solar_share<-Psol_int_sc1/Pccgt_sc1
    Total<-c(Pocgt_sc0,Pst_sc1,Pgt_sc0,Pccgt_sc1,Psol_int_sc1,Dumping_sc1,
             solar_share,lcoe,irr,npv)
  })
  Total<-data.frame(do.call(rbind,ttt))
  names(Total)<-c('Pocgt_sc0','Pst_sc1','Pgt_sc0','Pocgt_sc1','Psol_int_sc1','Dumping_sc1',
                  'solar_share','lcoe','irr','npv')
  return(Total)
})
tt
})
# each element of Annual is for a station, 7 steam turbine configurations and 61 different solar field configurations.
save(Annual,file='Annual.RData')

## Selecciono el número de lazos para el cual el dumping es mayor al solar integrable
K<-lapply(c(1:18),function(stat){
  x<-Annual[[stat]]
  M<-lapply(c(1:7),function(configu){
    y<-x[[configu]]
    t<-which(y$Dumping_sc1>y$Psol_int_sc1)[1]-1
    if(is.na(t)==FALSE){nloop<-t
                        nloopreal<-n_loops_max[t]; a<-data.frame(y[nloop,],nloopreal)}else{nloop<-max(n_loops_max); 
                                                                      a<-data.frame(y[length(n_loops_max),],nloop)}
    names(a)<-c('Pocgt','Pst','Pgt','Piscc','Psol_int','Dumping','solar_share','lcoe','irr','npv')
    a
  });M<-do.call(rbind,M)
  M
})
## Unidades
KK<-lapply(K,function(x){
  m<-data.frame(x[,c(1,4,5,6,2)]/1000,x[,c(7,8)]*100,x[,c(9,11)])
  names(m)<-c('Pocgt','Piscc','Psol','Dumping','Pst','sol_share','lcoe','irr','nloop')
  m
})

R_053<-data.frame(num,do.call(rbind,lapply(KK,function(x){x[1,]})))
R_055<-data.frame(num,do.call(rbind,lapply(KK,function(x){x[2,]})))
R_060<-data.frame(num,do.call(rbind,lapply(KK,function(x){x[3,]})))
library(xtable)
print(xtable(R_053,digits=c(1,1,1,1,1,1,1,2,2,2,0)),include.rownames=FALSE)
print(xtable(R_055,digits=c(1,1,1,1,1,1,1,2,2,2,0)),include.rownames=FALSE)
print(xtable(R_060,digits=c(1,1,1,1,1,1,1,2,2,2,0)),include.rownames=FALSE)

Pot<-data.frame(rbind(apply(R_053,2,sum)[-c(1,6,7,8,9)],apply(R_060,2,sum)[-c(1,6,7,8,9)],
                 apply(R_ccgt,2,sum)[c(2,5,6,7,11)]),row.names=c('OCGT-R053','OCGT-R060','CCGT'))
names(Pot)<-c('Y_origin','P_iscc','Psol','Dumping','nloop')
xtable(Pot,digits=c(1,1,1,1,1,0))

