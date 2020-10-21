library(tseries)
library(forecast)
library(magrittr)


set.seed(0)
#seed beállítása, hogy ugyanazok a pszeudorandom számaink legyenek

#1. Tárolófüggvény
dgp_parameterek<-function(tslength,nofts,noise,p,d,q){
  return(c(tslength,nofts,noise,p,d,q))
}
#eltárolja a szimuláció paramétereit



#2. Adatgeneráló függvény
ts_gen<-function(x){
  tslength<-x[1]
  nofts<-x[2]
  noise<-x[3]
  p<-x[4]
  d<-x[5]
  q<-x[6]
  #Kiolvassa az input vektorból az egyes paramétereket
  idosorok<-array(numeric(),c(tslength,nofts)) 
  #létrehoz egy most még üres mátrixot az idősoroknak
  #a mátrix tslength x nofts dimenziós
  #vagyis minden oszlop egy idősor lesz
  for (i in 1:nofts) {
    idosorok[1:tslength,i]<-arima.sim(model=list(order=c(p,d,q),
                                                 ar=runif(p,min=-0.5,max=0.5),ma=runif(q,min=-0.5,max=0.5)),
                                      n=tslength,sd=noise)}
  #legenerálom a megfelelő hosszúságú és számú idősort
  #a theták és phik a robusztusság miatt -0.5 és 0.5 közötti
  #egyenletes eloszlású véletlen számok
  return(idosorok)
}


#3. Illesztési paraméterek tárolása
fit_input<-function(maxp,maxq){
  return(c(maxp,maxq))
}

#4. Modellillesztés
fitmodel<-function(x,ts){
  maxp<-x[1]
  maxq<-x[2]
  akaike<-array(numeric(),c(maxp+1,maxq+1))
  schwarz<-array(numeric(),c(maxp+1,maxq+1))
  #definiálja az egyes modellekhez tartozó aic és bic értékek mátrixait
  #+1 sor és oszlop, mert a p és q 0 is lehet
  d <- ndiffs(ts,test="kpss")
  #meghatározza d-t
  for (i in 0:maxp) {
    for (j in 0:maxq) {
      fit<-arima(ts,order = c(i,d,j))
      akaike[i+1,j+1]<-AIC(fit)
      schwarz[i+1,j+1]<-BIC(fit)
    }
  }
  #minden (p,q) paraméterpárra fitteli az arimat
  #és kiszámolja erre az információs kritériumokat
  
  aic_opt <- which(akaike == min(akaike), arr.ind = TRUE)-1
  bic_opt <- which(schwarz == min(schwarz), arr.ind = TRUE)-1
  #megnézi, hogy az akaike és scwarz mátrixoknak melyik eleme a legkisebb
  #ennek az elemnek a sorát és oszlopát kimenti
  #mivel p=0 és q=0-val kezdünk ezért -1
  #ez lesz az aic illetve bic által optimálisnak talált (p,q) paraméterpár
  optimalis <- list(aic_opt,bic_opt)
  return(optimalis)
}


#########################################


#5. A szimulációs függvény
result_table<-function(pdq,max){
  aic_grid<-array(numeric(),c(5,5))
  bic_grid<-array(numeric(),c(5,5))
  noise<-c(0.01,0.1,1,10,100)
  tslength<-c(100,500,1000,5000,10000)
  #definiálom a két gridet és paramétereit
  for (i in 1:5) {
    for (j in 1:5) {
      parameterek <- c(tslength[i],1000,noise[j],pdq[1],pdq[2],pdq[3])
      #Összeszedem a grid i,j elemének megfelelő paraméterek értékeit egy vektorba
      ts <- ts_gen(parameterek)
      #legenerálom az idősorokat
      aicszamlalo <- 0
      bicszamlalo <- 0
      for (k in 1:1000){
        ered <- try(fitmodel(max,ts[,k]),silent=TRUE)
        #optimális modell keresése
        #a try() functionre az arima függvény numerikus megoldásakor 
        #adódó esetleges hibája miatt van szükség
        if(is.numeric(ered[[1]][1])==TRUE){
          #azt vizsgálja meg, hogy nincs-e hibaüzenet az outputban
          if(ered[[1]][1] == pdq[1] & ered[[1]][2] == pdq[3]){
            aicszamlalo <- aicszamlalo+1}
          if(ered[[2]][1] == pdq[1] & ered[[2]][2] == pdq[3]){
            bicszamlalo <- bicszamlalo+1}
          #megnézi, hogy aic, bic eltalálta-e az adatgeneráló folyamat rendjét 
        }
      }
      aic_grid[i,j] <- aicszamlalo/1000
      bic_grid[i,j] <- bicszamlalo/1000
      #az adott (zaj,mintaelemszám) specifikációra lefuttatott eredményt elmenti
    }
  }
  output <- list(aic_grid,bic_grid)
  #egyszerre jeleníti meg az aic és bic eredményeit
  return(output)
}



######### A függvény tesztelése ###########
#VIGYÁZAT! Egy-egy futtatás nagyon sok idő!

#1. specifikáció: sima egyszerű AR
spec1<-result_table(c(0,0,1),c(3,3))

#2. specifikáció: egyszerű ARMA
spec2<-result_table(c(1,0,1),c(3,3))

#3. specifikáció: komplex ARIMA
spec3<-result_table(c(2,1,1),c(4,4))

#4. specifikáció: még komplexebb ARIMA
spec4<-result_table(c(4,2,3),c(5,5))