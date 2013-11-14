library(ggplot2)
fns<-Sys.glob("*/mypredictionresults.csv")
allresults<-matrix( rep( NA , 12* length(fns) ) , ncol=12 )
ct<-1
for ( fn in fns ) {
  mycross<-read.csv(fn)$Cross;
  print(mean(mycross,na.rm=T)) ;
  plot( hist( mycross ) , main=fn) ;
  allresults[ct,  ]<-mycross
  ct<-ct+1
#  Sys.sleep(4)
}
subj1<-data.frame( RateCorrect=allresults[1, ], id='1' )
subj2<-data.frame( RateCorrect=allresults[2, ], id='2' )
subj3<-data.frame( RateCorrect=allresults[3, ], id='3' )
subj4<-data.frame( RateCorrect=allresults[4, ], id='4' )
subj5<-data.frame( RateCorrect=allresults[5, ], id='5' )
subj6<-data.frame( RateCorrect=allresults[6, ], id='6' )
#and combine into your new data frame vegLengths
subjs <- rbind(subj1, subj2, subj3, subj4, subj5, subj6 )
#now make your lovely plot
pdf('haxby_summary_results.pdf')
myplot <- ggplot(subjs, aes(RateCorrect, fill = id)) + geom_density(alpha = 0.2) + labs(title="Haxby 2001 Predictions - 6 Subjects, Majority Vote, 9 categories")
plot( myplot )
dev.off()
mydf<-data.frame(allresults)
colnames(mydf)<-paste("Run",c(1:12))
write.csv(mydf,"haxby_results_basic.csv")
