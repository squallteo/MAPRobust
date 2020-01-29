

decision = results[,1:4]
median_est = results[,5:8]
w_post = results[,9]


tt1 <- data.frame(TrueCtrlMean = muvec_c[m], Prop = colMeans(decision), Method = c("rMAP","Alt1", "Alt2", "Alt3"))
tt2 <- data.frame(TrueCtrlMean = muvec_c[m], Prop = colMeans(decision), Method = c("rMAP","Alt1", "Alt2", "Alt3"))


if(m==1) outdt_decision = ttt else outdt_decision = rbind(outdt_decision,ttt)

# 
# if(effsize != 0){
#   ggplot(data=outdt,aes(x=TrueCtrlMean,y=Prop,color=Method)) + geom_line(size = 1.5) + 
#     scale_x_continuous(breaks = muvec_c, name = "True Control Mean") + 
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Probability of Success") + 
#     geom_vline(xintercept = -49.9, linetype=2, size=1) + 
#     theme(axis.title = element_text(face="bold",size=15),
#           axis.text = element_text(size=12),
#           legend.title=element_text(size=15,face="bold"),
#           legend.text=element_text(size=12)
#     )
# }
# 
# if(effsize == 0){
#   ggplot(data=outdt,aes(x=TrueCtrlMean,y=Prop,color=Method)) + geom_line(size = 1.5) + 
#     scale_x_continuous(breaks = muvec_c, name = "True Control Mean") + 
#     scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.01), name = "Probability of Success") + 
#     geom_vline(xintercept = -49.9, linetype=2, size=1) + 
#     geom_hline(yintercept = 0.025, linetype=2, size=1) + 
#     theme(axis.title = element_text(face="bold",size=15),
#           axis.text = element_text(size=12),
#           legend.title=element_text(size=15,face="bold"),
#           legend.text=element_text(size=12)
#     )
# }
