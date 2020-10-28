library(ggplot2)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)

theme_Publication = function(base_size=20,base_family="Helvetica"){
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size,base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size=rel(1.2),hjust=0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black",size=1.5),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

fixgen = 250
dt = 10
step = 50

lower_bound = 0.0001
upper_bound = 0.0005

m = as.matrix(generations[,-1])
mat = as.data.frame(sweep(m,2,colSums(m,na.rm = TRUE),`/`))
colnames(mat)[1:(ncol(mat))]=seq(0,fixgen, dt)
mat$ID = row.names(mat)

mat[,"max"] = apply(mat[,-ncol(mat)],1, max)
sub_gen = mat[mat$max >= lower_bound,]

df = melt(sub_gen,id.vars = c('ID','max'))
grouped_df = df %>% group_by(max) %>% filter(max > upper_bound) %>% ungroup()
g = ggplot() + geom_line(aes(x=factor(variable),y=value,group=ID),data=df,colour="#CCCCCC",alpha=0.7) + geom_line(aes(x=factor(variable),y=value,group=ID,colour=ID),data=grouped_df,size=0.9) + scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) + theme_Publication() + scale_y_log10(limits=c(1e-7,1e0)) + scale_color_viridis_d(guide=FALSE,option = "D",direction = -1) + xlab("Generations") + ylab("Barcode count")
ggsave(g,filename = "Documents/log_geom_line.png",dpi=300,width = 10,height = 8)

tf = df[order(df$max),]
grouped_tf = tf %>% group_by(max) %>% filter(max > upper_bound) %>% ungroup()

n <- length(grouped_tf$ID)
cols <- hue_pal(h = c(0, 360) + 15,
c = 125, l = 72,
h.start = 0, direction = 1)(n)[order(sample(1:n, n))]

g = ggplot() + geom_area(aes(x=factor(variable),y=value,group=ID),data=tf,fill="#CCCCCC") + geom_area(aes(x=factor(variable),y=value,group=ID,fill=ID),data=grouped_tf) + scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) + theme_Publication() + scale_fill_manual(values = cols,guide=FALSE) + xlab("Generations") + ylab("Barcode fraction") + ylim(c(0,1))
ggsave(g,filename = "Documents/PROJECTS/SodaPop/out/TEST_DEV/graph/area.png",dpi=300,width = 10,height = 8)

