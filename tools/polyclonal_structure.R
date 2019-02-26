arg = commandArgs(trailingOnly=TRUE)

#### PARSE CURRENT SIMULATION DIRECTORY
src = arg[1]
dt=as.numeric(arg[2])
dir = getwd()

#### LOAD UTILITIES (INSTALL IF MISSING)
if(!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(ggplot2,reshape2,grid,ggthemes,viridis)
library(ggplot2)
library(reshape2)
library(viridis)

message("R: Importing time series data...")

#### IMPORT TIME SERIES DATA
generations = as.data.frame(read.csv(paste(dir,src,"ALL_generations.txt",sep=""), header = FALSE, sep =" "))

message("R: Extracting rows to keep...")
rows_to_keep = generations$V3 != 0
if(nrow(generations) > 100000){
	rows_to_keep = generations$V7 != 0
}
generations_filtered = generations[rows_to_keep,]

message("R: Importing average fitness trajectory...")
avg_fitness <- read.csv(paste(dir,src,"avg_fitness.txt",sep=""), header = FALSE)
avg_fitness$gen = as.numeric(row.names(avg_fitness))*dt-dt

message("R: Importing average stability trajectory...")
stabil_0 <- read.csv(paste(dir,src,"gene/stabil_0.txt",sep=""), header = FALSE)
#stabil_1 <- read.csv(paste(dir,src,"gene/stabil_1.txt",sep=""), header = FALSE)
#stabil_2 <- read.csv(paste(dir,src,"gene/stabil_2.txt",sep=""), header = FALSE)
#stabil_3 <- read.csv(paste(dir,src,"gene/stabil_3.txt",sep=""), header = FALSE)
#stabil_4 <- read.csv(paste(dir,src,"gene/stabil_4.txt",sep=""), header = FALSE)
#stabil_5 <- read.csv(paste(dir,src,"gene/stabil_5.txt",sep=""), header = FALSE)
#stabil_6 <- read.csv(paste(dir,src,"gene/stabil_6.txt",sep=""), header = FALSE)
#stabil_7 <- read.csv(paste(dir,src,"gene/stabil_7.txt",sep=""), header = FALSE)
#stabil_8 <- read.csv(paste(dir,src,"gene/stabil_8.txt",sep=""), header = FALSE)
#stabil_9 <- read.csv(paste(dir,src,"gene/stabil_9.txt",sep=""), header = FALSE)


fixgen = (ncol(generations_filtered) - 2)*dt
colnames(generations_filtered)[2:(ncol(generations_filtered))]=seq(0,fixgen, dt)
message("R: Melting dataframe...")
FLUX = melt(generations_filtered, id='V1')

fixgen_all = (ncol(generations) - 2)*dt
colnames(generations)[2:(ncol(generations))]=seq(0,fixgen_all, dt)
message("R: Melting dataframe...")
FLUX_ALL = melt(generations, id='V1')

#df = cbind(stabil_0,stabil_1,stabil_2,stabil_3,stabil_4,stabil_5,stabil_6,stabil_7,stabil_8,stabil_9)
#df = cbind(stabil_0,stabil_1,stabil_2,stabil_3,stabil_4,stabil_5,stabil_6,stabil_7,stabil_8)
#tf = as.data.frame(t(df))
#colnames(tf)[1:(ncol(tf))]=seq(0,fixgen, dt)
#tf$Gene = c(0,1,2,3,4,5,6,7,8,9)
#tf$Gene = c(0,1,2,3,4,5,6,7,8)
#tf$Gene = c(0)
#GENE = melt(tf, id='Gene')

message("R: Saving plot a to file...")

step = dt*10
if(fixgen/step > 12){
	step = dt*50
}

#"publication ready" plot theme
theme_Publication <- function(base_size=18, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size=1.5),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

a = ggplot(avg_fitness, aes(x=gen,y=V1)) + geom_line() + theme_Publication() + labs(x = "Generations",y="Average fitness")
ggsave("fitness.png", plot=last_plot(), path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

message("R: Saving plot b to file...")
b = ggplot(FLUX, aes(x=factor(variable),y=value,group=V1,colour=V1)) + geom_area(aes(fill=V1),alpha=0.5) + theme_Publication() + scale_color_discrete(guide=FALSE) + scale_fill_discrete(guide=FALSE) + scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) +
  labs(x = "Generations",y="Count")
b$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave("clonal_structure.png", plot=b, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

message("R: Saving plot c to file...")

c = ggplot(FLUX, aes(x=factor(variable),y=value,group=V1,colour=V1)) + geom_line() + theme_Publication() +  scale_color_discrete(guide=FALSE) + scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) + 
labs(x = "Generations",y="Count")
c$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave("clonal_trajectories.png", plot=c, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

message("R: Saving plot d to file...")

d = ggplot(FLUX_ALL, aes(x=factor(variable),y=log10(value),group=V1,colour=V1)) + geom_line() + theme_Publication() +  scale_color_discrete(guide=FALSE) + scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) + 
labs(x = "Generations",y="Log10(Count)")
d$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave("log_clonal_trajectories.png", plot=d, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

#e = ggplot(GENE, aes(x=factor(variable),y=value,group=factor(Gene),colour=factor(Gene))) + 
#    geom_line(size = 0.75) + scale_color_viridis_d(name = "Gene ID",option = "C") + theme_Publication() + 
#    scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) + labs(x = "Generations",y="Stability (∆G, kcal/mol)") + 
#    guides(color = guide_legend(order = 2, override.aes = list(shape = 15, size = 9)), shape = guide_legend(order = 1))
#e$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
#ggsave("stability_trajectories.png", plot=e, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

#ggsave("log_clonal_trajectories.eps", plot=d, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=600)
#ggsave("clonal_structure.eps", plot=b, device=cairo_ps, fallback_resolution = 300, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=600)
#ggsave("combo.eps", arrangeGrob(b, d,ncol=2), path = paste(dir,src,"graph/",sep=""), width = 12, height = 5, dpi=600)