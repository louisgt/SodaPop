arg = commandArgs(trailingOnly=TRUE)

#### PARSE CURRENT SIMULATION DIRECTORY
src = arg[1]
dt=as.numeric(arg[2])
dir = getwd()

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(ggplot2, reshape2)

print("Importing time series data...")

#### IMPORT TIME SERIES DATA
generations = as.data.frame(read.csv(paste(dir,src,"ALL_generations.txt",sep=""), header = FALSE, sep =" "))

print("Extracting rows to keep...")
rows_to_keep = generations$V3 != 0
generations = generations[rows_to_keep,]

print("Importing average fitness trajectory...")
avg_fitness <- read.csv(paste(dir,src,"avg_fitness.txt",sep=""), header = FALSE)
avg_fitness$gen = as.numeric(row.names(avg_fitness))*dt-dt

fixgen = (ncol(generations) - 2)*dt
colnames(generations)[2:(ncol(generations))]=seq(0,fixgen, dt)
print("Melting dataframe...")
FLUX = melt(generations, id='V1')

print("Saving plot a to file...")

step = dt*10
if(fixgen/step > 12){
	step = dt*50
}

a = ggplot(avg_fitness, aes(x=gen,y=V1)) + geom_line() + theme_bw() + labs(x = "Generations",y="Average fitness")
#a = ggplot(avg_fitness, aes(x=gen,y=V1)) + stat_smooth(linetype="dashed",color="black",size=0.6,span=0.1) + theme_bw() + labs(x = "Generations",y="Average fitness")
ggsave("fitness.png", plot=last_plot(), path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

print("Saving plot b to file...")
b = ggplot(FLUX, aes(x=factor(variable),y=value,group=V1,colour=V1)) + geom_area(aes(fill=V1),alpha=0.5) + theme_bw() + scale_color_discrete(guide=FALSE) + scale_fill_discrete(guide=FALSE) + scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) +
  labs(x = "Generations",y="Count")
b$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave("clonal_structure.png", plot=b, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

print("Saving plot c to file...")

c = ggplot(FLUX, aes(x=factor(variable),y=value,group=V1,colour=V1)) + geom_line() + theme_bw() +  scale_color_discrete(guide=FALSE) + scale_x_discrete(limits=0:fixgen, breaks = seq(0,fixgen,step)) + 
labs(x = "Generations",y="Count")
c$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave("clonal_trajectories.png", plot=c, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=300)

#ggsave("clonal_trajectories.eps", plot=c, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=600)
#ggsave("clonal_structure.eps", plot=b, device=cairo_ps, fallback_resolution = 300, path = paste(dir,src,"graph/",sep=""), width = 11, height = 8.5, dpi=600)
#ggsave("combo.eps", arrangeGrob(c, b,ncol=2), path = paste(dir,src,"graph/",sep=""), width = 12, height = 5, dpi=600)