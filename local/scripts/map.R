#### packages ####
rm(list=ls())
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(devtools)
library(BEDASSLE)
library(ggplot2)
library(vcfR)
library(fossil)
library(viridis)
library(patchwork)
library(sp)
library(raster)
library(ecodist)
library(radiator)
library(wesanderson)
library(gridExtra)
#devtools::install_github("dkahle/ggmap")
library(ggmap)
#register_google(key="AIzaSyAz8Urlhyb4VqYuiy_dBIv-ietj7eY4YBo", write = TRUE)
options(scipen = 999)
library(RgoogleMaps)
library(topoDistance)
library(elevatr)
library(googleway)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
require(mapdata)
library(ggspatial)
library(ggsn)
library(cowplot)

##### functions ##### 
#### 0. prep data ####
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/7_Assocation_studies/9_ANGSD/12.era.mel.altitude/")
bams <- read.csv("../../11_all.bams/all.bams.list.info.csv")
bams <- subset(bams, bams$to.use.pbs=="yes")
nrow(subset(bams, species=="erato"))
nrow(subset(bams, species=="melpomene"))

######### 1. zoom in maps 4 regions, both sp #########
bams$type.an.pop <- paste(substr(bams$species,0,3),".",
                          tolower(substr(bams$country,0,2)),".",
                          tolower(substr(bams$type,str_length(bams$type), str_length(bams$type))),".",
                          bams$pop.short, sep = "")
# remove species from group_bya (and ,type.an.pop) , just by pop.short
pops <- dplyr::summarise(dplyr::group_by(bams,  pop.short) ,
                               n=n(),
                         altitude=mean(altitude),
                         latitude=mean(latitude),
                         longitude=mean(longitude),
                         type=unique(substr(type,2,4)),
                         country.side=unique(substr(type.an.pop,5,8)),
                               depth.mean=mean(mean.depth.wg.withzeros),
                               depth.sd=sd(mean.depth.wg.withzeros),
                               depth.max=max(mean.depth.wg.withzeros),
                               depth.min=min(mean.depth.wg.withzeros)); pops

# pops <- dplyr::summarise(dplyr::group_by(bams, species pop.short) ,
#                          n=n(),
#                          altitude=mean(altitude),
#                          latitude=mean(latitude),
#                          longitude=mean(longitude),
#                          type=unique(substr(type,2,4)),
#                          country.side=unique(substr(type.an.pop,5,8)),
#                          depth.mean=mean(mean.depth.wg.withzeros),
#                          depth.sd=sd(mean.depth.wg.withzeros),
#                          depth.max=max(mean.depth.wg.withzeros),
#                          depth.min=min(mean.depth.wg.withzeros)); pops


# prep pops, do not subsample as in bedassle, as for pbs we use all, &n>4
co.e <- subset(pops,country.side=="co.e"); co.e
co.w <- subset(pops,country.side=="co.w"); co.w
ec.e <- subset(pops,country.side=="ec.e"); ec.e
ec.w <- subset(pops,country.side=="ec.w"); ec.w

# co.e <- subset(pops,country.side=="co.e"&species=="erato"); co.e
# co.w <- subset(pops,country.side=="co.w"&species=="erato"); co.w
# ec.e <- subset(pops,country.side=="ec.e"&species=="erato"); ec.e
# ec.w <- subset(pops,country.side=="ec.w"&species=="erato"); ec.w


# prep transparent colors
adjustcolor( "blue", alpha.f = 0.95)
adjustcolor( "green", alpha.f = 0.7)
adjustcolor( "darkgreen", alpha.f = 0.5)

## co.e #### 
mean(co.e[c(co.e$type=="hig"),]$latitude);mean(co.e[c(co.e$type=="hig"),]$longitude)

co.e.map <- get_googlemap(center=c(lon=-76, lat=1.41), zoom = 9, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
co.e.x <- c(-76.8, -75.4); co.e.y <- c(0.57, 1.97)
gc(); co.e.map.gg <- ggmap(co.e.map) +
  #theme_bw() +
  geom_jitter(data = co.e,
              aes(x = longitude, y = latitude, shape = type, fill=altitude),
              size=8, alpha=.8, width = .03, height = .03, stroke=1)+
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  #scale_fill_gradientn(colors = c( "#00FF0080", "#00640080","#0000FFCC"))+
  scale_colour_manual(values=c("black","black"))+
  scale_x_continuous(limits = c(-76.8, -75.4, expand = c(0, 0))) +
  scale_y_continuous(limits = c(0.57, 1.97), expand = c(0, 0)) +
  # this will add a scalebar of 50km
  scalebar(dist = 25, dist_unit = "km",location="bottomleft",
           st.bottom = TRUE, st.color = "black", height=0.01, st.dist = 0.3,
           x.min = co.e.x[1], x.max = co.e.x[2], size=1,box.fill=c("black", "black"),
           y.min = co.e.y[1]+0.04, y.max = co.e.y[2], 
           transform = TRUE, model = "WGS84") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  guides(size = guide_legend(order=2)); co.e.map.gg
ggsave("../../../22_pop.gen.paper/analyses/fig.1.map/cline.maps/co.e.png", scale = 0.7)   

## co.e leticia #### 
mean(co.e[c(co.e$type=="vlo"),]$latitude);mean(co.e[c(co.e$type=="vlo"),]$longitude)

co.e.2.map <- get_googlemap(center=c(lon=-70, lat=-3.95), zoom = 8, color="bw",  maptype = "terrain",
                          style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
co.e.2.x <- c(-70.55, -69.8); co.e.2.y <- c(-4.35, -3.6)
gc(); co.e.map2.gg <- ggmap(co.e.2.map) +
  #theme_bw() +
  geom_jitter(data = co.e,
              aes(x = longitude, y = latitude, shape = type, fill=altitude),
              size=8, alpha=.8, width = .03, height = .03, stroke=1)+
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  #scale_fill_gradientn(colors = c( "#00FF0080", "#00640080","#0000FFCC"))+
  scale_colour_manual(values=c("black","black")) +
  scale_x_continuous(limits = c(-70.55, -69.8, expand = c(0, 0))) +
  scale_y_continuous(limits = c(-4.35, -3.6), expand = c(0, 0)) +
  # this will add a scalebar of 50km
  scalebar(dist = 25, dist_unit = "km",location="bottomleft",
           st.bottom = TRUE, st.color = "black", height=0.01, st.dist = 0.3,
           x.min = co.e.2.x[1], x.max = co.e.2.x[2], size=1,box.fill=c("black", "black"),
           y.min = co.e.2.y[1]+0.04, y.max = co.e.2.y[2], 
           transform = TRUE, model = "WGS84") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  guides(size = guide_legend(order=2)); co.e.map2.gg
ggsave("../../../22_pop.gen.paper/analyses/fig.1.map/cline.maps/co.e.2.png", scale = 0.3)    


## co.w ####  
mean(co.w$latitude);mean(co.w$longitude)
devtools::install_github("clauswilke/scales")
co.w.map <- get_googlemap(center=c(lon=-76.15, lat=4.9), zoom = 7, color="bw",  maptype = "terrain",
                          style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
co.w.x <- c(-78.7, -75.35); co.w.y <- c(3.2, 6.55)
gc(); co.w.map.gg <- ggmap(co.w.map) +
  geom_jitter(data = co.w,
              aes(x = longitude, y = latitude, shape = type, fill=altitude),
              size=8, width = .1, height = .05, stroke=1)+
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  scale_colour_manual(values=c("black", "black", "black"))+
 # scale_fill_gradientn(colors = c( "#00FF0080", "#00640080","#0000FFCC"))+
  scale_x_continuous(limits = c(-78.7, -75.35, expand = c(0, 0))) +
  scale_y_continuous(limits = c(3.2, 6.55), expand = c(0, 0)) +
  # this will add a scalebar of 50km
  scalebar(dist = 25, dist_unit = "km",location="bottomleft",
           st.bottom = TRUE, st.color = "black", height=0.01, st.dist = 0.03,
           x.min = co.w.x[1], x.max = co.w.x[2], size=1,box.fill=c("black", "black"),
           y.min = co.w.y[1]+0.04, y.max = co.w.y[2], 
           transform = TRUE, model = "WGS84") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  guides(size = guide_legend(order=2)); co.w.map.gg
ggsave("../../../22_pop.gen.paper/analyses/fig.1.map/cline.maps/co.w.png", scale = 0.7)


## ec.e ####  
mean(ec.e$latitude);mean(ec.e$longitude)
ec.e.map <- get_googlemap(center=c(lon=-77.1, lat=-0.56), zoom = 9, color="bw",  maptype = "terrain",
                          style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

ec.e.x <- c(-77.9, -76.3); ec.e.y <- c(-1.3, 0.3)
gc(); ec.e.map.gg <- ggmap(ec.e.map) +
  #theme_bw() +
  geom_jitter(data = pops,
              aes(x = longitude, y = latitude, shape = type, fill=altitude),
              size=8, alpha=.8, width = .02, height = .02, stroke=1)+
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  #scale_fill_gradientn(colors = c( "#00FF0080", "#00640080","#0000FFCC"))+
  scale_colour_manual(values=c("black", "black", "black"))+
  scale_x_continuous(limits = c(-77.9, -76.3, expand = c(0, 0))) +
  scale_y_continuous(limits = c(-1.3, 0.3), expand = c(0, 0)) +
  # this will add a scalebar of 50km
  scalebar(dist = 25, dist_unit = "km",location="bottomleft",
           st.bottom = TRUE, st.color = "black", height=0.01, st.dist = 0.03,
           x.min = ec.e.x[1], x.max = ec.e.x[2], size=1,box.fill=c("black", "black"),
           y.min = ec.e.y[1]+0.04, y.max = ec.e.y[2], 
           transform = TRUE, model = "WGS84") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  guides(size = guide_legend(order=2)); ec.e.map.gg
ggsave("../../../22_pop.gen.paper/analyses/fig.1.map/cline.maps/ec.e.png", scale = 0.7)

## ec.w
mean(ec.w$latitude);mean(ec.w$longitude)
devtools::install_github("clauswilke/scales")
ec.w.map <- get_googlemap(center=c(lon=-79.3, lat=0.29), zoom = 9, color="bw",  maptype = "terrain",
                          style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
ec.w.x <- c(-80.05, -78.7); ec.w.y <- c(-0.3, 1.05)
gc(); ec.w.map.gg <- ggmap(ec.w.map) +
  #theme_bw() +
  geom_jitter(data = ec.w,
              aes(x = longitude, y = latitude, shape = type, fill=altitude),
              size=8, width = .02, height = .02, stroke=1)+
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "#00FF00B3", high = "#0000FFF2")+
  scale_colour_manual(values=c("black", "black", "black")) +
  # scale_fill_gradientn(colors = c( "#00FF0080", "#00640080","#0000FFCC"))+
  scale_x_continuous(limits = c(-80.05, -78.7, expand = c(0, 0))) +
  scale_y_continuous(limits = c(-0.3, 1.05), expand = c(0, 0)) +
  # this will add a scalebar of 50km
  scalebar(dist = 25, dist_unit = "km",location="bottomleft",
           st.bottom = TRUE, st.color = "black", height=0.01, st.dist = 0.06,
           x.min = ec.w.x[1], x.max = ec.w.x[2], size=1,box.fill=c("black", "black"),
           y.min = ec.w.y[1]+0.04, y.max = ec.w.y[2], 
           transform = TRUE, model = "WGS84") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  guides(size = guide_legend(order=2)); ec.w.map.gg
ggsave("../../../22_pop.gen.paper/analyses/fig.1.map/cline.maps/ec.w.png", scale = 0.7)


plot_grid(co.w.map.gg, co.e.map.gg, ec.w.map.gg, ec.e.map.gg, ncol=2)

## all SA boxes ## 
extend(ec.e.map)
GetMap.bbox(ec.e.map)
map <- get_googlemap(center=c(lon=-74.54, lat=-0.2602158), zoom = 5, color="bw",  maptype = "terrain",
                          style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
dev.off()
gc(); p1 <- ggmap(map) +
  geom_jitter(data = subset(pops, n>2),
              aes(x = longitude, y = latitude, shape = type),
              size=.2, alpha=.8, width = .02, height = .02, stroke=1)+
  scale_x_continuous(limits = c(-88, -65), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-9, 13), expand = c(0, 0)) +
  scale_shape_manual(values = c(24,22,21)) +
  scale_colour_manual(values=c("black", "black", "white"))+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_rect(mapping = aes(xmin=ec.e.x[1], xmax=ec.e.x[2], ymin=ec.e.y[1], ymax=ec.e.y[2]),fill=NA,color="black", alpha=.5)+
  geom_rect(mapping = aes(xmin=ec.w.x[1], xmax=ec.w.x[2], ymin=ec.w.y[1], ymax=ec.w.y[2]),fill=NA,color="black", alpha=.5)+
  geom_rect(mapping = aes(xmin=co.w.x[1], xmax=co.w.x[2], ymin=co.w.y[1], ymax=co.w.y[2]),fill=NA,color="black", alpha=.5)+
  geom_rect(mapping = aes(xmin=co.e.x[1], xmax=co.e.x[2], ymin=co.e.y[1], ymax=co.e.y[2]),fill=NA,color="black", alpha=.5)+
  geom_rect(mapping = aes(xmin=co.e.2.x[1], xmax=co.e.2.x[2], ymin=co.e.2.y[1], ymax=co.e.2.y[2]),fill=NA,color="black", alpha=.5)+
  guides(size = guide_legend(order=2)); p1
ggsave("../../../22_pop.gen.paper/analyses/fig.1.map/cline.maps/SA.png")


#### alt types plot ####

# mean distance between high and low pops
pop.type.ec.e <- dplyr::summarise(dplyr::group_by(ec.e, type, country.side),
          altitude=mean(altitude),
          latitude=mean(latitude),
          longitude=mean(longitude))
era.ec.e.dist <- earth.dist(as.data.frame(pop.type.ec.e[c("longitude", "latitude")]), dist = TRUE); era.ec.e.dist
era.ec.e.dist<- as.matrix(era.ec.e.dist); colnames(era.ec.e.dist) <- c(as.character(pop.type.ec.e$type)); rownames(era.ec.e.dist) <- c(as.character(pop.type.ec.e$type))

pop.type.ec.w <- dplyr::summarise(dplyr::group_by(ec.w, type, country.side),
                                  altitude=mean(altitude),
                                  latitude=mean(latitude),
                                  longitude=mean(longitude)); pop.type.ec.w
era.ec.w.dist <- earth.dist(as.data.frame(pop.type.ec.w[c("longitude", "latitude")]), dist = TRUE); era.ec.w.dist
era.ec.w.dist<- as.matrix(era.ec.w.dist); colnames(era.ec.w.dist) <- c(as.character(pop.type.ec.w$type)); rownames(era.ec.w.dist) <- c(as.character(pop.type.ec.w$type))

pop.type.co.e <- dplyr::summarise(dplyr::group_by(co.e, type, country.side),
                                  altitude=mean(altitude),
                                  latitude=mean(latitude),
                                  longitude=mean(longitude))
era.co.e.dist <- earth.dist(as.data.frame(pop.type.co.e[c("longitude", "latitude")]), dist = TRUE); era.co.e.dist
era.co.e.dist<- as.matrix(era.co.e.dist); colnames(era.co.e.dist) <- c(as.character(pop.type.co.e$type)); rownames(era.co.e.dist) <- c(as.character(pop.type.co.e$type))

pop.type.co.w <- dplyr::summarise(dplyr::group_by(co.w, type, country.side),
                                  altitude=mean(altitude),
                                  latitude=mean(latitude),
                                  longitude=mean(longitude)); pop.type.co.w
era.co.w.dist <- earth.dist(as.data.frame(pop.type.co.w[c("longitude", "latitude")]), dist = TRUE); era.co.w.dist
era.co.w.dist<- as.matrix(era.co.w.dist); colnames(era.co.w.dist) <- c(as.character(pop.type.co.w$type)); rownames(era.co.w.dist) <- c(as.character(pop.type.co.w$type))

mean.dist.all <- apply(data.frame(era.ec.e.dist[,1],era.ec.w.dist[,1],era.co.e.dist[,1],era.co.w.dist[,1]), 1, mean) 
pop.type <- dplyr::summarise(dplyr::group_by(pops, type),
                                  altitude=mean(altitude)); pop.type 
pop.type$dist.to.hig <- c(mean.dist.all); pop.type

ggplot(aes(x=dist.to.hig, y=altitude, shape=type, fill=altitude), data=pop.type)+
  #geom_line(aes(x=dist.to.hig, y=altitude), inherit.aes = FALSE, size=8, colour="grey24")+
  geom_point( size=30, stroke=4)+
  ylim(-50,1350)+
  xlim(-10,350)+
  #geom_text(aes(label=c("High", "Low", "Low distant")),size=10, fontface = "bold", nudge_x = c(15,15,-30),hjust = 0, nudge_y = c(20, 85, 80))+
  scale_shape_manual(values = c(24,22,21), guide=FALSE) +
  scale_fill_gradient( low = "green", high = "blue")+
  theme_classic()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c("right"), axis.text = element_text(size=80, face="bold", colour = "black"),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_line(size = 3, colour = "black"),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave("../../../22_pop.gen.paper/analyses/fig.1/cline.maps/cline.points.no.line.png",  bg = "transparent", scale=1, width = 16, height = 8)

## era.e only for fig4
ggplot(aes(x=dist.to.hig, y=altitude, shape=type, fill=altitude), data=pop.type)+
  geom_line(aes(x=dist.to.hig, y=altitude), inherit.aes = FALSE, size=1, colour="grey24")+
  geom_point( size=8, stroke=1)+
  #geom_text(aes(label=c("High", "Low", "Low distant")),size=10, fontface = "bold", nudge_x = c(15,15,-30),hjust = 0, nudge_y = c(20, 85, 80))+
  scale_shape_manual(values = c(24,22,21), guide="none") +
  scale_fill_gradient( low = "green", high = "blue")+
  scale_x_continuous(expand = c(0.2,0.2))+
  scale_y_continuous(expand = c(0.2,0.5), limits = c(0,1300))+
  theme_classic()+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c("right"), axis.text = element_text(size=14,  colour = "black"),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_line(size = 1, colour = "black"),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave("../../../22_pop.gen.paper/figures/fig4.inv.local.pca/cline.points.png",  bg = "transparent", scale=1, width = 4, height = 2)



