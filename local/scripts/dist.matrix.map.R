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
#register_google(key="XXX", write = TRUE) # replace with your own key
options(scipen = 999)
library(RgoogleMaps)
library(topoDistance)
library(elevatr)

##### functions ##### 
# set palette
scale_fill_wes <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(wes_palette("BottleRocket2", 15, type = "continuous"), levels(era.ec.e$type)), 
    ...)}

scale_color_wes <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(wes_palette("BottleRocket2", 15, type = "continuous"), levels(era.ec.e$type)), 
    ...)}

#### 0. prep data ####
setwd("")
bam.list <- read.csv("02.info/all.bams.list.info.csv")
names(bam.list)

era.ec.e <- subset(bam.list, species=="erato"&to.use.pbs=="yes"&(type=="EhighE"|type=="ElowE"|type=="EvlowE"))
era.ec.e$altitude <- as.numeric(as.character(era.ec.e$altitude))
era.ec.w <- subset(bam.list, species=="erato"&to.use.pbs=="yes"&(type=="EhighW"|type=="ElowW"|type=="EvlowW"))
era.co.e <- subset(bam.list, species=="erato"&to.use.pbs=="yes"&(type=="ChighE"|type=="ClowE"|type=="CvlowE"))
era.co.w <- subset(bam.list, species=="erato"&to.use.pbs=="yes"&(type=="ChighW"|type=="ClowW"|type=="CvlowW"))
era.e <- subset(bam.list, species=="erato"&to.use.pbs=="yes"&(type=="EhighE"|type=="ElowE"|type=="EvlowE"|type=="ChighE"|type=="ClowE"|type=="CvlowE" ))
era.w <- subset(bam.list, species=="erato"&to.use.pbs=="yes"&(type=="EhighW"|type=="ElowW"|type=="EvlowW"|type=="ChighW"|type=="ClowW"|type=="CvlowW" ))
era <- subset(bam.list, species=="erato"&to.use.pbs=="yes")

### make pop.name.full dists ####
# instead of indiv
##### geo.dist ####
era.ec.e.loc <- era.ec.e[,c("pop.name", "latitude", "longitude")]; era.ec.e.loc
era.ec.w.loc <- era.ec.w[,c("pop.name", "latitude", "longitude")]; era.ec.w.loc
era.e.loc <- era.e[,c("pop.name", "latitude", "longitude")]; era.e.loc
era.w.loc <- era.w[,c("pop.name", "latitude", "longitude")]; era.w.loc
era.loc <- era[,c("pop.name", "latitude", "longitude")]; era.loc

era.ec.e.loc <- dplyr::summarise(dplyr::group_by(era.ec.e.loc, pop.name),
                          latitude=mean(latitude),
                          longitude=mean(longitude))
era.ec.w.loc <- dplyr::summarise(dplyr::group_by(era.ec.w.loc, pop.name),
                          latitude=mean(latitude),
                          longitude=mean(longitude))

era.e.loc <- dplyr::summarise(dplyr::group_by(era.e.loc, pop.name),
                                 latitude=mean(latitude),
                                 longitude=mean(longitude))
era.w.loc <- dplyr::summarise(dplyr::group_by(era.w.loc, pop.name),
                              latitude=mean(latitude),
                              longitude=mean(longitude))
era.loc <- dplyr::summarise(dplyr::group_by(era.loc, pop.name),
                              latitude=mean(latitude),
                              longitude=mean(longitude))

era.ec.e.dist <- earth.dist(as.data.frame(era.ec.e.loc[c("longitude", "latitude")]), dist = TRUE); era.ec.e.dist
era.ec.w.dist <- earth.dist(as.data.frame(era.ec.w.loc[c("longitude", "latitude")]), dist = TRUE); era.ec.w.dist
era.e.dist <- earth.dist(as.data.frame(era.e.loc[c("longitude", "latitude")]), dist = TRUE); era.e.dist
era.w.dist <- earth.dist(as.data.frame(era.w.loc[c("longitude", "latitude")]), dist = TRUE); era.w.dist
era.dist <- earth.dist(as.data.frame(era.loc[c("longitude", "latitude")]), dist = TRUE); era.dist

era.ec.e.dist<- as.matrix(era.ec.e.dist); colnames(era.ec.e.dist) <- c(as.character(era.ec.e.loc$pop.name)); rownames(era.ec.e.dist) <- c(as.character(era.ec.e.loc$pop.name))
era.ec.w.dist<- as.matrix(era.ec.w.dist); colnames(era.ec.w.dist) <- c(as.character(era.ec.w.loc$pop.name)); rownames(era.ec.w.dist) <- c(as.character(era.ec.w.loc$pop.name))
era.e.dist<- as.matrix(era.e.dist); colnames(era.e.dist) <- c(as.character(era.e.loc$pop.name)); rownames(era.e.dist) <- c(as.character(era.e.loc$pop.name))
era.w.dist<- as.matrix(era.w.dist); colnames(era.w.dist) <- c(as.character(era.w.loc$pop.name)); rownames(era.w.dist) <- c(as.character(era.w.loc$pop.name))
era.dist<- as.matrix(era.dist); colnames(era.dist) <- c(as.character(era.loc$pop.name)); rownames(era.dist) <- c(as.character(era.loc$pop.name))

##### alt dist ######
era.ec.e.alt <- era.ec.e[,c("pop.name", "altitude")]; era.ec.e.alt
era.ec.w.alt <- era.ec.w[,c("pop.name", "altitude")]; era.ec.w.alt
era.e.alt <- era.e[,c("pop.name", "altitude")]; era.e.alt
era.w.alt <- era.w[,c("pop.name", "altitude")]; era.w.alt
era.alt <- era[,c("pop.name", "altitude")]; era.alt

era.ec.e.alt <- dplyr::summarise(dplyr::group_by(era.ec.e.alt, pop.name),
                          altitude=mean(altitude))
era.ec.w.alt <- dplyr::summarise(dplyr::group_by(era.ec.w.alt, pop.name),
                          altitude=mean(altitude))
era.e.alt <- dplyr::summarise(dplyr::group_by(era.e.alt, pop.name),
                                 altitude=mean(altitude))
era.w.alt <- dplyr::summarise(dplyr::group_by(era.w.alt, pop.name),
                              altitude=mean(altitude))
era.alt <- dplyr::summarise(dplyr::group_by(era.alt, pop.name),
                              altitude=mean(altitude))

era.ec.e.alt.dist <- as.matrix(dist(era.ec.e.alt$altitude)); colnames(era.ec.e.alt.dist) <- c(as.character(era.ec.e.alt$pop.name)); rownames(era.ec.e.alt.dist) <- c(as.character(era.ec.e.alt$pop.name))
era.ec.w.alt.dist <- as.matrix(dist(era.ec.w.alt$altitude)); colnames(era.ec.w.alt.dist) <- c(as.character(era.ec.w.alt$pop.name)); rownames(era.ec.w.alt.dist) <- c(as.character(era.ec.w.alt$pop.name))
era.e.alt.dist <- as.matrix(dist(era.e.alt$altitude)); colnames(era.e.alt.dist) <- c(as.character(era.e.alt$pop.name)); rownames(era.e.alt.dist) <- c(as.character(era.e.alt$pop.name))
era.w.alt.dist <- as.matrix(dist(era.w.alt$altitude)); colnames(era.w.alt.dist) <- c(as.character(era.w.alt$pop.name)); rownames(era.w.alt.dist) <- c(as.character(era.w.alt$pop.name))
era.alt.dist <- as.matrix(dist(era.alt$altitude)); colnames(era.alt.dist) <- c(as.character(era.alt$pop.name)); rownames(era.alt.dist) <- c(as.character(era.alt$pop.name))

## pop n
era.ec.e.n <- dplyr::summarise(dplyr::group_by(era.ec.e, pop.name),
                        n=n());era.ec.e.n
era.ec.w.n <- dplyr::summarise(dplyr::group_by(era.ec.w, pop.name),
                        n=n()); era.ec.w.n

era.e.n <- dplyr::summarise(dplyr::group_by(era.e, pop.name),
                               n=n());era.e.n

era.w.n <- dplyr::summarise(dplyr::group_by(era.w, pop.name),
                            n=n());era.w.n

era.n <- dplyr::summarise(dplyr::group_by(era.w, pop.name),
                            n=n());era.n
## save matrices
write.table(era.ec.e.dist, "local/output/dist.vars/era.ec.e.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.ec.w.dist, "local/output/dist.vars/era.ec.w.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.ec.e.alt.dist, "local/output/e.vars/era.ec.e.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.ec.w.alt.dist, "local/output/e.vars/era.ec.w.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.ec.e.n, "local/output/dist.vars/era.ec.e.dist.pop.n.txt", row.names = TRUE, col.names = TRUE)
write.table(era.ec.w.n, "local/output/dist.vars/era.ec.w.dist.pop.n.txt", row.names = TRUE, col.names = TRUE)
write.table(era.e.dist, "local/output/dist.vars/era.e.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.e.alt.dist, "local/output/e.vars/era.e.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.w.dist, "local/output/dist.vars/era.w.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.w.alt.dist, "local/output/e.vars/era.w.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.dist, "local/output/dist.vars/era.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(era.alt.dist, "local/output/e.vars/era.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
##### maps ####
mean(era.ec.e.loc$latitude);mean(era.ec.e.loc$longitude)
era.ec.e.map <- get_googlemap(center=c(lon=-77.2615, lat=-0.700245), zoom = 9, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")


gc(); p1 <- ggmap(era.ec.e.map) +
  #theme_bw() +
  geom_jitter(data = era.ec.e, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.ec.e.png")

## west 

mean(era.ec.w.loc$latitude);mean(era.ec.w.loc$longitude)
era.ec.w.map <- get_googlemap(center=c(lon=-78.96914, lat=0.21412), zoom = 9, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(era.ec.w.map) +
  #theme_bw() +
  geom_jitter(data = era.ec.w, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.ec.w.png")


mean(era.co.e$latitude);mean(era.co.e$longitude)
era.co.e.map <- get_googlemap(center=c(lon= -76.68982, lat=1.214120), zoom = 8, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(era.co.e.map) +
  #theme_bw() +
  geom_jitter(data = era.co.e, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.co.e.png")

## west 

mean(era.co.w$latitude);mean(era.co.w$longitude)
era.co.w.map <- get_googlemap(center=c(lon=-77.08683, lat=3.87885), zoom = 8, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(era.co.w.map) +
  #theme_bw() +
  geom_jitter(data = era.co.w, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.co.w.png")



### all east ##

mean(era.e.loc$latitude);mean(era.e.loc$longitude)
era.e.map <- get_googlemap(center=c(lon=-75.67008, lat=-0.4237762), zoom = 6, color="bw",  maptype = "terrain",
                           style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
era.e.map

gc(); p1 <- ggmap(era.e.map) +
  #theme_bw() +
  geom_jitter(data = era.e, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21,24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.e.png")

### all west ##

mean(era.w.loc$latitude);mean(era.w.loc$longitude)
era.w.map <- get_googlemap(center=c(lon=-78.49544, lat=2.2), zoom = 7, color="bw",  maptype = "terrain",
                           style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(era.w.map) +
  #theme_bw() +
  geom_jitter(data = era.w, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21,24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.w.png")


### all era ##

mean(era.loc$latitude);mean(era.loc$longitude)
era.map <- get_googlemap(center=c(lon=-74.5, lat=0.5), zoom = 6, color="bw",  maptype = "terrain",
                         style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

era$alt.type <- str_sub(era$type, 2, 4)
gc(); p1 <- ggmap(era.map) +
  #theme_bw() +
  geom_jitter(data = era, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.png")


#### by pop #####
names(era)
era.pop <-  summarise(group_by(era, pop.name, pop.short),
          n=n(),
          latitude=mean(latitude),
          longitude=mean(longitude),
          altitude=mean(altitude),
          alt.type=unique(str_sub(type, 2, 4)))

gc(); p1 <- ggmap(era.map) +
  #theme_bw() +
  geom_jitter(data = era.pop, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  xlab("Longitude") + ylab("Latitude") +
  geom_text(data=era.pop, aes(x=longitude, y=latitude,label=n),hjust=1.5, vjust=1.5, inherit.aes = FALSE)+
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.pops.png")

gc(); p1 <- ggmap(era.map) +
  #theme_bw() +
  geom_jitter(data = era.pop, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  xlab("Longitude") + ylab("Latitude") +
  geom_text(data=era.pop, aes(x=longitude, y=latitude,label=n),hjust=1.5, vjust=1.5, inherit.aes = FALSE)+
  scale_shape_manual(values = c(24,22,21))+
  ylim(-1.5,5.7)+xlim(-80.1,-75.3)+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.pops.zoom.png")

gc(); p1 <- ggmap(era.map) +
  #theme_bw() +
  geom_jitter(data = era.pop, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  xlab("Longitude") + ylab("Latitude") +
  geom_text(data=era.pop, aes(x=longitude, y=latitude,label=n),hjust=1.5, vjust=1.5, inherit.aes = FALSE)+
  scale_shape_manual(values = c(24,22,21))+
  ylim(-1.5,2.1)+xlim(-80.1,-75.3)+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/era.pops.zoom2.png")

############# mel ####
mel.ec.e <- subset(bam.list, species=="melpomene"&to.use.pbs=="yes"&(type=="EhighE"|type=="ElowE"|type=="EvlowE"))
mel.ec.e$altitude <- as.numeric(as.character(mel.ec.e$altitude))
mel.ec.w <- subset(bam.list, species=="melpomene"&to.use.pbs=="yes"&(type=="EhighW"|type=="ElowW"|type=="EvlowW"))
mel.co.e <- subset(bam.list, species=="melpomene"&to.use.pbs=="yes"&(type=="ChighE"|type=="ClowE"|type=="CvlowE"))
mel.co.w <- subset(bam.list, species=="melpomene"&to.use.pbs=="yes"&(type=="ChighW"|type=="ClowW"|type=="CvlowW"))
mel.e <- subset(bam.list, species=="melpomene"&to.use.pbs=="yes"&(type=="EhighE"|type=="ElowE"|type=="EvlowE"|type=="ChighE"|type=="ClowE"|type=="CvlowE" ))
mel.w <- subset(bam.list, species=="melpomene"&to.use.pbs=="yes"&(type=="EhighW"|type=="ElowW"|type=="EvlowW"|type=="ChighW"|type=="ClowW"|type=="CvlowW" ))
mel <- subset(bam.list, species=="melpomene"&to.use.pbs=="yes")

### make pop.name.full dists ####
# instead of indiv
##### geo.dist ####
mel.ec.e.loc <- mel.ec.e[,c("pop.name", "latitude", "longitude")]; mel.ec.e.loc
mel.ec.w.loc <- mel.ec.w[,c("pop.name", "latitude", "longitude")]; mel.ec.w.loc
mel.e.loc <- mel.e[,c("pop.name", "latitude", "longitude")]; mel.e.loc
mel.w.loc <- mel.w[,c("pop.name", "latitude", "longitude")]; mel.w.loc
mel.loc <- era[,c("pop.name", "latitude", "longitude")]; mel.loc

mel.ec.e.loc <- dplyr::summarise(dplyr::group_by(mel.ec.e.loc, pop.name),
                                 latitude=mean(latitude),
                                 longitude=mean(longitude))
mel.ec.w.loc <- dplyr::summarise(dplyr::group_by(mel.ec.w.loc, pop.name),
                                 latitude=mean(latitude),
                                 longitude=mean(longitude))

mel.e.loc <- dplyr::summarise(dplyr::group_by(mel.e.loc, pop.name),
                              latitude=mean(latitude),
                              longitude=mean(longitude))
mel.w.loc <- dplyr::summarise(dplyr::group_by(mel.w.loc, pop.name),
                              latitude=mean(latitude),
                              longitude=mean(longitude))
mel.loc <- dplyr::summarise(dplyr::group_by(mel.loc, pop.name),
                            latitude=mean(latitude),
                            longitude=mean(longitude))

mel.ec.e.dist <- earth.dist(as.data.frame(mel.ec.e.loc[c("longitude", "latitude")]), dist = TRUE); mel.ec.e.dist
mel.ec.w.dist <- earth.dist(as.data.frame(mel.ec.w.loc[c("longitude", "latitude")]), dist = TRUE); mel.ec.w.dist
mel.e.dist <- earth.dist(as.data.frame(mel.e.loc[c("longitude", "latitude")]), dist = TRUE); mel.e.dist
mel.w.dist <- earth.dist(as.data.frame(mel.w.loc[c("longitude", "latitude")]), dist = TRUE); mel.w.dist
mel.dist <- earth.dist(as.data.frame(mel.loc[c("longitude", "latitude")]), dist = TRUE); mel.dist

mel.ec.e.dist<- as.matrix(mel.ec.e.dist); colnames(mel.ec.e.dist) <- c(as.character(mel.ec.e.loc$pop.name)); rownames(mel.ec.e.dist) <- c(as.character(mel.ec.e.loc$pop.name))
mel.ec.w.dist<- as.matrix(mel.ec.w.dist); colnames(mel.ec.w.dist) <- c(as.character(mel.ec.w.loc$pop.name)); rownames(mel.ec.w.dist) <- c(as.character(mel.ec.w.loc$pop.name))
mel.e.dist<- as.matrix(mel.e.dist); colnames(mel.e.dist) <- c(as.character(mel.e.loc$pop.name)); rownames(mel.e.dist) <- c(as.character(mel.e.loc$pop.name))
mel.w.dist<- as.matrix(mel.w.dist); colnames(mel.w.dist) <- c(as.character(mel.w.loc$pop.name)); rownames(mel.w.dist) <- c(as.character(mel.w.loc$pop.name))
mel.dist<- as.matrix(mel.dist); colnames(mel.dist) <- c(as.character(mel.loc$pop.name)); rownames(mel.dist) <- c(as.character(mel.loc$pop.name))

##### alt dist ######
mel.ec.e.alt <- mel.ec.e[,c("pop.name", "altitude")]; mel.ec.e.alt
mel.ec.w.alt <- mel.ec.w[,c("pop.name", "altitude")]; mel.ec.w.alt
mel.e.alt <- mel.e[,c("pop.name", "altitude")]; mel.e.alt
mel.w.alt <- mel.w[,c("pop.name", "altitude")]; mel.w.alt
mel.alt <- era[,c("pop.name", "altitude")]; mel.alt

mel.ec.e.alt <- dplyr::summarise(dplyr::group_by(mel.ec.e.alt, pop.name),
                                 altitude=mean(altitude))
mel.ec.w.alt <- dplyr::summarise(dplyr::group_by(mel.ec.w.alt, pop.name),
                                 altitude=mean(altitude))
mel.e.alt <- dplyr::summarise(dplyr::group_by(mel.e.alt, pop.name),
                              altitude=mean(altitude))
mel.w.alt <- dplyr::summarise(dplyr::group_by(mel.w.alt, pop.name),
                              altitude=mean(altitude))
mel.alt <- dplyr::summarise(dplyr::group_by(mel.alt, pop.name),
                            altitude=mean(altitude))

mel.ec.e.alt.dist <- as.matrix(dist(mel.ec.e.alt$altitude)); colnames(mel.ec.e.alt.dist) <- c(as.character(mel.ec.e.alt$pop.name)); rownames(mel.ec.e.alt.dist) <- c(as.character(mel.ec.e.alt$pop.name))
mel.ec.w.alt.dist <- as.matrix(dist(mel.ec.w.alt$altitude)); colnames(mel.ec.w.alt.dist) <- c(as.character(mel.ec.w.alt$pop.name)); rownames(mel.ec.w.alt.dist) <- c(as.character(mel.ec.w.alt$pop.name))
mel.e.alt.dist <- as.matrix(dist(mel.e.alt$altitude)); colnames(mel.e.alt.dist) <- c(as.character(mel.e.alt$pop.name)); rownames(mel.e.alt.dist) <- c(as.character(mel.e.alt$pop.name))
mel.w.alt.dist <- as.matrix(dist(mel.w.alt$altitude)); colnames(mel.w.alt.dist) <- c(as.character(mel.w.alt$pop.name)); rownames(mel.w.alt.dist) <- c(as.character(mel.w.alt$pop.name))
mel.alt.dist <- as.matrix(dist(mel.alt$altitude)); colnames(mel.alt.dist) <- c(as.character(mel.alt$pop.name)); rownames(mel.alt.dist) <- c(as.character(mel.alt$pop.name))

## pop n
mel.ec.e.n <- dplyr::summarise(dplyr::group_by(mel.ec.e, pop.name),
                               n=n());mel.ec.e.n
mel.ec.w.n <- dplyr::summarise(dplyr::group_by(mel.ec.w, pop.name),
                               n=n()); mel.ec.w.n

mel.e.n <- dplyr::summarise(dplyr::group_by(mel.e, pop.name),
                            n=n());mel.e.n

mel.w.n <- dplyr::summarise(dplyr::group_by(mel.w, pop.name),
                            n=n());mel.w.n

mel.n <- dplyr::summarise(dplyr::group_by(mel.w, pop.name),
                          n=n());mel.n
## save matrices
write.table(mel.ec.e.dist, "local/output/dist.vars/mel.ec.e.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.ec.w.dist, "local/output/dist.vars/mel.ec.w.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.ec.e.alt.dist, "local/output/e.vars/mel.ec.e.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.ec.w.alt.dist, "local/output/e.vars/mel.ec.w.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.ec.e.n, "local/output/dist.vars/mel.ec.e.dist.pop.n.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.ec.w.n, "local/output/dist.vars/mel.ec.w.dist.pop.n.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.e.dist, "local/output/dist.vars/mel.e.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.e.alt.dist, "local/output/e.vars/mel.e.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.w.dist, "local/output/dist.vars/mel.w.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.w.alt.dist, "local/output/e.vars/mel.w.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.dist, "local/output/dist.vars/mel.dist.pop.txt", row.names = TRUE, col.names = TRUE)
write.table(mel.alt.dist, "local/output/e.vars/mel.alt.dist.pop.txt", row.names = TRUE, col.names = TRUE)


##### maps ####
mean(mel.ec.e.loc$latitude);mean(mel.ec.e.loc$longitude)
mel.ec.e.map <- get_googlemap(center=c(lon=-77.44, lat=-0.5286628), zoom = 9, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")


gc(); p1 <- ggmap(mel.ec.e.map) +
  #theme_bw() +
  geom_jitter(data = mel.ec.e, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
               size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
    panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.ec.e.png")

## west 

mean(mel.ec.w.loc$latitude);mean(mel.ec.w.loc$longitude)
mel.ec.w.map <- get_googlemap(center=c(lon=-78.96914, lat=0.21412), zoom = 9, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(mel.ec.w.map) +
  #theme_bw() +
  geom_jitter(data = mel.ec.w, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.ec.w.png")


mean(mel.co.e$latitude);mean(mel.co.e$longitude)
mel.co.e.map <- get_googlemap(center=c(lon= -76.68982, lat=1.214120), zoom = 8, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(mel.co.e.map) +
  #theme_bw() +
  geom_jitter(data = mel.co.e, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.co.e.png")

## west 

mean(mel.co.w$latitude);mean(mel.co.w$longitude)
mel.co.w.map <- get_googlemap(center=c(lon=-77.08683, lat=3.87885), zoom = 8, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(mel.co.w.map) +
  #theme_bw() +
  geom_jitter(data = mel.co.w, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.co.w.png")



### all east ##

mean(mel.e.loc$latitude);mean(mel.e.loc$longitude)
mel.e.map <- get_googlemap(center=c(lon=-75.67008, lat=-0.4237762), zoom = 6, color="bw",  maptype = "terrain",
                              style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(mel.e.map) +
  #theme_bw() +
  geom_jitter(data = mel.e, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21,24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.e.png")

### all west ##

mean(mel.w.loc$latitude);mean(mel.w.loc$longitude)
mel.w.map <- get_googlemap(center=c(lon=-78.49544, lat=2.2), zoom = 7, color="bw",  maptype = "terrain",
                           style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
gc(); p1 <- ggmap(mel.w.map) +
  #theme_bw() +
  geom_jitter(data = mel.w, colour="black",
              aes(x = longitude, y = latitude, shape = type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21,24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.w.png")


### all mel ##

mean(mel.loc$latitude);mean(mel.loc$longitude)
mel.map <- get_googlemap(center=c(lon=-74.5, lat=0.5), zoom = 6, color="bw",  maptype = "terrain",
                           style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

gc(); p1 <- ggmap(mel.map) +
  #theme_bw() +
  geom_jitter(data = era, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  #scale_color_wes(name="type") +
  xlab("Longitude") + ylab("Latitude") +
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1

ggsave("local/output/plots/maps/mel.png")

#### by pop #####
names(mel)
mel.pop <-  summarise(group_by(mel, pop.name, pop.short),
                      n=n(),
                      latitude=mean(latitude),
                      longitude=mean(longitude),
                      alt.type=unique(str_sub(type, 2, 4)),
                      altitude.max=max(altitude),
                      altitude.min=min(altitude),
                      altitude=mean(altitude));mel.pop

gc(); p1 <- ggmap(mel.map) +
  #theme_bw() +
  geom_jitter(data = mel.pop, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  xlab("Longitude") + ylab("Latitude") +
  geom_text(data=mel.pop, aes(x=longitude, y=latitude,label=n),hjust=1.5, vjust=1.5, inherit.aes = FALSE)+
  scale_shape_manual(values = c(24,22,21))+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.pops.png")

gc(); p1 <- ggmap(mel.map) +
  #theme_bw() +
  geom_jitter(data = mel.pop, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  xlab("Longitude") + ylab("Latitude") +
  geom_text(data=mel.pop, aes(x=longitude, y=latitude,label=n),hjust=1.5, vjust=1.5, inherit.aes = FALSE)+
  scale_shape_manual(values = c(24,22,21))+
  ylim(-1.5,6.4)+xlim(-80.1,-75.3)+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.pops.zoom.png")

gc(); p1 <- ggmap(mel.map) +
  #theme_bw() +
  geom_jitter(data = mel.pop, colour="black",
              aes(x = longitude, y = latitude, shape = alt.type, fill=altitude, colour=altitude),
              size=8, alpha=.8, width = .02, height = .02)+
  xlab("Longitude") + ylab("Latitude") +
  geom_text(data=mel.pop, aes(x=longitude, y=latitude,label=n),hjust=1.5, vjust=1.5, inherit.aes = FALSE)+
  scale_shape_manual(values = c(24,22,21))+
  ylim(-1.5,2.1)+xlim(-80.1,-75.3)+
  scale_fill_gradient( low = "green", high = "blue")+
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(size = guide_legend(order=2)); p1
ggsave("local/output/plots/maps/mel.pops.zoom2.png")

