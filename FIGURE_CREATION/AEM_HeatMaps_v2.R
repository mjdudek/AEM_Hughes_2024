## FLAVRS Builder
# Build FLAVRS tbrke for meteor

# Developed by: Marissa Dudek (1)
# (1) University of North Carolina at Chapel Hill, Geological Sciences

# ---------------------------------------------------------------------

# Import library
library(ggplot2)
library(treemap)
library(ggplot2)

# ---------------------------------------------------------------------

# Set column names
colnames <- c('Density.kg/m3','IntDiam.m','IntMass.kg','IntVel.kms','IntAng.deg','FinDiam.m','FinMass.kg','MassR.perc','FinVel.kms','FinAng.deg','CraterDiamt.m','CraterDiamf.m','CraterDep.m','PanFact','AltBurst.m')

# Import lunar data
L_stat <- as.data.frame(read.csv("statsL_ablation.csv")); colnames(L_stat) <- colnames

# Import ablation data
M_abl <- read.csv("statsM_ablation.csv"); colnames(M_abl) <- colnames
E_abl <- read.csv("statsE_ablation.csv"); colnames(E_abl) <- colnames
V_abl <- read.csv("statsV_ablation.csv"); colnames(V_abl) <- colnames

# Import ablation + breakup data
M_brk <- read.csv("statsM_ablation-break.csv"); colnames(M_brk) <- colnames
E_brk <- read.csv("statsE_ablation-break.csv"); colnames(E_brk) <- colnames
V_brk <- read.csv("statsV_ablation-break.csv"); colnames(V_brk) <- colnames

# Remove zeros

# Plot styles
#p1 = ggplot(as.data.frame(cbind(x1,y1)), aes(x=x1, y=y1) ) +
#  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') +
#  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#  scale_fill_distiller(palette= "Spectral", direction=-1, limits=c(0,0.03)) +
#  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
#  scale_x_continuous(labels=scales::comma) +
#  theme_bw ()
#p1
#p1 = ggplot(as.data.frame(cbind(x1,y1)), aes(x=x1, y=y1) ) +
#  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') +
#  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
#  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
#  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
#  scale_x_continuous(limits=c(0,45), labels=scales::comma) +
#  theme_bw ()
#p1

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Initial Diameter vs. Crater Diameter 

# Mars
# Pull data
x1 = na.omit(M_abl$IntVel.kms/1000)
y1 = na.omit(M_abl$CraterDiam.m/1000)
# Plot 2d density
p1 = ggplot(as.data.frame(cbind(x1,y1)), aes(x=x1, y=y1) ) +
  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,45), labels=scales::comma) +
  theme_bw ()
p1

# Earth
# Pull data
x2 = na.omit(E_abl$IntVel.kms/1000)
y2 = na.omit(E_abl$CraterDiam.m/1000)
# Plot 2d density
p2 = ggplot(as.data.frame(cbind(x2,y2)), aes(x=x2, y=y2) ) +
  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,45), labels=scales::comma) +
  theme_bw ()
p2

# Venus
# Pull data
x3 = na.omit(V_abl$IntVel.kms/1000)
y3 = na.omit(V_abl$CraterDiam.m/1000)
# Plot 2d density
p3 = ggplot(as.data.frame(cbind(x3,y3)), aes(x=x3, y=y3) ) + ylim(0,350) +
  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + 
  geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,45), labels=scales::comma) +
  theme_bw ()
p3

# ---------------------------------------------------------------------
# Mass vs. Crater Diameter 

# Mars
# Pull data
x4 = na.omit(M_abl$IntMass.kg)
y4 = na.omit(M_abl$CraterDiam.m/1000)
# Plot 2d density
p4 = ggplot(as.data.frame(cbind(x4,y4)), aes(x=x4, y=y4) ) + ylim(0,350) +
  labs(x='Initial Projectile Mass (kg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + 
  geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(trans='log10', labels = function(x) format(x, scientific = TRUE)) +
  theme_bw ()
p4

# Earth
# Pull data
x5 = na.omit(E_abl$IntMass.kg)
y5 = na.omit(E_abl$CraterDiam.m/1000)
# Plot 2d density
p5 = ggplot(as.data.frame(cbind(x5,y5)), aes(x=x5, y=y5) ) +
  labs(x='Initial Projectile Mass (kg)', y='Final Crater Diameter (km)') + 
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + 
  geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(trans='log10', labels = function(x) format(x, scientific = TRUE)) +
  theme_bw ()
p5

# Venus
# Pull data
x6 = na.omit(V_abl$IntMass.kg)
y6 = na.omit(V_abl$CraterDiam.m/1000)
# Plot 2d density
p6 = ggplot(as.data.frame(cbind(x6,y6)), aes(x=x6, y=y6) ) + ylim(0,350) +
  labs(x='Initial Projectile Mass (kg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(trans='log10', labels = function(x) format(x, scientific = TRUE)) +
  theme_bw ()
p6

# ---------------------------------------------------------------------
# Angle vs. Crater Diameter 

# Mars
# Pull data
x7 = na.omit(M_abl$IntAng.deg)
y7 = na.omit(M_abl$CraterDiam.m/1000)
# Plot 2d density
p7 = ggplot(as.data.frame(cbind(x7,y7)), aes(x=x7, y=y7) ) + ylim(0,350) +
  labs(x='Initial Projectile Angle (deg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.075, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,90), labels=scales::comma) +
  theme_bw ()
p7

# Earth
# Pull data
x8 = na.omit(E_abl$IntAng.deg)
y8 = na.omit(E_abl$CraterDiam.m/1000)
# Plot 2d density
p8 = ggplot(as.data.frame(cbind(x8,y8)), aes(x=x8, y=y8) ) +
  labs(x='Initial Projectile Angle (deg)', y='Final Crater Diameter (km)') + 
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.075, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,95), labels=scales::comma) +
  theme_bw ()
p8

# Venus
# Pull data
x9 = na.omit(V_abl$IntAng.deg)
y9 = na.omit(V_abl$CraterDiam.m/1000)
# Plot 2d density
p9 = ggplot(as.data.frame(cbind(x9,y9)), aes(x=x9, y=y9) ) + ylim(0,350) +
  labs(x='Initial Projectile Angle (deg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.075, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,95), labels=scales::comma) +
  theme_bw ()
p9

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Misc. plots

# Venus
# Pull data
x10 = na.omit(V_abl$IntVel.kms/1000)
y10 = na.omit(V_abl$IntMass.kg)
z10 = na.omit(V_abl$CraterDiam.m/1000)
# Plot 2d density
p10 = ggplot(as.data.frame(cbind(x10,y10)), aes(x=x10, y=y10) ) +
  labs(x='Initial Projectile Velocity (km/s)', y='Initial Projectile Mass (kg)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(aes(size = z10), alpha=0.075, colour = 'black') +
  scale_y_continuous(trans='log10', labels=scales::comma, labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(labels=scales::comma) +
  theme_bw ()
p10

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Initial Diameter vs. Crater Diameter 

# Mars
# Pull data
x1 = na.omit(M_brk$IntVel.kms/1000)
y1 = na.omit(M_brk$CraterDiam.m/1000)
# Plot 2d density
p1 = ggplot(as.data.frame(cbind(x1,y1)), aes(x=x1, y=y1) ) +
  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,45), labels=scales::comma) +
  theme_bw ()
p1

# Earth
# Pull data
x2 = na.omit(E_brk$IntVel.kms/1000)
y2 = na.omit(E_brk$CraterDiam.m/1000)
# Plot 2d density
p2 = ggplot(as.data.frame(cbind(x2,y2)), aes(x=x2, y=y2) ) +
  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') + 
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,45), labels=scales::comma) +
  theme_bw ()
p2

# Venus
# Pull data
x3 = na.omit(V_brk$IntVel.kms/1000)
y3 = na.omit(V_brk$CraterDiam.m/1000)
# Plot 2d density
p3 = ggplot(as.data.frame(cbind(x3,y3)), aes(x=x3, y=y3) ) + ylim(0,350) +
  labs(x='Initial Projectile Velocity (km/s)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,45), labels=scales::comma) +
  theme_bw ()
p3

# ---------------------------------------------------------------------
# Mass vs. Crater Diameter 

# Mars
# Pull data
x4 = na.omit(M_brk$IntMass.kg)
y4 = na.omit(M_brk$CraterDiam.m/1000)
# Plot 2d density
p4 = ggplot(as.data.frame(cbind(x4,y4)), aes(x=x4, y=y4) ) + ylim(0,350) +
  labs(x='Initial Projectile Mass (kg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(trans='log10', labels = function(x) format(x, scientific = TRUE)) +
  theme_bw ()
p4

# Earth
# Pull data
x5 = na.omit(E_brk$IntMass.kg)
y5 = na.omit(E_brk$CraterDiam.m/1000)
# Plot 2d density
p5 = ggplot(as.data.frame(cbind(x5,y5)), aes(x=x5, y=y5) ) +
  labs(x='Initial Projectile Mass (kg)', y='Final Crater Diameter (km)') + 
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(trans='log10', labels = function(x) format(x, scientific = TRUE)) +
  theme_bw ()
p5

# Venus
# Pull data
x6 = na.omit(V_brk$IntMass.kg)
y6 = na.omit(V_brk$CraterDiam.m/1000)
# Plot 2d density
p6 = ggplot(as.data.frame(cbind(x6,y6)), aes(x=x6, y=y6) ) + ylim(0,350) +
  labs(x='Initial Projectile Mass (kg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.025, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(trans='log10', labels = function(x) format(x, scientific = TRUE)) +
  theme_bw ()
p6

# ---------------------------------------------------------------------
# Angle vs. Crater Diameter 

# Mars
# Pull data
x7 = na.omit(M_brk$IntAng.deg)
y7 = na.omit(M_brk$CraterDiam.m/1000)
# Plot 2d density
p7 = ggplot(as.data.frame(cbind(x7,y7)), aes(x=x7, y=y7) ) + ylim(0,350) +
  labs(x='Initial Projectile Angle (deg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.075, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,90), labels=scales::comma) +
  theme_bw ()
p7

# Earth
# Pull data
x8 = na.omit(E_brk$IntAng.deg)
y8 = na.omit(E_brk$CraterDiam.m/1000)
# Plot 2d density
p8 = ggplot(as.data.frame(cbind(x8,y8)), aes(x=x8, y=y8) ) +
  labs(x='Initial Projectile Angle (deg)', y='Final Crater Diameter (km)') + 
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.075, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,95), labels=scales::comma) +
  theme_bw ()
p8

# Venus
# Pull data
x9 = na.omit(V_brk$IntAng.deg)
y9 = na.omit(V_brk$CraterDiam.m/1000)
# Plot 2d density
p9 = ggplot(as.data.frame(cbind(as.factor(x9),y9)), aes(x=x9, y=y9) ) + ylim(0,350) +
  labs(x='Initial Projectile Angle (deg)', y='Final Crater Diameter (km)') +
  stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1, limits=c(0,1)) + geom_point(alpha=0.075, colour = 'black') +
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  scale_x_continuous(limits=c(0,95), labels=scales::comma) +
  theme_bw ()
p9

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

d.colors <- c("2700" = 'darkblue', "3300" = 'chartreuse4', "4500" = 'chocolate1', "7800" = 'brown4')

x10 = na.omit(L_stat$IntAng.deg)
y10 = na.omit(L_stat$CraterDiam.m/1000)

df.angle <- data.frame(
  "x" = factor(c(x10, x7, x8, x9)),
  "y" = c(y10, y7, y8, y9),
  "p" = factor(
    c(
      rep("Moon", length(y10)),
      rep("Mars", length(y7)), 
      rep("Earth", length(y8)), 
      rep("Venus", length(y9))), 
    levels = c("Moon", "Mars", "Earth", "Venus"),
    labels = c("Moon", "Mars", "Earth", "Venus"))
)

# Basic violin plot

p <- ggplot(df.angle, aes(x=x, y=y, fill = p) )  + 
  theme_bw()+
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  geom_violin(position = 'identity', alpha = 0.25, colour = 'black')+
  scale_fill_manual(values= c("Moon" = 'gray', "Mars" = 'red', "Earth" = 'blue', "Venus" = 'orange'), 
                      labels = levels(df.velma$p))+
  theme(legend.position = "none", aspect.ratio = 0.75)+
  facet_wrap(~p, nrow = 4)+
  labs(x = "Angle (deg.)", y = "Final Crater Diameter (km)")
p

# ---------------------------------------------------------------------

x12 <- na.omit(L_stat$IntMass.kg)

df.mass <- data.frame(
  "x" = c(x12, x4, x5, x6),
  "y" = c(y10, y7, y8, y9),
  "p" = factor(
    c(
      rep("Moon", length(y10)),
      rep("Mars", length(y7)), 
      rep("Earth", length(y8)), 
      rep("Venus", length(y9))), 
    levels = c("Moon", "Mars", "Earth", "Venus"),
    labels = c("Moon", "Mars", "Earth", "Venus")),
  "d" = factor(
    c(
    na.omit(L_stat$`Density.kg/m3`),
    na.omit(M_brk$`Density.kg/m3`),
    na.omit(E_brk$`Density.kg/m3`),
    na.omit(V_brk$`Density.kg/m3`))
  )
)

p.mass <- ggplot(df.mass, aes(x=x, y=y, colour = p) )  + 
  theme_bw()+
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma)+
  scale_x_continuous(trans= 'log10', limits = c(10^6, 10^16))+
  geom_point(alpha = 0.08)+
  scale_colour_manual(values= c("Moon" = 'gray', "Mars" = 'red', "Earth" = 'blue', "Venus" = 'orange'), 
                    labels = levels(df.velma$p))+
  theme(aspect.ratio = 0.75,
        axis.text.y = element_blank(), 
        legend.position = 'none')+
  facet_wrap(~p, nrow = 4)+
  labs(x = "Mass (kg)", y = NULL, colour = expression(paste("Density (kg/m"^"3", ")")))
p.mass

# ---------------------------------------------------------------------

x13 <- na.omit(L_stat$IntVel.kms/1000)

df.angle <- data.frame(
  "x" = factor(c(x13, x1, x2, x3)),
  "y" = c(y10, y1, y2, y3),
  "p" = factor(
    c(
      rep("Moon", length(y10)),
      rep("Mars", length(y1)), 
      rep("Earth", length(y2)), 
      rep("Venus", length(y3))), 
    levels = c("Moon", "Mars", "Earth", "Venus"),
    labels = c("Moon", "Mars", "Earth", "Venus"))
)

# Basic violin plot

p <- ggplot(df.angle, aes(x=x, y=y, fill = p) )  + 
  theme_bw()+
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma) +
  geom_violin(position = 'identity', alpha = 0.25, colour = 'black')+
  scale_fill_manual(values= c("Moon" = 'gray', "Mars" = 'red', "Earth" = 'blue', "Venus" = 'orange'), 
                    labels = levels(df.velma$p))+
  theme(legend.position = "none", aspect.ratio = 0.75)+
  facet_wrap(~p, nrow = 4)+
  labs(x = "Velocity (km/s)", y = "Final Crater Diameter (km)")
p

# ---------------------------------------------------------------------

x14 <- na.omit(L_stat$PanFact)
x15 = na.omit(M_brk$PanFact)
x16 = na.omit(E_brk$PanFact)
x17 = na.omit(M_brk$PanFact)


df.panfact <- data.frame(
  "x" = c(x14, x15, x16, x17),
  "y" = c(y10, y1, y2, y3),
  "p" = factor(
    c(
      rep("Moon", length(y10)),
      rep("Mars", length(y7)), 
      rep("Earth", length(y8)), 
      rep("Venus", length(y9))), 
    levels = c("Moon", "Mars", "Earth", "Venus"),
    labels = c("Moon", "Mars", "Earth", "Venus")),
  "d" = factor(
    c(
      na.omit(L_stat$PanFact),
      na.omit(M_brk$PanFact),
      na.omit(E_brk$PanFact),
      na.omit(M_brk$PanFact)
  )
))

p.panfact <- ggplot(df.panfact, aes(x=x, y=y, colour = p) )  + 
  theme_bw()+
  scale_y_continuous(trans='log10', limits=c(0.01,350), labels=scales::comma)+
  scale_x_continuous(limits = c(0, 50))+
  geom_point(alpha = 0.25)+
  scale_colour_manual(values= c("Moon" = 'gray', "Mars" = 'red', "Earth" = 'blue', "Venus" = 'orange'), 
                      labels = levels(df.velma$p))+
  theme(aspect.ratio = 0.75,
        axis.text.y = element_blank(), 
        legend.position = 'none')+
  facet_wrap(~p, nrow = 4)+
  labs(x = "Pancake Factor", y = NULL, colour = expression(paste("Density (kg/m"^"3", ")")))
p.panfact

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Final diameter
# ablation
mean(L_stat$FinDiam.m[L_stat$FinDiam.m != 0])
mean(M_abl$FinDiam.m[M_abl$FinDiam.m != 0])
mean(E_abl$FinDiam.m[E_abl$FinDiam.m != 0])
mean(V_abl$FinDiam.m[V_abl$FinDiam.m != 0])
# ablation and breakup
mean(L_stat$FinDiam.m[L_stat$FinDiam.m != 0])
mean(M_brk$FinDiam.m[M_brk$FinDiam.m != 0])
mean(E_brk$FinDiam.m[E_brk$FinDiam.m != 0])
mean(V_brk$FinDiam.m[V_brk$FinDiam.m != 0])

# Final angle
# ablation
mean(L_stat$FinAng.deg[L_stat$FinAng.deg != 0])
mean(M_abl$FinAng.deg[M_abl$FinAng.deg != 0])
mean(E_abl$FinAng.deg[E_abl$FinAng.deg != 0])
mean(V_abl$FinAng.deg[V_abl$FinAng.deg != 0])
# ablation and breakup
mean(L_stat$FinAng.deg[L_stat$FinAng.deg != 0])
mean(M_brk$FinAng.deg[M_brk$FinAng.deg != 0])
mean(E_brk$FinAng.deg[E_brk$FinAng.deg != 0])
mean(V_brk$FinAng.deg[V_brk$FinAng.deg != 0])

# Final velocity
# ablation
mean(L_stat$FinVel.kms[L_stat$FinVel.kms != 0])
mean(M_abl$FinVel.kms[M_abl$FinVel.kms != 0])
mean(E_abl$FinVel.kms[E_abl$FinVel.kms != 0])
mean(V_abl$FinVel.kms[V_abl$FinVel.kms != 0])
# ablation and breakup
mean(L_stat$FinVel.kms[L_stat$FinVel.kms != 0])
mean(M_brk$FinVel.kms[M_brk$FinVel.kms != 0])
mean(E_brk$FinVel.kms[E_brk$FinVel.kms != 0])
mean(V_brk$FinVel.kms[V_brk$FinVel.kms != 0])

# Final mass
# ablation
mean(L_stat$FinMass.kg[L_stat$FinMass.kg != 0]/L_stat$IntMass.kg)
mean(M_abl$FinMass.kg[M_abl$FinMass.kg != 0]/M_abl$IntMass.kg[M_abl$FinMass.kg != 0])
mean(E_abl$FinMass.kg[E_abl$FinMass.kg != 0]/E_abl$IntMass.kg[E_abl$FinMass.kg != 0])
mean(V_abl$FinMass.kg[V_abl$FinMass.kg != 0]/V_abl$IntMass.kg[V_abl$FinMass.kg != 0])
# ablation and breakup
mean(L_stat$FinMass.kg[L_stat$FinMass.kg != 0]/L_stat$IntMass.kg)
mean(M_brk$FinMass.kg[M_brk$FinMass.kg != 0]/M_brk$IntMass.kg[M_brk$FinMass.kg != 0])
mean(E_brk$FinMass.kg[E_brk$FinMass.kg != 0]/E_brk$IntMass.kg[E_brk$FinMass.kg != 0])
mean(V_brk$FinMass.kg[V_brk$FinMass.kg != 0]/V_brk$IntMass.kg[V_brk$FinMass.kg != 0])

# Crater stats
# ablation
summary(L_stat$CraterDiamf.m[L_stat$CraterDiamf.m != 0])
summary(M_abl$CraterDiamf.m[M_abl$CraterDiamf.m != 0])
summary(E_abl$CraterDiamf.m[E_abl$CraterDiamf.m != 0])
summary(V_abl$CraterDiamf.m[V_abl$CraterDiamf.m != 0])
# ablation and breakup
summary(L_stat$CraterDiamf.m[L_stat$CraterDiamf.m != 0])
summary(M_brk$CraterDiamf.m[M_brk$CraterDiamf.m != 0])
summary(E_brk$CraterDiamf.m[E_brk$CraterDiamf.m != 0])
summary(V_brk$CraterDiamf.m[V_brk$CraterDiamf.m != 0])

# Density
# ablation
summary(L_stat$`Density.kg/m3`[L_stat$CraterDiamf.m != 0])
summary(M_abl$`Density.kg/m3`[M_abl$CraterDiamf.m != 0])
summary(E_abl$`Density.kg/m3`[E_abl$CraterDiamf.m != 0])
summary(V_abl$`Density.kg/m3`[V_abl$CraterDiamf.m != 0])
# ablation and breakup
summary(L_stat$`Density.kg/m3`[L_stat$CraterDiamf.m != 0])
summary(M_brk$`Density.kg/m3`[M_brk$CraterDiamf.m != 0])
summary(E_brk$`Density.kg/m3`[E_brk$CraterDiamf.m != 0])
summary(V_brk$`Density.kg/m3`[V_brk$CraterDiamf.m != 0])

# Initial diameter
# ablation
unique(L_stat$IntDiam.m[L_stat$CraterDiamf.m != 0])
unique(M_abl$IntDiam.m[M_abl$CraterDiamf.m != 0])
unique(E_abl$IntDiam.m[E_abl$CraterDiamf.m != 0])
unique(V_abl$IntDiam.m[V_abl$CraterDiamf.m != 0])
# ablation and breakup
unique(L_stat$IntDiam.m[L_stat$CraterDiamf.m != 0])
unique(M_brk$IntDiam.m[M_brk$CraterDiamf.m != 0])
unique(E_brk$IntDiam.m[E_brk$CraterDiamf.m != 0])
unique(V_brk$IntDiam.m[V_brk$CraterDiamf.m != 0])

# Number of zeros
colSums(V_brk==0)


