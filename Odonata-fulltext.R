# This script is to process the decisions made for the 62 articles screened
# at full text and generate plots of the metadata

# Functions and libraries ------------------------------------------------------
library(ggplot2) # for internal plotting functions
dragon_palette <- c("#aeb6ba", "#29216d", "#068775", "#16bb83", "#80c75d", "#1598c0")

# Read in the results of full text screening -----------------------------------

# Final 45 studies after independent screening at full text and adding in
# articles from Phase 1 full text screening
odonata <- read.csv("./included_studies.csv")

# Length of studies ------------------------------------------------------------
odonata$total_length <- as.numeric(odonata$LastYear) - as.numeric(odonata$FirstYear) + 1

lengths <- factor(append(odonata$total_length, 2:max(odonata$total_length)))[1:length(odonata$total_length)]
maxtimes <- factor(append(as.numeric(odonata$MaxTimeSeries), 1:max(odonata$total_length)))[1:length(odonata$total_length)]

par(pty="s", las=1)
barplot(as.numeric(maxtimes)[order(odonata$total_length, decreasing = T)], cex.names = .5, add=F, axes=F,
        col=paste(dragon_palette[3], "ff", sep=""), border = F, ylim=c(0,40), xlim=c(0,54), ylab="Time series length (years)", xlab="Studies")
z <- barplot(sort(odonata$total_length, decreasing = T), col=paste(dragon_palette[3], "55", sep=""), 
             border=F, add=T, axes=F)
legend("topright", legend=c("Maximum continuous time series", "Total length of time series including breaks"), bty="n",
       col=paste(dragon_palette[c(3,3)], c("ff", "55"), sep=""), pch=15, pt.cex=2)
axis(2, at=c(0,5,10,20,30,40))


# Study locations --------------------------------------------------------------
library(maptools)
data(wrld_simpl)

world <- spTransform(wrld_simpl, CRS("+proj=robin"))
world <- world[!world$NAME == "Antarctica", ]
z <- raster::extent(world)
z[3] <- -6084340
z[4] <- 8238578
worldclip <- raster::crop(world, z)

studysites <-
  unique(data.frame(as.numeric(odonata$mapLat), as.numeric(odonata$mapLong)))
names(studysites) <- c("latitude", "longitude")
sp::coordinates(studysites) <- ~ longitude + latitude
sp::proj4string(studysites) <- sp::proj4string(wrld_simpl)
studysites <- sp::spTransform(studysites, sp::proj4string(world))


blue <- paste(dragon_palette[6], "55", sep="")

plot(worldclip, bg=blue, border=blue, lty=0,
     col="white")
plot(worldclip, border = blue, lty=1, lwd=.25,
     col= paste(dragon_palette[1], "15", sep=""), add=T)

points(studysites,
       pch = 21,
       cex = 2, lwd=1, col="black",
       bg = paste(colorRampPalette(c("#16bb83", dragon_palette[2]))(4)[ceiling(log(odonata$total_length))], "aa", sep=""))

legend("topright",
       legend = c(36, NA, NA, NA, NA, 2),
       fill = rev(colorRampPalette(c("#16bb83", dragon_palette[2]))(6)),
       border = NA,
       y.intersp = 0.5,
       cex = 2, bty="n")
legend("topleft", "Study length", bty="n")




# Habitats ---------------------------------------------------------------------

habitat_classification <- read.csv("./odonata-habitats.csv")

habitat_dictionary <- topictagger::create_dictionary(topictagger::fill_rows(habitat_classification))

habitat_tags <- topictagger::tag_strictly(odonata$HabitatDescription, habitat_dictionary)

habitats <- names(habitat_dictionary)
habitat_toplevel <- array(dim=c(nrow(odonata), length(habitats)))
for(i in 1:length(habitats)){
  habitat_toplevel[,i] <- rowSums(habitat_tags[,grep(habitats[i], colnames(habitat_tags))])
}
colnames(habitat_toplevel) <- habitats

odonata$habitats <- apply(habitat_toplevel, 1, function(x){
  paste(names(x)[x>0], collapse=";")
})

odonata$habitats[odonata$habitats==""] <- "z_unknown"
odonata$habitats <- gsub("Forest", "z_unknown", odonata$habitats) # too few, set as other for plot

# Sampling methods -------------------------------------------------------------

# ############################################################################ #
# NOTE TO SELF: ADD SAMPLE METHODS AS LOCAL DATA -------------------------------
# and include extra terms for sweep netting that are currently missing
# and microhabitat samples for drop box
# ############################################################################ #

samplemethods <- read.csv("./insect-sampling-methods.csv")
samplemethods <- topictagger::fill_rows(samplemethods)

method_dictionary <- topictagger::create_dictionary(samplemethods, 
                                                    return_dictionary = F)

# Before creating the dictionary, we want to allow for stemming and word forms
# For example, we want "Surber sampler" to match "Surber samplers"
# For this, we need to use an asterisk
for(i in 1:length(method_dictionary)){
  tmp <- method_dictionary[[i]]
  for(j in 1:length(tmp)){
    tmpj <- tmp[[j]]
    for(k in 1:length(tmpj)){
      method_dictionary[[i]][[j]][[k]] <- gsub("\\*\\*", "*", paste(tmpj[[k]], "*", sep=""))
    }
    
  }
}

methods_tagged <- topictagger::tag_strictly(odonata$SamplingMethod, method_dictionary)

samplemethods <- tolower(names(method_dictionary))
methods_lookup <- array(dim=c(nrow(odonata), length(samplemethods)))


for(i in 1:length(samplemethods)){
  indexed_vals <- methods_tagged[,grep(samplemethods[i], tolower(colnames(methods_tagged))), drop=F]
  methods_lookup[,i] <- rowSums(indexed_vals)
}
odonata$SamplingMethod[rowSums(methods_lookup)==0]

methods_lookup <- data.frame(methods_lookup)
colnames(methods_lookup) <- samplemethods

odonata$Methods <- apply(methods_lookup, 1, function(x){paste(colnames(methods_lookup)[x>0], collapse="; ")})
odonata[odonata==""] <- NA
odonata$Methods[is.na(odonata$Methods)] <- "z_unknown"

# Heatmap of insect orders vs sampling methods ---------------------------------
tmpdat <- data.frame(cbind(odonata$Methods, odonata$habitats))

merged_dat <- unlist(apply(tmpdat, 1, function(x){
  x1 <- unlist(strsplit(x[1], ";"))
  y1 <- unlist(strsplit(x[2], ";"))
  paste(rep(x1, each=length(y1)), y1, sep=";;;")
}))


tabdat <- (data.frame(matrix(trimws(unlist(strsplit(trimws(merged_dat), ";;;"))), ncol=2, byrow=T)))
tabdat[,1] <- paste(tabdat[,1], " (", table(tabdat[,1])[factor(tabdat[,1])], ")", sep="")
tabdat[,2] <- paste(tabdat[,2], " (", table(tabdat[,2])[factor(tabdat[,2])], ")", sep="")
tabdat[,2] <- gsub("z_unknown", "z_other", tabdat[,2])


tmp <- table(tabdat)
tmp[tmp == 0] <- NA
tmp <- data.frame(tmp)

print(ggplot(tmp, aes(X1, X2, fill = Freq)) +
        geom_tile() +
        scale_fill_continuous(high = dragon_palette[2],
                              low = "#16bb8388",
                              na.value = "grey96") + # fill with color
        theme(
          line = element_blank(),
          # remove the background, tickmarks, etc
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 12
          ),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 9)
        ) +
        coord_equal())

