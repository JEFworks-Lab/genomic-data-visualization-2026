data <- read.csv('~/Desktop/genomic-data-visualization-2026/data/Visium-IRI-ShamR_matrix.csv.gz')
pos <- data[,c('x', 'y')]
rownames(pos) <- data[,1]
gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
head(pos)
gexp[1:5,1:5]
dim(gexp)

library(ggplot2)
library(gganimate)
df <- data.frame(pos, gene=gexp[, 'Gpx3'])
head(df)
ggplot(df, aes(x=x, y=y, col=gene)) + 
  geom_point(size=0.5) +
  coord_fixed()

head(pos)

points <- data.frame(pos)

# Compute the center (mean of x and y coordinates)
center_x <- mean(points$x)
center_y <- mean(points$y)

?cat
cat("Center coordinates:\n")
cat("x:", center_x, "\n")
cat("y:", center_y, "\n\n")

# Compute distance from each point to the center
distances <- sqrt((points$x - center_x)^2 + (points$y - center_y)^2)
length(distances)
dim(points)

# Add distances to the dataframe
points$distance_to_center <- distances
head(points)

df <- data.frame(pos, gene=gexp[, 'Gpx3'])
head(df)
g1 <- ggplot(df, aes(x=x, y=y, col=gene)) + 
  geom_point(size=0.5) +
  coord_fixed()

df <- data.frame(points, gene=gexp[, 'Gpx3'])
head(df)
g2 <- ggplot(df, aes(x=distance_to_center, y=gene, col=gene)) + 
  geom_point(size=0.5) 

library(patchwork)
g1 + g2

## AI output
library(ggplot2)
library(gganimate)
library(dplyr)

# State 1: spatial plot
df1 <- data.frame(pos, gene=gexp[, 'Gpx3']) %>%
  mutate(state = "Spatial View",
         x_plot = x,
         y_plot = y)
head(df1)

# State 2: distance plot
df2 <- data.frame(points, gene=gexp[, 'Gpx3']) %>%
  mutate(state = "Distance View",
         x_plot = distance_to_center,
         y_plot = gene)
head(df2)

# Combine datasets
df_combined <- bind_rows(df1, df2)
head(df_combined)

# Create animated plot
anim <- ggplot(df_combined, aes(x=x_plot, y=y_plot, col=gene)) + 
  geom_point(size=0.5) +
  transition_states(state, 
                    transition_length = 2,
                    state_length = 1) +
  labs(title = '{closest_state}',
       x = NULL,
       y = NULL) +
  ease_aes('cubic-in-out') +
  view_follow(fixed_y = FALSE, fixed_x = FALSE) +
  theme_minimal()

# Render animation
animate(anim, nframes = 100, fps = 10, width = 300, height = 300)

?gifski
## If gifski does not work for you and instead you get a bunch of PNGs not stitched together to a gif
## try the ImageMagick renderer instead
#install.packages('magick')
library(magick)
animate(anim, nframes = 100, fps = 10, width = 300, height = 300, 
        renderer = magick_renderer())

## debugging by splitting in two
## first make ggplot
anim <- ggplot(df_combined, aes(x=x_plot, y=y_plot, col=gene)) + 
  geom_point(size=2) 
anim ## preview both plots on top of each other
## then add animation 
anim +
  transition_states(state, 
                    transition_length = 2,
                    state_length = 1) +
  labs(title = '{closest_state}',
       x = NULL,
       y = NULL) +
  ease_aes('bounce-in') + ## how to transition between the two plots
  view_follow(fixed_y = FALSE, fixed_x = FALSE) +
  theme_minimal()

## if taking a long time, subsampling or downsampling data
## figure out what you like, then apply to full data

