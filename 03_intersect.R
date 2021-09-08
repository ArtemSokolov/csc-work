library( tidyverse )
library( seriation )

nPop <- 40
clExcl <- c(-1,3,14)

## Load cluster assignments and gating results
X <- read_csv( "data/hdbscan.csv", col_types=cols() ) %>%
    rename( HDBSCAN = hdbscan ) %>%
    filter( !(HDBSCAN %in% clExcl) )
Y <- read_csv( "output/gating.csv", col_types=cols() ) %>%
    filter( Label %in% str_c("Pop",1:nPop) )
XY <- inner_join(X, Y, by="CellID")

## Isolate intersections that pass the cell count filter
Z <- XY %>% group_by(HDBSCAN, Label) %>%
    summarize( Count = n(), across(CD3:Vimentin, unique), .groups='drop' ) %>%
    filter( Count >= 15 ) %>% arrange( HDBSCAN, desc(Count) )

## Common theme elements
mytheme <- function()
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.minor = element_blank())

fplot <- function(.df, clusRange) {
    ZZ <- .df %>% filter(HDBSCAN %in% clusRange) %>%
        mutate( Slice = str_c(HDBSCAN, '_', Label) ) %>%
        mutate( Slice = factor(Slice,Slice) )

    ZZ1 <- ZZ %>% select(Slice, HDBSCAN, `# Cells` = Count)
    gg1 <- ggplot(ZZ1, aes(x=Slice, y=`# Cells`)) + theme_minimal() +
        geom_bar( stat='identity', color='darkgray', fill='lightgray' ) +
        scale_y_log10(labels=scales::comma, breaks=c(10,100,1000,10000,100000)) +
        mytheme() + facet_grid( ~HDBSCAN, scales='free', space='free' ) +
        theme(panel.grid.major.x = element_blank(),
              plot.margin=unit(c(0,0.5,0,0.5),unit="cm"))

    ZZ2 <- ZZ %>% select(-Label, -Count) %>%
        pivot_longer( -c(HDBSCAN, Slice), names_to="Marker" )
    gg2 <- ggplot( ZZ2, aes(x=Slice, y=Marker, color=value) ) +
        theme_minimal() + geom_point() +
        scale_color_manual(values=c(pos="black", neg="lightgray"), guide=FALSE ) +
        mytheme() + facet_grid( ~HDBSCAN, scales='free', space='free' ) +
        theme(panel.grid.major = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              plot.margin=unit(c(0,0.5,0,0.5),unit="cm"))

    egg::ggarrange( gg1, gg2, ncol=1, heights=c(1,2) )
}

gg1 <- fplot( Z, 0:11 )
gg2 <- fplot( Z, 12:18 )
ggsave( "plots/intersect1.png", gg1, width=8, height=3 )
ggsave( "plots/intersect2.png", gg2, width=8, height=3 )
