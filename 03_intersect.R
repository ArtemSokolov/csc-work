library( tidyverse )

## Load cluster assignments and gating results
X <- read_csv( "data/hdbscan.csv", col_types=cols() ) %>%
    rename( HDBSCAN = hdbscan ) %>% filter( HDBSCAN != -1 )
Y <- read_csv( "output/gating-80.csv", col_types=cols() )

## Compute the representation of each subpopulation in each cluster
XY <- inner_join(X, Y, by="CellID") %>%
    count(across(-CellID), name="nCells") %>%
    group_by(HDBSCAN) %>%
    mutate( Frac = nCells / sum(nCells),
            Pct = scales::percent(Frac, accuracy=0.01)) %>%
    select( nCells, Frac, Pct, everything() ) %>% ungroup()

## Isolate intersections that pass the 1% filter
Z <- XY %>% filter(Frac >= 0.01, nCells >= 5) %>%
    arrange( HDBSCAN, desc(nCells) )

## Common theme elements
mytheme <- function()
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.minor = element_blank())

## Plots a subset of clusters
fplot <- function(.df, clusRange) {
    ZZ <- .df %>% filter(HDBSCAN %in% clusRange) %>%
        mutate(Slice = 1:n(),
               Slice = factor(Slice,Slice))

    gg1 <- ggplot(ZZ, aes(x=Slice, y=nCells)) +
        theme_minimal() + ylab("# Cells") +
        geom_bar( stat='identity', color='darkgray', fill='lightgray' ) +
        geom_text( aes(label=Pct), size=2, angle=90, hjust=1 ) +
        scale_y_log10(labels=scales::comma, breaks=c(10,100,1000,10000,100000)) +
        mytheme() + facet_grid( ~HDBSCAN, scales='free', space='free' ) +
        theme(panel.grid.major.x = element_blank(),
              plot.margin=unit(c(0,0.5,0,0.5),unit="cm"))

    ZZ2 <- ZZ %>% select(-Frac, -nCells, -Pct) %>%
        pivot_longer( -c(HDBSCAN, Slice), names_to="Marker" )
    gg2 <- ggplot( ZZ2, aes(x=Slice, y=Marker, color=value, fill=value) ) +
        theme_minimal() + geom_point(shape=21) +
        scale_fill_manual(values=c(pos="black", amb="lightgray", neg="white"), guide=FALSE) +
        scale_color_manual(values=c(pos="black", amb="darkgray", neg="lightgray"), guide=FALSE) +
        mytheme() + facet_grid( ~HDBSCAN, scales='free', space='free' ) +
        theme(panel.grid.major = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              plot.margin=unit(c(0,0.5,0,0.5),unit="cm"))

    egg::ggarrange( gg1, gg2, ncol=1, heights=c(1,2) )
}

ggsave( "plots/intersect1.png", fplot( Z, 0:6 ), width=8, height=3 )
ggsave( "plots/intersect2.png", fplot( Z, 7:13 ), width=8, height=3 )
ggsave( "plots/intersect3.png", fplot( Z, 14:18 ), width=8, height=3 )
