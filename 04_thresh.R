library( tidyverse )

load("models/GMMs.RData")

mypanel <- c("CD3", "CD45RO", "aSMA", "CD163", "CD20", "CD4", "CD45", "CD68", "CD8a",
             "Desmin", "Ecad", "FOXP3", "Keratin", "PCNA", "PD1", "PDL1", "Vimentin" )

P <- GMMs %>% filter(Marker %in% mypanel) %>% naivestates::GMMreshape()

## A marker is considered expressed if p >= thresh
##   not expressed if p <= thresh
##   and ambiguous otherwise
gate <- function( .df, thresh=0.8 ) {
    .df %>% mutate( across(-CellID, ~case_when(
                                       .x >= thresh ~ "pos",
                                       .x <= (1-thresh) ~ "neg",
                                       TRUE ~ "amb")) )
}

## Counts the fraction of ambiguous calls in each marker for a given treshold
fAmb <- function( .df, thresh ) {
    gate(.df, thresh) %>%
        summarize( across(-CellID, ~sum(.x == "amb")/length(.x)) ) %>%
        pivot_longer( everything(), names_to="Marker", values_to="FracAmb" )
}

A <- c(0.5, 0.6, 0.7, 0.8, 0.9) %>% set_names() %>%
    map( ~fAmb(P, .x) ) %>% bind_rows( .id="Threshold" ) %>%
    mutate( across(Threshold, as.numeric) )

ggplot( A, aes(x=Threshold, y=FracAmb) ) +
    theme_bw() + geom_line() +
    facet_wrap( ~Marker, nrow=3, ncol=6 ) +
    ylab( "Fraction Ambiguous" ) +
    theme(panel.grid.minor.x = element_blank(),
          strip.background = element_blank()) +
    ggsave( "threshold.png", width=8, height=5 )

G80 <- P %>% gate(0.8)
G80 %>% select(-CellID) %>% filter( if_all(.fns= ~.x == "amb") )

checkDesmin <- function() {
    X <- read_csv("data/WD-76845-097-ij_subtracted_50_qc.csv",
              col_types=cols( CellID=col_integer() ))

    Y <- read_csv( "data/hdbscan.csv", col_types=cols() ) %>%
        rename( HDBSCAN = hdbscan )
    XY <- inner_join(X, Y) %>% filter(HDBSCAN==4)
    Z <- bind_rows( All = select(X, Desmin_555_cellRingMask),
                   Clus4 = select(XY, Desmin_555_cellRingMask),
                   .id = "Source" )
    
    ggplot( Z, aes(x=Desmin_555_cellRingMask, color=Source) ) +
        geom_density(adjust=0.1) + theme_bw() +
        scale_color_manual( values=c(All = 'black', Clus4 = 'tomato'),
                           labels=c(All = 'All (962,343 cells)',
                                    Clus4 = 'Cluster4 (2,300 cells)') ) +
        ggsave("Desmin.png", width=8, height=5 )
}
