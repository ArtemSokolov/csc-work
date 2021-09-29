library( tidyverse )

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
