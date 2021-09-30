library( tidyverse )

X <- read_csv("data/WD-76845-097-ij_subtracted_50_qc.csv",
              col_types=cols( CellID=col_integer() ))
Y <- read_csv( "data/hdbscan.csv", col_types=cols() ) %>%
    rename( HDBSCAN = hdbscan ) %>% filter( HDBSCAN != -1 )

XY <- inner_join(X, Y, by="CellID")
S <- XY %>% count(HDBSCAN, name="nCells") %>%
    mutate( Lbl = str_c("#Cells:", scales::comma(nCells), "  ") )

etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
ebl <- element_blank

fplot <- function( marker, cl, nc ) {
    Z <- XY %>% filter( HDBSCAN %in% cl )
    S1 <- S %>% filter( HDBSCAN %in% cl )

    ggplot( Z, aes(x={{marker}}) ) +
        geom_density( data=X, color='darkgray', lwd=1.2 ) +
        geom_density( color='red', lwd=1 ) + theme_bw() +
        geom_text( data=S1, aes(label=Lbl), x=Inf, y=Inf, hjust=1, vjust=1.1 ) +
        facet_wrap( ~HDBSCAN, scales='free_y', ncol=nc ) +
        ylab( "Density" ) + xlab( "Normalized Expression" ) +
        theme( strip.background = ebl(), strip.text = etxt(12),
              axis.ticks.y = ebl(), axis.text.y = ebl(),
              axis.title = etxt(14), axis.text.x = etxt(12) )
}


## Init the plot directory
dir.create("plots/04", showWarnings=FALSE, recursive=TRUE)

## Desmin - all
ncl <- unique(XY$HDBSCAN) %>% length
nc <- ceiling( sqrt(ncl) )
fplot( Desmin_555_cellRingMask, 0:max(XY$HDBSCAN), nc ) +
    ggsave( "plots/04/Desmin-all.png" )

## Desmin - specific clusters
fplot( Desmin_555_cellRingMask, c(3,4,7,15), 2 ) +
    ggsave( "plots/04/Desmin-3-4-7-15.png" )

