library(tidyverse)

## Fit bi-modal GMMs across all markers
if( !file.exists("models/GMMs.RData") ) {

    ## Load QCd expression
    X <- read_csv("data/consensus_clustering.csv.gz",
                  col_types=cols( CellID=col_integer() ))

    ## Columns to fit
    vcol <- list(CD3         = "anti_CD3_nucleiRingMask",
                 CD45RO      = "anti_CD45RO_nucleiRingMask",
                 aSMA        = "aSMA_660_nucleiRingMask",
                 Keratin     = "Keratin_570_nucleiRingMask",
                 CD4         = "CD4_488_nucleiRingMask",
                 CD45        = "CD45_PE_nucleiRingMask",
                 PD1         = "PD1_647_nucleiRingMask",
                 CD20        = "CD20_488_nucleiRingMask",
                 CD68        = "CD68_555_nucleiRingMask",
                 CD8a        = "CD8a_660_nucleiRingMask",
                 CD163       = "CD163_488_nucleiRingMask",
                 FOXP3       = "FOXP3_570_nucleiRingMask",
                 PDL1        = "PDL1_647_nucleiRingMask",
                 Ecad        = "Ecad_488_nucleiRingMask",
                 Vimentin    = "Vimentin_555_nucleiRingMask",
                 CDX2        = "CDX2_647_nucleiRingMask",
                 LaminABC    = "LaminABC_488_nucleiRingMask",
                 Desmin      = "Desmin_555_nucleiRingMask",
                 CD31        = "CD31_647_nucleiRingMask",
                 PCNA        = "PCNA_488_nucleiRingMask",
                 CollagenIV  = "CollagenIV_647_nucleiRingMask")

    ## Compute fits and cache to disk
    GMMs <- naivestates::GMMfit(X, CellID, !!!vcol)
    dir.create("models", showWarnings=FALSE)
    save( GMMs, file="models/GMMs.RData" )
    
} else {
    load("models/GMMs.RData")
}

## Bold element_text of desired font size s
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
ebl <- element_blank

## Init the plot directory
dir.create("plots/01/markers", showWarnings=FALSE, recursive=TRUE)

## Downsample to 100k points
npt <- nrow(GMMs$Values[[1]])
vs <- sample(1:npt, 1e5)

## Expression overview
Y <- GMMs %>% mutate_at( "Values", map, slice, vs ) %>%
    select( Marker, Values ) %>% unnest( Values ) %>%
    filter( AdjVal >= 0, AdjVal <= 1 )
gg_expr <- ggplot( Y, aes(x=AdjVal) ) + theme_bw() +
    ylab( "Density" ) + xlab( "Normalized Expression" ) +
    geom_density(linewidth=1.1) + facet_wrap( ~Marker, scales="free_y", nrow=4, ncol=6 ) +
    scale_x_continuous( breaks=c(0,0.5,1) ) +
    theme( axis.text.x = etxt(12), axis.title = etxt(14), strip.text = etxt(12),
          axis.text.y = ebl(), strip.background = ebl(), axis.ticks.y=ebl() )
ggsave( "plots/01/expression.png", gg_expr, width=12, height=7 )

## Fit overview
ggf <- naivestates::plotFitOverview( GMMs )
ggsave( "plots/01/overview.png", ggf, width=8, height=8 )

## Fit of individual markers
walk(GMMs$Marker, ~ggsave(
                     str_c("plots/01/markers/", .x, ".png"),
                     naivestates::plotMarker(GMMs, .x)
                 ))

## Marker panel of interest
mypanel <- c("CD3", "CD45RO", "aSMA", "CD163", "CD20", "CD4", "CD45", "CD68", "CD8a",
             "Desmin", "Ecad", "FOXP3", "Keratin", "PCNA", "PD1", "PDL1", "Vimentin" )
stopifnot( all(mypanel %in% GMMs$Marker) )

## Extract probabilities of expression
dir.create("output", showWarnings=FALSE)
P <- GMMs %>% filter(Marker %in% mypanel) %>%
    naivestates::GMMreshape()
write_csv( P, "output/probs.csv.gz" )

