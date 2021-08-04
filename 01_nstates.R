library(tidyverse)

## Load QCd expression
X <- read_csv("data/WD-76845-097-ij_subtracted_50_qc.csv",
              col_types=cols( CellID=col_integer() ))

## Columns to fit
vcol <- list(anti_CD3    = "anti_CD3_cytoRingMask",
             anti_CD45RO = "anti_CD45RO_cytoRingMask",
             aSMA        = "aSMA_660_cellRingMask",
             Keratin     = "Keratin_570_cellRingMask",
             CD4         = "CD4_488_cytoRingMask",
             CD45        = "CD45_PE_cytoRingMask",
             PD1         = "PD1_647_cytoRingMask",
             CD20        = "CD20_488_cytoRingMask",
             CD68        = "CD68_555_cellRingMask",
             CD8a        = "CD8a_660_cytoRingMask",
             CD163       = "CD163_488_cellRingMask",
             FOXP3       = "FOXP3_570_nucleiRingMask",
             PDL1        = "PDL1_647_cytoRingMask",
             Ecad        = "Ecad_488_cellRingMask",
             Vimentin    = "Vimentin_555_cellRingMask",
             CDX2        = "CDX2_647_cellRingMask",
             LaminABC    = "LaminABC_488_nucleiRingMask",
             Desmin      = "Desmin_555_cellRingMask",
             CD31        = "CD31_647_nucleiRingMask",
             PCNA        = "PCNA_488_nucleiRingMask",
             CollagenIV  = "CollagenIV_647_cellRingMask")

## Marker panels of interest
panel1 <- c("anti_CD3", "anti_CD45RO", "CD4", "CD45", "PD1", "CD20", "CD68", "CD8a", "CD163", "FOXP3")
panel2 <- c(panel1, "CD31", "aSMA", "Vimentin", "Keratin", "PDL1")

if( !file.exists("models/GMMs.RData") ) {
    ## Fit to the markers of interest only
    GMMs <- naivestates::GMMfit(X, CellID, !!!vcol)
    save( GMMs, file="models/GMMs.RData" )

    ## Plot the general overview of the fit quality
    ggf <- naivestates::plotFitOverview( GMMs )
    ggsave( "plots/overview.png", ggf, width=8, height=8 )
    walk(names(vcol), ~ggsave(
                         str_c("plots/markers/", .x, ".png"),
                         naivestates::plotMarker(GMMs, .x)
                     ))
} else {
    load("models/GMMs.RData")
}

## Extract probabilities of expression
P <- naivestates::GMMreshape(GMMs)

## A marker is considered expressed if p >= 0.8
##   not expressed if p <= 0.2
##   and ambiguous otherwise
M <- P %>% mutate( across(-CellID, ~case_when(
                                     .x >= 0.8 ~ "pos",
                                     .x <= 0.2 ~ "neg",
                                     TRUE ~ "amb")) )

## Determine unique subpopulations in the marker panels of interest
topPop <- function( panel, nPop ) {
    Y <- M %>% select( all_of(panel) ) %>%
        filter( if_all(.fns= ~.x != "amb") ) %>%
        group_by( !!!syms(panel) ) %>% tally() %>%
        ungroup() %>% arrange( desc(n) ) %>%
        mutate( Population = 1:n() )

    ## Count the total number and reduce to the top populations
    nTotal <- Y %>% summarize( nTotal = sum(n) ) %>%
        mutate( Label = str_c("Total # Cells: ", scales::comma(nTotal)) )
    Y <- Y %>% slice( 1:nPop )

    Y1 <- Y %>% select(Population, `# Cells` = n)
    gg1 <- ggplot( Y1, aes(x=Population, y=`# Cells`) ) + theme_minimal() +
        geom_bar( stat='identity', color="darkgray", fill="lightgray" ) +
        scale_x_continuous(limits=c(0,(nPop+1)), expand=c(0,0)) +
        scale_y_log10(labels=scales::comma, breaks=c(10,100,1000,10000,100000)) + 
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank() ) +
        geom_text( data=nTotal, aes(label=Label), x=Inf, y=Inf, hjust=1, vjust=1 )

    Y2 <- Y %>% select( -n ) %>%
        pivot_longer( -c(Population), names_to="Marker" )
    gg2 <- ggplot( Y2, aes(x=Population, y=Marker, color=value) ) +
        theme_minimal() + geom_point() +
        scale_color_manual(values=c(pos="black", neg="lightgray"), guide=FALSE ) +
        scale_x_continuous(limits=c(0,(nPop+1)), expand=c(0,0)) +
        theme(axis.text.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank() )

    egg::ggarrange( gg1, gg2, ncol=1 )
}

gg1 <- topPop( panel1, 20 )
ggsave( "immune.png", gg1, width=8, height=5 )

gg2 <- topPop( panel2, 20 )
ggsave( "all.png", gg2, width=8, height=5 )

