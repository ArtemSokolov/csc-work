library(tidyverse)

## Load QCd expression
X <- read_csv("data/WD-76845-097-ij_subtracted_50_qc.csv",
              col_types=cols( CellID=col_integer() ))

## Columns to fit
vcol <- list(CD3         = "anti_CD3_cytoRingMask",
             CD45RO      = "anti_CD45RO_cytoRingMask",
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

## Marker panel of interest:
## CD3, CD45RO, aSMA, CD163, CD20, CD4, CD45, CD68, CD8a,
## Desmin, Ecad, FOXP3, Keratin, PCNA, PD1, PDL1, Vimentin.
mypanel <- c("CD3", "CD45RO", "aSMA", "CD163", "CD20", "CD4", "CD45", "CD68", "CD8a",
             "Desmin", "Ecad", "FOXP3", "Keratin", "PCNA", "PD1", "PDL1", "Vimentin" )
             
if( !file.exists("models/GMMs.RData") ) {
    ## Fit to the markers of interest only
    GMMs <- naivestates::GMMfit(X, CellID, !!!vcol)
    dir.create("models", showWarnings=FALSE)
    save( GMMs, file="models/GMMs.RData" )

    ## Plot the general overview of the fit quality
    ggf <- naivestates::plotFitOverview( GMMs )
    dir.create("plots/markers", showWarnings=FALSE, recursive=TRUE)
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

## Summarizes cell populations based on combinations of positive and negative markers
popSummary <- function(.df, panel) {
    .df %>% select( CellID, all_of(panel) ) %>%
        filter( if_all(.fns= ~.x != "amb") ) %>%
        group_by( !!!syms(panel) ) %>%
        summarize( nCells=n(), CellIDs=list(CellID), .groups="drop" ) %>%
        arrange( desc(nCells) ) %>%
        mutate( PopIndex = 1:n(), Population = str_c("Pop", PopIndex) )
}

## Determine unique subpopulations in the marker panels of interest
topPop <- function(M, panel, nPop) {
    Y <- popSummary(M, panel) %>% slice( 1:nPop ) %>%
        mutate( Population = factor(Population, Population) )
    
    ## Count the total number and reduce to the top populations
    nTotal <- Y %>% summarize( nTotal = sum(nCells) ) %>%
        mutate( Label = str_c("Total # Cells: ", scales::comma(nTotal)) )

    Y1 <- Y %>% select(Population, `# Cells` = nCells)
    gg1 <- ggplot( Y1, aes(x=Population, y=`# Cells`) ) + theme_minimal() +
        geom_bar( stat='identity', color="darkgray", fill="lightgray" ) +
        scale_y_log10(labels=scales::comma, breaks=c(10,100,1000,10000,100000)) + 
        theme(axis.text.x = element_text(angle=90, vjust=0.5),
              axis.title.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank() ) +
        geom_text( data=nTotal, aes(label=Label), x=Inf, y=Inf, hjust=1, vjust=1 )

    Y2 <- Y %>% select( -nCells, -PopIndex, -CellIDs ) %>%
        pivot_longer( -c(Population), names_to="Marker" )
    gg2 <- ggplot( Y2, aes(x=Population, y=Marker, color=value) ) +
        theme_minimal() + geom_point() +
        scale_color_manual(values=c(pos="black", neg="lightgray"), guide=FALSE ) +
        theme(axis.text.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank() )

    egg::ggarrange( gg1, gg2, ncol=1 )
}

Y <- popSummary(M, mypanel) %>% filter( nCells > 100 ) %>%
    select(Population, nCells, everything())
(ggplot( Y, aes(x=PopIndex, y=nCells) ) +
    geom_point() + theme_bw() +
    scale_y_log10( labels=scales::comma )) %>%
    plotly::ggplotly() %>% htmlwidgets::saveWidget("popsize.html")

popSummary(M, mypanel) %>% filter( nCells > 100 ) %>%
    select(-nCells, -PopIndex) %>% unnest(CellIDs) %>%
    select( CellID=CellIDs, Label=Population, everything() ) %>%
    write_csv( "output/gating.csv" )

gg <- topPop(M, mypanel, 40)
ggsave( "output/gating_summary.png", gg, width=8, height=5 )
