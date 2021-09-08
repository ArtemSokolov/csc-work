library( tidyverse )
library( seriation )

nPop <- 40

## Load cluster assignments and gating results
X <- read_csv( "data/hdbscan.csv", col_types=cols() ) %>%
    rename( HDBSCAN = hdbscan )
Y <- read_csv( "output/gating.csv", col_types=cols() ) %>%
    select( CellID, Label ) %>%
    filter( Label %in% str_c("Pop",1:nPop) )

## Given a set of labeled cells, counts the number of clusters that don't have any labels
## dfClus - data frame mappping CellID to cluster index
## dfLbl  - data frame mapping CellID to population labels
## cl     - column in .dfClus that contains cluster indices
countUnlab <- function( dfClus, dfLbl, cl ) {
    left_join( dfClus, dfLbl, by="CellID" ) %>%
        count( {{cl}}, Label, name="Count" ) %>%
        count( {{cl}}, name="nLabels" ) %>%
        filter( nLabels == 1 ) %>% nrow()
}

## countUnlab( X, Y, HDBSCAN )

## Compute intersections of all clusters by labels
XY <- left_join(X, Y, by="CellID") %>%
    mutate_at( "Label", replace_na, "Other" ) %>%
    count( HDBSCAN, Label, name="Count" ) %>%
    pivot_wider( names_from=Label, values_from=Count, values_fill=0 ) %>%
    pivot_longer( -HDBSCAN, names_to="Label", values_to="Count" )

## Plot raw counts
(function() {
    XYp <- XY %>%
        mutate(Txt = ifelse(Count == 0, "", scales::label_number_si()(Count)),
               HDBSCAN = factor(HDBSCAN, unique(HDBSCAN)),
               Label = factor(Label, c(str_c("Pop",1:nPop),"Other")))

    ggplot(XYp, aes(y=Label, x=HDBSCAN)) + theme_bw() +
        geom_tile( color="gray", fill="white" ) +
        geom_text(aes(label=Txt)) +
        ggsave( "plots/hdbscan/raw.png", width=10, height=10 ) +
        theme( panel.grid.major=element_blank(),
              panel.grid.minor=element_blank() )
})()

## Normalize counts by a given column
fnorm <- function(.df, .col) {
    .df %>% filter( HDBSCAN != -1, Label != "Other" ) %>%
        group_by({{.col}}) %>%
        mutate(Frac = Count / sum(Count),
               Frac = ifelse(is.nan(Frac),0,Frac)) %>%
        ungroup() %>% select(-Count)
}

## Compute optimal leaf order for a given column .col
##  by clustering on the distances across features in column .feat
foptlo <- function(.df, .col, .feat) {
    DM <- .df %>% pivot_wider( names_from=all_of(.feat), values_from="Frac" ) %>%
        as.data.frame() %>% column_to_rownames(.col) %>% dist()
    hclust(DM) %>% reorder(DM) %>% dendextend::order.hclust() %>% labels(DM)[.]
}

## Plot normalized and reordered rows and columns
fplot <- function(.df, clOrd, lblOrd, fnOut) {
    .dfp <- .df %>%
        mutate(Txt = scales::label_number(accuracy=0.01)(Frac),
               Txt = ifelse(Txt == "0.00", "", Txt),
               HDBSCAN = factor(HDBSCAN, clOrd),
               Label = factor(Label, lblOrd))

    pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))[4:7]
    ggplot(.dfp, aes(y=Label, x=HDBSCAN)) + theme_bw() +
        geom_tile( aes(fill=Frac), color="gray" ) +
        geom_text(aes(label=Txt)) +
        scale_fill_gradientn( colors=pal, limits=c(0,1), guide=FALSE ) +
        theme( panel.grid.major=element_blank(),
              panel.grid.minor=element_blank() ) +
        ggsave( fnOut, width=10, height=10 )
}

## Normalize and reorder by cluster
ZC <- fnorm(XY, HDBSCAN)
fplot( ZC, foptlo(ZC, "HDBSCAN", "Label"),
      str_c("Pop", 1:nPop), "plots/hdbscan/clusters.png" )

## Normalize and reorder by population
ZP <- fnorm(XY, Label)
fplot( ZP, 0:max(X$HDBSCAN), foptlo(ZP, "Label", "HDBSCAN"),
      "plots/hdbscan/pops.png" )
