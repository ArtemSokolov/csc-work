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

## Inspect the distribution of intersection values
Z <- XY %>% count( HDBSCAN, Label, name="Count" ) %>%
    arrange( desc(Count) ) %>% filter( Count >= 15 )

Z %>% group_by(HDBSCAN) %>% summarize( MX = max(Count), MN = min(Count) )

ggplot(Z, aes(x=Label, y=Count)) + theme_bw() +
    geom_bar( stat='identity' ) +
    facet_wrap(~HDBSCAN, scales='free') +
    scale_y_log10(labels=scales::comma) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    ggsave("plots/intersect.png", width=8, height=6)

Z1 <- Z %>% 
    mutate( Intersection = str_c(HDBSCAN, "_", Label) ) %>%
    mutate( Intersection = factor(Intersection, Intersection) )
(ggplot( Z1, aes(x=Intersection, y=Count) ) +
 theme_bw() + geom_point() +
 scale_y_log10(labels=scales::comma) +
 theme(axis.text.x = element_text(angle=90, vjust=0.5))) %>%
    plotly::ggplotly() %>% htmlwidgets::saveWidget("plots/intersect.html")

