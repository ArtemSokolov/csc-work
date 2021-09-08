library( tidyverse )

## Load cluster assignments and gating results
X <- read_csv( "data/hdbscan.csv", col_types=cols() ) %>%
    rename( HDBSCAN = hdbscan ) %>%
    filter( HDBSCAN != -1 )
Y <- read_csv( "output/gating.csv", col_types=cols() )
XY <- inner_join(X, Y, by="CellID")

## Isolate intersections that pass the cell count filter
Z <- XY %>% group_by(HDBSCAN, Label) %>%
    summarize( Count = n(), across(CD3:Vimentin, unique), .groups='drop' ) %>%
    group_by(HDBSCAN) %>% mutate( Pct = Count / sum(Count) ) %>% ungroup() %>%
    filter( Count >= 10, Pct >= 0.05 )

finspect <- function() {
    Z %>% group_by(HDBSCAN) %>% summarize( MX = max(Count), MN = min(Count) )
    
    Z1 <- Z %>% arrange(desc(Count)) %>%
        mutate( Intersection = str_c(HDBSCAN, "_", Label) ) %>%
        mutate( Intersection = factor(Intersection, Intersection) )
    (ggplot( Z1, aes(x=Intersection, y=Count) ) +
     theme_bw() + geom_point() +
     scale_y_log10(labels=scales::comma) +
     theme(axis.text.x = element_text(angle=90, vjust=0.5))) %>%
        plotly::ggplotly() %>% htmlwidgets::saveWidget("intersect.html")
}

## Common theme elements
mytheme <- function()
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.minor = element_blank())

fplot <- function(.df, clusRange) {
    ZZ <- .df %>% filter(HDBSCAN %in% clusRange) %>%
        arrange( HDBSCAN, desc(Count) ) %>%
        mutate( Slice = str_c(HDBSCAN, '_', Label) ) %>%
        mutate( Slice = factor(Slice,Slice) )

    ZZ1 <- ZZ %>% select(Slice, HDBSCAN, `# Cells` = Count)
    gg1 <- ggplot(ZZ1, aes(x=Slice, y=`# Cells`)) + theme_minimal() +
        geom_bar( stat='identity', color='darkgray', fill='lightgray' ) +
        scale_y_log10(labels=scales::comma, breaks=c(10,100,1000,10000,100000)) +
        mytheme() + facet_grid( ~HDBSCAN, scales='free', space='free' ) +
        theme(panel.grid.major.x = element_blank(),
              plot.margin=unit(c(0,0.5,0,0.5),unit="cm"))

    ZZ2 <- ZZ %>% select(-Label, -Count, -Pct) %>%
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

ggsave( "plots/intersect1.png", fplot( Z, 0:9 ), width=8, height=3 )
ggsave( "plots/intersect2.png", fplot( Z, 10:18 ), width=8, height=3 )

