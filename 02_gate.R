library( tidyverse )

P <- read_csv( "output/probs.csv",
              col_types=cols(CellID=col_integer()) )

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

dir.create("plots/02", showWarnings=FALSE, recursive=TRUE)
ggplot( A, aes(x=Threshold, y=FracAmb) ) +
    theme_bw() + geom_line() +
    facet_wrap( ~Marker, nrow=3, ncol=6 ) +
    ylab( "Fraction Ambiguous" ) +
    theme(panel.grid.minor.x = element_blank(),
          strip.background = element_blank()) +
    ggsave( "plots/02/threshold.png", width=8, height=5 )

G80 <- P %>% gate(0.8)
G80 %>% select(-CellID) %>% filter( if_all(.fns= ~.x == "amb") )

G80 %>% write_csv( "output/gating-80.csv" )
