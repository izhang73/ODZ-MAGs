library(UpSetR)
library(xlsx)
library(ComplexUpset)

##set working directory
setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs") 

AMAL_upsetR <- read_excel("Bess-Xin-Amal-MAGs.xlsx", 
                                sheet="Presence-absence")

AMAL_Arabian <- read_excel("Bess-Xin-Amal-MAGs.xlsx", 
                          sheet="Presence-absence-Arabian")

AMAL_ETNP<- read_excel("Bess-Xin-Amal-MAGs.xlsx", 
                           sheet="Presence-absence-ETNP")

AMAL_ETNP_oxy<- read_excel("Bess-Xin-Amal-MAGs.xlsx", 
                           sheet="Presence-absence-ETNP-oxycline")

AMAL_ETNP_ODZ<- read_excel("Bess-Xin-Amal-MAGs.xlsx", 
                           sheet="Presence-absence-ETNP-ODZ")

AMAL.df <- as.data.frame(AMAL_upsetR)
AMAL.df <- AMAL.df[,-1]
AMAL_no_nap.df <- AMAL.df[,-1]

AMAL_Arabian.df <- as.data.frame(AMAL_Arabian)
AMAL_Arabian.df <- AMAL_Arabian.df[,-1]

AMAL_ETNP.df <- as.data.frame(AMAL_ETNP)
AMAL_ETNP.df <- AMAL_ETNP.df[,-1]

AMAL_oxycline.df <- as.data.frame(AMAL_ETNP_oxy)
AMAL_oxycline.df <- AMAL_oxycline.df[,-1]

AMAL_ODZ.df <- as.data.frame(AMAL_ETNP_ODZ)
AMAL_ODZ.df <- AMAL_ODZ.df[,-1]

all_sets <- upset(AMAL.df,
      sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
      order.by="freq", main.bar.color = "maroon", sets.bar.color ="blue", 
      nintersects = 40, point.size = 2.5, keep.order=TRUE)

genes = colnames(AMAL.df)

genes_no_nap <- genes[-1]

complex_all_sets <- upset(AMAL.df,genes, 
      sort_intersections_by=c('degree'),
      sort_sets=FALSE,
      base_annotations=list(
        'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
        'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
        )
      )

complex_all_sets_inclusive <- upset(AMAL.df,genes, 
                          sort_intersections_by=c('degree'),
                          sort_sets=FALSE,
                          base_annotations=list(
                            'Intersection ratio'=intersection_ratio(mode = "intersect", text_mapping=aes(label=!!upset_text_percentage(digits = 0, sep = "", mode = "intersect"))),
                            'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                          )
)


complex_ETNP_sets <- upset(AMAL_ETNP.df,genes, 
                          sort_intersections_by=c('degree'),
                          sort_sets=FALSE,
                          base_annotations=list(
                            'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                            'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                          )
)


complex_ETNP_sets_inclusive <- upset(AMAL_ETNP.df,genes, 
                                    sort_intersections_by=c('degree'),
                                    sort_sets=FALSE,
                                    base_annotations=list(
                                      'Intersection ratio'=intersection_ratio(mode = "intersect", text_mapping=aes(label=!!upset_text_percentage(digits = 0, sep = "", mode = "intersect"))),
                                      'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                    )
)

complex_Arabian_sets <- upset(AMAL_Arabian.df,genes, 
                           sort_intersections_by=c('degree'),
                           sort_sets=FALSE,
                           base_annotations=list(
                             'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                             'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                           )
)

complex_Arabian_sets_inclusive <- upset(AMAL_Arabian.df,genes, 
                                     sort_intersections_by=c('degree'),
                                     sort_sets=FALSE,
                                     base_annotations=list(
                                       'Intersection ratio'=intersection_ratio(mode = "intersect", text_mapping=aes(label=!!upset_text_percentage(digits = 0, sep = "", mode = "intersect"))),
                                       'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                     )
)

complex_ETNP_oxy_sets <- upset(AMAL_oxycline.df,genes, 
                              sort_intersections_by=c('degree'),
                              sort_sets=FALSE,
                              base_annotations=list(
                                'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                                'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                              )
)

complex_ETNP_oxy_sets_inclusive <- upset(AMAL_oxycline.df,genes, 
                                        sort_intersections_by=c('degree'),
                                        sort_sets=FALSE,
                                        base_annotations=list(
                                          'Intersection ratio'=intersection_ratio(mode = "intersect", text_mapping=aes(label=!!upset_text_percentage(digits = 0, sep = "", mode = "intersect"))),
                                          'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                        )
)

complex_ETNP_ODZ_sets <- upset(AMAL_ODZ.df,genes, 
                               sort_intersections_by=c('degree'),
                               sort_sets=FALSE,
                               base_annotations=list(
                                 'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                                 'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                               )
)

complex_ETNP_ODZ_sets_inclusive <- upset(AMAL_ODZ.df,genes, 
                                         sort_intersections_by=c('degree'),
                                         sort_sets=FALSE,
                                         base_annotations=list(
                                           'Intersection ratio'=intersection_ratio(mode = "intersect", text_mapping=aes(label=!!upset_text_percentage(digits = 0, sep = "", mode = "intersect"))),
                                           'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                         )
)


complex_all_sets_no_nap <- upset(AMAL.df,genes_no_nap, 
                          sort_intersections_by=c('degree'),
                          sort_sets=FALSE,
                          base_annotations=list(
                            'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                            'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                          )
)

complex_ETNP_sets_no_nap <- upset(AMAL_ETNP.df,genes_no_nap, 
                                 sort_intersections_by=c('degree'),
                                 sort_sets=FALSE,
                                 base_annotations=list(
                                   'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                                   'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                 )
)

complex_Arabian_sets_no_nap <- upset(AMAL_Arabian.df,genes_no_nap, 
                                  sort_intersections_by=c('degree'),
                                  sort_sets=FALSE,
                                  base_annotations=list(
                                    'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                                    'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                  )
)

complex_ETNP_oxy_sets_no_nap <- upset(AMAL_oxycline.df,genes_no_nap, 
                                     sort_intersections_by=c('degree'),
                                     sort_sets=FALSE,
                                     base_annotations=list(
                                       'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                                       'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                     )
)

complex_ETNP_ODZ_sets_no_nap <- upset(AMAL_ODZ.df,genes_no_nap, 
                                      sort_intersections_by=c('degree'),
                                      sort_sets=FALSE,
                                      base_annotations=list(
                                        'Intersection ratio'=intersection_ratio(mode = "distinct", text_mapping=aes(label=!!upset_text_percentage())),
                                        'Size'=(intersection_size(counts=FALSE, mode='inclusive_intersection'))
                                      )
)

all_sets_ETNP <- upset(AMAL_ETNP.df,
                           sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
                           order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                           nintersects = 40, point.size = 2.5, keep.order=TRUE)

no_nap_sets_ETNP <- upset(AMAL_ETNP.df,
                       sets = c("narG", "nirK", "nirS", "norB", "nosZ"),
                       order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                       nintersects = 40, point.size = 2.5, keep.order=TRUE)


all_sets_ETNP_oxy <- upset(AMAL_oxycline.df,
                  sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
                  order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                  nintersects = 40, point.size = 2.5, keep.order=TRUE)

no_nap_sets_ETNP_oxy <- upset(AMAL_oxycline.df,
                           sets = c("narG", "nirK", "nirS", "norB", "nosZ"),
                           order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                           nintersects = 40, point.size = 2.5, keep.order=TRUE)

all_sets_ETNP_ODZ <- upset(AMAL_ODZ.df,
                           sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
                           order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                           nintersects = 40, point.size = 2.5, keep.order=TRUE)

no_nap_sets_ETNP_ODZ <- upset(AMAL_ODZ.df,
                              sets = c("narG", "nirK", "nirS", "norB", "nosZ"),
                              order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                              nintersects = 40, point.size = 2.5, keep.order=TRUE)

all_sets_Arabian <- upset(AMAL_Arabian.df,
                  sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
                  order.by="freq", main.bar.color = "maroon", sets.bar.color ="blue", 
                  nintersects = 40, point.size = 2.5, keep.order=TRUE)

all_sets_Arabian_degree <- upset(AMAL_Arabian.df,
                          sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
                          order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                          nintersects = 40, point.size = 2.5, keep.order=TRUE)

no_nap_sets_Arabian <- upset(AMAL_Arabian.df,
                     sets = c("narG", "nirK", "nirS", "norB", "nosZ"),
                     order.by = "degree", main.bar.color = "maroon", sets.bar.color ="blue",
                     nintersects = 40, point.size = 2.5, keep.order = TRUE)

all_sets_degree <- upset(AMAL.df,
                  sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
                  order.by="degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                  nintersects = 40, point.size = 2.5, keep.order=TRUE)

all_sets_group <- upset(AMAL.df,
                         sets = c("napA", "narG", "nirK", "nirS", "norB", "nosZ"),
                         group.by = "degree", main.bar.color = "maroon", sets.bar.color ="blue", 
                         nintersects = 40, point.size = 2.5, keep.order=TRUE)

no_nap_sets <- upset(AMAL.df,
                  sets = c("narG", "nirK", "nirS", "norB", "nosZ"),
                  order.by = "freq", main.bar.color = "maroon", sets.bar.color ="blue",
                  nintersects = 40, point.size = 2.5, keep.order = TRUE)

no_nap_sets_degree <- upset(AMAL.df,
                     sets = c("narG", "nirK", "nirS", "norB", "nosZ"),
                     order.by = "degree", main.bar.color = "maroon", sets.bar.color ="blue",
                     nintersects = 40, point.size = 2.5, keep.order = TRUE)

no_nap_sets_group <- upset(AMAL.df,
                            sets = c("narG", "nirK", "nirS", "norB", "nosZ"),
                            group.by = "degree", main.bar.color = "maroon", sets.bar.color ="blue",
                            nintersects = 40, point.size = 2.5, keep.order = TRUE)
