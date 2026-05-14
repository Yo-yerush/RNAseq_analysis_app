Revigo_plots <- function(GO_df_up, GO_df_down,
                         GO_type, #direction,
                         show_labels=TRUE,
                         specific_terms_lables_up=NULL, # specific by term column
                         specific_terms_lables_down=NULL, # specific by term column
                         specific_parents_lables_up=NULL, # specific by parents column
                         specific_parents_lables_down=NULL, # specific by parents column
                         specific_parents_lables_show=NULL # if null, will take the default terms label
                         ) {
  
  library(rrvgo)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # prepare revigo data frame
  revigo_df_fun <- function(GO_df) {
    
    if (!is.numeric(GO_df$pValue)) {
      GO_df$pValue = gsub("< 1e-30", 1e-30, GO_df$pValue) %>% as.numeric()
    }
    
    ### rrvgo pipeline
    simMatrix <- calculateSimMatrix(GO_df$GO.ID,
                                    orgdb="org.At.tair.db",
                                    ont=GO_type,
                                    method="Rel")
    
    scores <- setNames(-log10(GO_df$pValue), GO_df$GO.ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.At.tair.db")
    
    tryCatch({
      x <<- scatterPlot(simMatrix, reducedTerms, algorithm = "umap")$data
    }, error = function(e){
      message("using 'PCA' instead of 'UMAP'")
      x <<- scatterPlot(simMatrix, reducedTerms, algorithm = "pca")$data
    })
    
    xx = data.frame(term_ID = row.names(x),
                           description = as.character(x$term),
                           plot_X = as.numeric(x$V1),
                           plot_Y = as.numeric(x$V2),
                           value = as.numeric(x$score),
                           parent = as.character(x$parent),
                           parentTerm = as.character(x$parentTerm))
    return(list(plot_df = xx,
                reducedTerms = reducedTerms))
  }

  revigo_up_0 = revigo_df_fun(GO_df_up)
  revigo_up = revigo_up_0$plot_df
  revigo_up$direction = "up"
  
  revigo_down_0 = revigo_df_fun(GO_df_down)
  revigo_down = revigo_down_0$plot_df
  revigo_down$direction = "down"
  
  
  # bind 'up and 'down' data frames
  revigo_df = rbind(revigo_up,revigo_down) %>% arrange(parentTerm)# %>% arrange(desc(direction))
  
  
  # set colors for different parent terms
  n.unique_parents = length(unique(revigo_df$parentTerm))
  if (n.unique_parents > 17) {
    color_palette = rainbow(n.unique_parents)
  } else {
    color_palette <- suppressWarnings(brewer.pal(n = n.unique_parents, "Set1"))
    if (length(unique(revigo_df$parentTerm)) > 9) {
      color_palette = suppressWarnings(c(color_palette, brewer.pal(n = n.unique_parents-9, "Set2")))
    }
  }
  
  # parents and color index
  parent_col = data.frame(parentTerm = unique(revigo_df$parentTerm), hex_col = color_palette)
  revigo_df = merge(revigo_df, parent_col, by = "parentTerm")
  revigo_df_up = revigo_df[revigo_df$direction == "up",]
  revigo_df_down = revigo_df[revigo_df$direction == "down",]
  
  
  # max n.char for margin argument
  terms_character = paste(specific_terms_lables_up,
                          specific_terms_lables_down,
                          specific_parents_lables_up,
                          specific_parents_lables_down, sep = "|")
  if (length(terms_character) != 0) {
    lables_nchar = max(nchar(strsplit(terms_character, split = "\\|")[[1]]))
    lables_nchar_index = lables_nchar*0.15
  } else {
    lables_nchar_index = ifelse(show_labels, 6, 0.1)
  }

  
  ##############################
  ##############################
  ##############################
  
  scatterPlot_function <- function(df, specific_terms_lables, specific_parents_lables, margin_index = lables_nchar_index) {
    
    one.x_range = (max(df$plot_X) - min(df$plot_X))
    one.y_range = (max(df$plot_Y) - min(df$plot_Y))
    
    #### labels - filter terms - 4 options
    # by specific terms
    if (!is.null(specific_terms_lables)) {
      filtered_data = df[grep(specific_terms_lables, df$description),]%>% # keep specific terms
        dplyr::distinct(description, .keep_all = TRUE) # remove duplicate terms
      filtered_data$parentTerm = filtered_data$description
      
    } else if (!is.null(specific_parents_lables)) {
      filtered_data = df[grep(specific_parents_lables, df$parentTerm),]%>% # keep specific terms
        dplyr::distinct(parentTerm, .keep_all = TRUE) # remove duplicate terms
      
      # change to costum label
      if (!is.null(specific_parents_lables_show)) {
        # Replace values in filtered_data$parentTerm
        lookup <- setNames(specific_parents_lables_show$b, specific_parents_lables_show$a)
        filtered_data$parentTerm <- ifelse(filtered_data$parentTerm %in% names(lookup),
                                           lookup[filtered_data$parentTerm],
                                           filtered_data$parentTerm)
      }
    } else {
      
      
      n.n = 8
      n.unique = length(unique(df$parentTerm))
      # less than n.n
      if (n.unique > n.n) {
        
        filtered_data <- df %>%
          dplyr::group_by(parentTerm) %>%
          dplyr::filter(n() >= 2) %>% # filter term that appears less than 3 times
          dplyr::ungroup() %>%
          dplyr::distinct(parentTerm, .keep_all = TRUE) %>% # remove duplicate terms
          as.data.frame()
        
        # more than n.n
        if (nrow(filtered_data) < n.n) {
          filtered_data = df[!duplicated(df$parentTerm),]
        } else {
          filtered_data <- df %>%
            dplyr::group_by(parentTerm) %>%
            dplyr::summarize(count = n()) %>% # keep top n.n terms
            dplyr::arrange(desc(count)) %>%
            dplyr::slice_max(order_by = count, n = n.n) %>%
            dplyr::slice_head(n = n.n) %>% 
            dplyr::select(parentTerm) %>%
            dplyr::inner_join(df, by = "parentTerm") %>%
            dplyr::ungroup()  %>%
            dplyr::distinct(parentTerm, .keep_all = TRUE) %>%
            as.data.frame() 
        }
        
      } else {
        filtered_data = df[!duplicated(df$parentTerm),] # keep all terms
      }
    }
    
    
    scatterPlot_vec = df %>% ggplot2::ggplot(ggplot2::aes(x = plot_X, y = plot_Y, color = parentTerm)) + 
      ggplot2::geom_point(ggplot2::aes(size = value), alpha = 0.5, show.legend = FALSE) + 
      
      ggplot2::scale_color_manual(values = unique(df$hex_col), guide = "none") + 
      ggplot2::geom_point(ggplot2::aes(plot_X, plot_Y, size = value), 
                          shape = 21, fill = "transparent", colour = I(ggplot2::alpha("black", 0.6)),
                          show.legend = FALSE) + 
      
      ggplot2::theme_classic() + 
      ggplot2::geom_hline(yintercept=0, colour="grey50", linetype = "dashed") + 
      ggplot2::geom_vline(xintercept = 0, colour="grey50", linetype = "dashed") +
      
      ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA, size=0.75),
                     axis.text.y = element_text(size = 9),
                     axis.text.x = element_text(size = 9),
                     axis.line = element_blank(),
                     plot.margin =  unit(c(0.1, margin_index, 0.1, 0.1), "cm")) +
      
      ggrepel::geom_text_repel(ggplot2::aes(label = parentTerm), 
                               data = filtered_data,
                               direction = "y",
                               hjust = 0,
                               #alpha = 0.65,
                               #min.segment.length = 0,
                               #family = 'Times',
                               fontface = 'bold',
                               #box.padding = 1,
                               #box.padding = 0.35, 
                               #point.padding = 0.5,
                               #point.padding = NA,
                               #label.padding = .2,
                               size = 3,
                               #max.time = 10,
                               #max.iter = 1000000,
                               #force = 5,
                               colour="black",
                               xlim = c(max(df$plot_X)+one.x_range/7.5, Inf)) + 
      
      ggplot2::labs (y = "semantic space y", x = "semantic space x") +#, size = "-log10(p)") + 
      #ggplot2::guides(color = guide_legend(override.aes = list(size = 4, shape = 21, colour = "black", fill = unique(df$hex_col)))) +
      #ggplot2::theme(legend.key = ggplot2::element_blank(), legend.box = "horizontal") +
      ggplot2::coord_cartesian(clip = "off",
                               xlim = c(min(df$plot_X)-one.x_range/10,max(df$plot_X)+one.x_range/10),
                               ylim = c(min(df$plot_Y)-one.y_range/10,max(df$plot_Y)+one.y_range/10))
    
    return(scatterPlot_vec)
  }
  scatterPlot_up = scatterPlot_function(revigo_df_up, specific_terms_lables_up, specific_parents_lables_up)
  scatterPlot_down = scatterPlot_function(revigo_df_down, specific_terms_lables_down, specific_parents_lables_down)
  
 
  return(list(scatterPlot_up = scatterPlot_up,
              scatterPlot_down = scatterPlot_down,
              scatterPlot_legend = parent_col,
              treemapPlot_up = revigo_up_0$reducedTerms,
              treemapPlot_down = revigo_down_0$reducedTerms,
              revigo_df = revigo_df))
}