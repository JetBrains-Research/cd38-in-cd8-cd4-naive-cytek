Downsampling_FlowSet <- function(x, samplesize , replace=TRUE, prob=NULL, seed_number=seed_number){
  if(missing(samplesize)) samplesize <- min(flowCore::fsApply(x,nrow))
  flowCore::fsApply(x, function(ff){
    set.seed(seed_number)
    i <- sample(nrow(ff), size = samplesize, replace=replace, prob)
    ff[i,]
  })
}

plotDR_2 <- function(x, dr = NULL, 
                     color_by = "condition", facet_by = NULL, ncol = NULL,
                     assay = "exprs", scale = TRUE, q = 0.01, dims = c(1, 2),
                     k_pal = CATALYST:::.cluster_cols, 
                     a_pal = hcl.colors(10, "Viridis"), color_range = NULL, used_in_clustering) {
  
  # check validity of input arguments
  stopifnot(
    is(x, "SingleCellExperiment"),
    .check_assay(x, assay),
    length(reducedDims(x)) != 0,
    is.logical(scale), length(scale) == 1,
    is.numeric(q), length(q) == 1, q >= 0, q < 0.5)
  .check_pal(a_pal)
  .check_cd_factor(x, facet_by)
  
  if (!is.null(ncol)) 
    stopifnot(is.numeric(ncol), length(ncol) == 1, ncol %% 1 == 0)
  
  if (is.null(dr)) {
    dr <- reducedDimNames(x)[1]
  } else {
    stopifnot(
      is.character(dr), length(dr) == 1, 
      dr %in% reducedDimNames(x))
  }
  stopifnot(is.numeric(dims), length(dims) == 2, 
            dims %in% seq_len(ncol(reducedDim(x, dr))))
  
  if (!all(color_by %in% rownames(x))) {
    stopifnot(length(color_by) == 1)
    if (!color_by %in% names(colData(x))) {
      .check_sce(x, TRUE)
      .check_pal(k_pal)
      .check_k(x, color_by)
      kids <- cluster_ids(x, color_by)
      nk <- nlevels(kids)
      if (length(k_pal) < nk)
        k_pal <- colorRampPalette(k_pal)(nk)
    } else kids <- NULL
  }
  
  # construct data.frame of reduced dimensions & relevant cell metadata
  xy <- reducedDim(x, dr)[, dims]
  colnames(xy) <- c("x", "y")
  df <- data.frame(colData(x), xy, check.names = FALSE)
  if (all(color_by %in% rownames(x))) {
    es <- as.matrix(assay(x, assay))
    es <- es[color_by, , drop = FALSE]
    if (scale) 
      es <- .scale_exprs(es, 1, q)
    df <- reshape2::melt(
      cbind(df, t(es)), 
      id.vars = colnames(df))
    l <- switch(assay, exprs = "expression", assay)
    l <- paste0("scaled\n"[scale], l)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    thm <- guide <- NULL
    color_by <- "value"
    facet <- NULL
  } else if (is.numeric(df[[color_by]])) {
    if (scale) {
      vs <- as.matrix(df[[color_by]])
      df[[color_by]] <- .scale_exprs(vs, 2, q)
    }
    l <- paste0("scaled\n"[scale], color_by)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    color_by <- sprintf("`%s`", color_by)
    facet <- thm <- guide <- NULL
  } else {
    facet <- NULL
    if (!is.null(kids)) {
      df[[color_by]] <- kids
      scale <- scale_color_manual(values = k_pal)
    } else scale <- NULL
    n <- nlevels(droplevels(factor(df[[color_by]])))
    guide <- guides(col = guide_legend(
      ncol = ifelse(n > 12, 2, 1),
      override.aes = list(alpha = 1, size = 3))) 
    thm <- theme(legend.key.height = unit(0.8, "lines"))
  }
  
  # set axes equal for linear dimension reductions
  if (dr %in% c("PCA", "MDS")) {
    asp <- coord_equal()
  } else asp <- NULL
  
  # get axes labels
  if (dr == "PCA") {
    labs <- paste0("PC", dims)
  } else labs <- paste0(dr, "_", dims)
  
  # remove cells for which no reduced dimensions are available
  df <- df[!(is.na(df$x) | is.na(df$y)), ]
  
  colors <- c('navyblue', 'darkmagenta', 'darkorange1')
  # %>% arrange(df[[color_by]])
  p <- ggplot(df, aes_string("x", "y", color = color_by)) +
    geom_point(size = 0.4, alpha = 0.8) + 
    labs(x = labs[1], y = labs[2]) +
    scale_color_gradientn(colours=a_pal, limits=color_range) +
    #+ scale
    facet  + guide + asp + ggtitle(paste(levels(df$variable), used_in_clustering)) +
    theme_minimal() + thm + theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      aspect.ratio = if (is.null(asp)) 1 else NULL)
  
  if (is.null(facet_by)) 
    return(p)
  
  if (is.null(facet)) {
    p + facet_wrap(facet_by, ncol = ncol)
  } else {
    if (nlevels(df$variable) == 1) {
      p + facet_wrap(facet_by, ncol = ncol) + 
        ggtitle(paste(levels(df$variable), used_in_clustering))
    } else {
      fs <- c("variable", facet_by)
      ns <- vapply(df[fs], nlevels, numeric(1))
      if (ns[2] > ns[1]) fs <- rev(fs)
      p + facet_grid(reformulate(fs[1], fs[2]))
    }
  }
}
.check_assay <- function(x, y) {
  stopifnot(
    length(y) == 1, 
    is.character(y),
    sum(y == assayNames(x)) == 1)
  return(TRUE)
}
.check_pal <- function(x, n = 2) {
  if (is.null(x)) 
    return(TRUE)
  stopifnot(
    length(x) >= n,
    is.character(x))
  foo <- tryCatch(col2rgb(x),
                  error = function(e) {})
  if (is.null(foo)) {
    arg_nm <- deparse(substitute(x))
    stop(sprintf("'%s' is invalid.", arg_nm))
  }
  return(TRUE)
}
.check_cd_factor <- function(x, y, n = 1) {
  if (is.null(y))
    return(TRUE)
  if (!is.null(n))
    stopifnot(length(y) == n)
  stopifnot(
    is.character(y), 
    all(y %in% names(colData(x))),
    !vapply(colData(x)[y], is.numeric, logical(1)))
  return(TRUE)
}
.check_sce <- function(x, y = FALSE) {
  stopifnot(
    is(x, "SingleCellExperiment"), 
    !is.null(x$sample_id))
  if (y) 
    stopifnot(
      !is.null(x$cluster_id),
      !is.null(metadata(x)$cluster_codes))
}
.check_k <- function(x, k) {
  kids <- names(cluster_codes(x))
  if (is.null(k)) return(kids[1])
  stopifnot(length(k) == 1, is.character(k))
  if (!k %in% kids)
    stop("Clustering ", dQuote(k), 
         " doesnt't exist; valid are",
         " 'names(cluster_codes(x))'.")
  return(k)
}

gauss <- function(x, mu, sigma, A) {
  return(A * dnorm(x, mu, sigma))
}

bimodal <- function(x, mu1, sigma1, A1, mu2, sigma2, A2) {
  return(gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2))
}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

fit_gaussians_and_find_thresholds <- function(df, feature, title_info = "", xmin = -Inf, xmax = Inf) {
  samples_thresholds <- c()
  
  normaldens_first <- data.frame()
  normaldens_second <- data.frame()
  annotation_df <- data.frame()
  grid <- seq(min(df[,feature]), max(df[,feature]), length = 100)
  sample_id <- unique(df$sample_id)[13]
  for(sample_id in unique(df$sample_id)) {
    wait <- df[df$sample_id == sample_id, feature]
    
    wait.kmeans <- kmeans(wait, 2)
    wait.kmeans.cluster <- wait.kmeans$cluster
    
    wait.df <- data_frame(x = wait, cluster = wait.kmeans.cluster)
    
    wait.summary.df <- wait.df %>%
      group_by(cluster) %>%
      dplyr::summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    
    #wait.summary.df %>%
    #  select(cluster, mu, variance, std)
    
    wait.summary.df <- wait.summary.df %>%
      mutate(alpha = size / sum(size))
    
    #wait.summary.df %>%
    #  select(cluster, size, alpha)
    
    DensityFaithfulY <- density(wait)$y
    DensityFaithfulX <- density(wait)$x
    
    m <- nlsLM(DensityFaithfulY ~ bimodal(DensityFaithfulX, mu1, sigma1, A1, mu2, sigma2, A2),
               start = list(mu1 = wait.summary.df[["mu"]][1], sigma1  = wait.summary.df[["std"]][1], A1 = wait.summary.df[["alpha"]][1], mu2  = wait.summary.df[["mu"]][2], sigma2 = wait.summary.df[["std"]][2], A2 = wait.summary.df[["alpha"]][2]), weights = 1 / (DensityFaithfulX - 2 * DensityFaithfulX[1] + DensityFaithfulX[2]))
    
    m_res <- coef(m)
    
    min_index <- (which.min(m_res[c(1, 4)]) - 1) * 3 + 1
    max_index <- (which.max(m_res[c(1, 4)]) - 1) * 3 + 1
    
    normaldens_first <- rbind(normaldens_first, data.frame(sample_id = sample_id, feature = grid, density = plot_mix_comps(grid, m_res[1], m_res[2], m_res[3])))
    normaldens_second <- rbind(normaldens_second, data.frame(sample_id = sample_id, feature = grid, density = plot_mix_comps(grid, m_res[4], m_res[5], m_res[6])))
    
    threshold_value <- tryCatch({
      uniroot(function(x) m_res[3] * stats::dnorm(x, mean = m_res[1], sd = m_res[2]) - m_res[6] * stats::dnorm(x, mean = m_res[4], sd = m_res[5]) , c(min(m_res[c(1, 4)]), max(m_res[c(1, 4)])))$root
    }, error=function(cond) {
      tryCatch({
        uniroot(function(x) m_res[3] * stats::dnorm(x, mean = m_res[1], sd = m_res[2]) - m_res[6] * stats::dnorm(x, mean = m_res[4], sd = m_res[5]) , c(m_res[min_index] - 1 * m_res[min_index + 1], m_res[max_index] + 1 * m_res[max_index + 1]))$root
      }, error=function(cond) {
        max(df[[feature]])
      })
    })
    
    annotation_df <- rbind(annotation_df, data.frame(sample_id = sample_id, x=0.1,  y=0.95, lab=paste("Threshold value:", format(round(threshold_value, 2), nsmall = 2))))
    
    #samples_thresholds <- c(samples_thresholds, max(m_res[c(1,4)]))
    samples_thresholds <- c(samples_thresholds, threshold_value)
  }
  
  samples_thresholds_df <- as.data.frame(samples_thresholds)
  samples_thresholds_df$sample_id <- unique(df$sample_id)
  
  sample.labs <- str_replace(str_replace(unique(df$sample_id), "pb037-01-03-Unmixed2-pb037-36-color-", ""), "_Unmixed.fcs", "")
  names(sample.labs) <- unique(df$sample_id)
  
  colnames(normaldens_first)[2] <- feature
  colnames(normaldens_second)[2] <- feature
  
  division_plot <- ggplot(df[, c("sample_id", feature)], aes_string(x = feature)) +
    geom_histogram(aes(y = ..density..), colour = "black", 
                   fill = "white", bins = 50) +
    geom_line(aes(y = density), data = normaldens_first, colour = "red") +
    geom_line(aes(y = density), data = normaldens_second, colour = "blue") +
    facet_wrap(~ sample_id, labeller = labeller(sample_id = sample.labs))+
    geom_vline(data = samples_thresholds_df, aes(xintercept=samples_thresholds)) +
    ylab("Density") +
    xlab("Values") +
    ggtitle(paste(feature, "division", title_info)) +
    geom_text(data = annotation_df, aes(x = x,  y = y, label = lab))  + 
    xlim(xmin, xmax)
  
  return(list("thresholds" = samples_thresholds_df, "plot" = division_plot))
}

fit_one_gaussian_and_find_thresholds <- function(df, feature, title_info = "", xmin = -Inf, xmax = Inf, shift = 3) {
  samples_thresholds <- c()
  
  normaldens <- data.frame()
  annotation_df <- data.frame()
  grid <- seq(min(df[,feature]), max(df[,feature]), length = 100)
  sample_id <- unique(df$sample_id)[1]
  for(sample_id in unique(df$sample_id)) {
    wait <- df[df$sample_id == sample_id, feature]
    
    wait.kmeans <- kmeans(wait, 1)
    wait.kmeans.cluster <- wait.kmeans$cluster
    
    wait.df <- data_frame(x = wait, cluster = wait.kmeans.cluster)
    
    wait.summary.df <- wait.df %>%
      group_by(cluster) %>%
      dplyr::summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    
    #wait.summary.df %>%
    #  select(cluster, mu, variance, std)
    
    wait.summary.df <- wait.summary.df %>%
      mutate(alpha = size / sum(size))
    
    #wait.summary.df %>%
    #  select(cluster, size, alpha)
    
    DensityFaithfulY <- density(wait)$y
    DensityFaithfulX <- density(wait)$x
    
    m <- nlsLM(DensityFaithfulY ~ gauss(DensityFaithfulX, mu, sigma, A), 
               start = list(mu = wait.summary.df[["mu"]][1], sigma  = wait.summary.df[["std"]][1], A = wait.summary.df[["alpha"]][1]), weights = 1 / (DensityFaithfulX - 2 * DensityFaithfulX[1] + DensityFaithfulX[2]))
    
    m_res <- coef(m)
    
    normaldens <- rbind(normaldens, data.frame(sample_id = sample_id, feature = grid, density = plot_mix_comps(grid, m_res[1], m_res[2], m_res[3])))
    
    threshold_value <- m_res[1] + shift * m_res[2]
    
    annotation_df <- rbind(annotation_df, data.frame(sample_id = sample_id, x=0.1,  y=0.95, lab=paste("Threshold value:", format(round(threshold_value, 2), nsmall = 2))))
    
    samples_thresholds <- c(samples_thresholds, threshold_value)
  }
  
  samples_thresholds_df <- as.data.frame(samples_thresholds)
  samples_thresholds_df$sample_id <- unique(df$sample_id)
  
  sample.labs <- str_replace(str_replace(unique(df$sample_id), "pb037-01-03-Unmixed2-pb037-36-color-", ""), "_Unmixed.fcs", "")
  names(sample.labs) <- unique(df$sample_id)
  
  colnames(normaldens)[2] <- feature
  
  division_plot <- ggplot(df[, c("sample_id", feature)], aes_string(x = feature)) +
    geom_histogram(aes(y = ..density..), colour = "black", 
                   fill = "white", bins = 50) +
    geom_line(aes(y = density), data = normaldens, colour = "red") +
    facet_wrap(~ sample_id, labeller = labeller(sample_id = sample.labs))+
    geom_vline(data = samples_thresholds_df, aes(xintercept=samples_thresholds)) +
    ylab("Density") +
    xlab("Values") +
    ggtitle(paste(feature, "division", title_info)) +
    geom_text(data = annotation_df, aes(x = x,  y = y, label = lab)) +
    xlim(xmin, xmax)
  return(list("thresholds" = samples_thresholds_df, "plot" = division_plot))
}