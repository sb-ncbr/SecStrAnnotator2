library(dplyr)  # Tested with dplyr 0.8.5
library(tidyr)  # Tested with tidyr 1.0.2
library(readr)  # Tested with readr 1.3.1
library(ggplot2)  # Tested with ggplot2 3.3.0

#############################################################################
# CytochromeP450-specific constants and functions

HELIX_ORDER = c("A'", "A", "B", "B'", "B''", "C", "D", "E", "F", "F'", "G'", "G", "H", "I", "J", "J'", "K", "K'", "K''", "L", "L'")
SHEET_ORDER = c("1-0", "1-1", "1-2", "1-3", "1-4", "1-5", "2-1", "2-2", "3-1", "3-2", "3-3", "4-1", "4-2", "5-1", "5-2", "6-1", "6-2")
OUR_SSES = c(HELIX_ORDER, SHEET_ORDER)

# Decide if the SSE label corresponds to major, minor helix, or strand
sse_category = function(name) {
  name = as.character(name)
  n = nchar(name)
  if (is_digit(substr(name, 1, 1))) 'strand'
  else if (substr(name, n, n)=='\'') 'minor h.' 
  else 'major h.'
}

#############################################################################
# General functions

NONRETARDED_PAIRED_PALETTE = c("#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")  # Colors from palette Paired, ordered consistently with Set1
NONRETARDED_PAIRED_PALETTE_SEPARATED = c("#FB9A99", "#A6CEE3", "#B2DF8A", "#E31A1C", "#1F78B4", "#33A02C")  # Colors from palette Paired, ordered consistently with Set1, first light then dark
NONRETARDED_PAIRED_PALETTE_DARK = c("#FB6A69", "#C30A0C", "#769EE3", "#0F68A4", "#82CF5A", "#23901C")  # Colors from palette Paired, ordered consistently with Set1, first light then dark, darkened
NONRETARDED_PAIRED_PALETTE_SEPARATED_DARK = c("#FB6A69", "#769EE3", "#82CF5A", "#C30A0C", "#0F68A4", "#23901C")  # Colors from palette Paired, ordered consistently with Set1, first light then dark, darkened
HELIX_TYPE_PALETTE = c("#B080C8", "#30A028", "#F8A870", "#CCCCCC") # bright violet (310-helix), dark green (alpha-helix), orange (pi-helix), grey (other)


# Interpret path as starting in DATADIR
full_path = function(path) { 
  paste(DATADIR, path, sep='/') 
}


# Get list of domain names and their UniProtIDs. Optionally summarize domains for each UniProtID to one row.
get_domains = function(sse_annotation_df, summarize_by_UniProt=FALSE) {
  domains = sse_annotation_df %>% select('UniProt', 'Domain') %>% distinct() %>% arrange(UniProt, Domain)
  if (summarize_by_UniProt) {
    domains = domains %>% group_by(UniProt) %>% summarize(Count = length(Domain), Domains = paste(Domain, collapse=','))
  }
  domains
}

# Select only one (first) domain for each PDB ID
keep_one_domain_per_pdb = function(sse_annotation_df) {
  selected_domains = sse_annotation_df %>% group_by(PDB) %>% summarize(Domain=first(Domain)) %>% select(Domain)
  filtered_df = right_join(sse_annotation_df, selected_domains)
  filtered_df
}

# Compare two SSE annotation dataframes
compare_annotations = function(dfX, dfY){
  oldies = anti_join(dfX %>% select(Domain, label, start, end, type, length), dfY %>% select(Domain, label, start, end, type, length))
  newbies = anti_join(dfY %>% select(Domain, label, start, end, type, length), dfX %>% select(Domain, label, start, end, type, length))
  diff = full_join(oldies, newbies, by=c('Domain', 'label'))
  diff
}

# Calculate absolute and relative occurrence for each SSE label + confidence intervals
table_sse_occurrence = function(df, alpha = 0.05) {
  if (!is.null(df$Domain)) {
    n_pdb = n_distinct(df$Domain)
  } else {
    warning('Column Domain not found, using column PDB to count domains')
    n_pdb = n_distinct(df$PDB)
  }
  data = df %>% filter(length > 0) %>% group_by(label) %>% count() %>% 
    mutate(freq = n/n_pdb, type = sse_category(label), freq_min = agresti_coull_confidence_interval(n, n_pdb, alpha)[1], freq_max = agresti_coull_confidence_interval(n, n_pdb, alpha)[2])
  data
}

# Calculate confidence intervals for binomial distribution
agresti_coull_confidence_interval = function(successes, total, alpha = 0.05) {
  z = qnorm(1-alpha/2)
  n_ = total + z^2
  p_ = 1/n_ * (successes + z^2/2)
  d = z * sqrt(p_ * (1-p_) / n_)
  c(p_-d, p_+d)
}

# Decide if character is a letter
is_letter = function(x) {
  grepl("[[:alpha:]]", x)
}

# Decide if character is a digit
is_digit = function(x) {
  grepl("[[:digit:]]", x)
}

# Change first letter of each string to uppercase, rest to lowercase
capitalize = function(x) { 
  paste(toupper(substring(x, 1, 1)), tolower(substring(x, 2)), sep='')
}

# Barchart with relative occurrence for each SSE label
plot_sse_occurrence = function(df, show_confidence = TRUE, alpha = 0.05, title = NULL, stagger_labels = FALSE, turn_labels = FALSE) {
  g = ggplot(table_sse_occurrence(df, alpha)) +
    aes(x = label, y = 100*freq, fill = type) +
    geom_col() + 
    geom_vline(xintercept = seq(1.5, length(unique(df$label))-0.5, 1), lwd = 0.4, colour = 'gray90') +
    scale_fill_brewer(palette = 'Set1') + 
    coord_cartesian(ylim = c(0, 105)) +
    scale_y_continuous(expand = expand_scale(add = 0), breaks = seq(0, 100, by = 20)) +
    labs(x = 'SSE', y = 'Occurrence [%]', fill = 'Type', title = title) +
    order_x_labels(stagger = stagger_labels) +
    theme_bw() +
    customize_theme(stagger_labels = stagger_labels, turn_labels = turn_labels)
  if (show_confidence) { 
    g = g + geom_errorbar(aes(ymin = 100*freq_min, ymax = 100*freq_max), width = 0.4, position = position_dodge(.9)) 
  }
  g
}

# Barchart with relative occurrence for each SSE label, comparing multiple datasets
plot_sse_occurrence_multi = function(..., show_confidence = TRUE, alpha = 0.05, title = NULL, stagger_labels = FALSE, turn_labels = FALSE) {
  tables = lapply(pairlist(...), table_sse_occurrence, alpha)
  table = bind_rows(tables, .id = 'set')
  g = ggplot(table) +
    aes(x = label, y = 100*freq, fill = paste(set, ' (', type, ')', sep = '')) +
    geom_col(position = 'dodge') + 
    geom_vline(xintercept = seq(1.5, length(unique(table$label))-0.5, 1), lwd = 0.4, colour = 'gray90') +
    scale_fill_manual(values = NONRETARDED_PAIRED_PALETTE_SEPARATED) + 
    coord_cartesian(ylim = c(0,105)) +
    scale_y_continuous(expand = expand_scale(add = 0), breaks = seq(0, 100, by = 20)) +
    labs(x = 'SSE', y = 'Occurrence [%]', fill = 'Type', title = title) +
    order_x_labels(stagger = stagger_labels) +
    theme_bw() +
    customize_theme(stagger_labels = stagger_labels, turn_labels = turn_labels)
  if (show_confidence) { 
    g = g + geom_errorbar(aes(ymin = 100*freq_min, ymax = 100*freq_max), width = 0.4, position = position_dodge(.9)) 
  }
  g
}

# Barchart with relative occurrence for each beta-bulge label
plot_bulge_occurrence = function(df, bulge_df, show_confidence = TRUE, alpha = 0.05, title = NULL, include_parent_strands = TRUE, turn_labels = FALSE, stagger_labels = FALSE) {
  g = ggplot(table_bulge_occurrence(df, bulge_df, alpha = alpha, include_parent_strands = include_parent_strands)) +
    aes(x = label, y = 100*freq, fill = type) +
    geom_col() + 
    geom_vline(xintercept = seq(1.5, length(unique(df$label))-0.5, 1), lwd = 0.4, colour = 'gray90') +
    scale_fill_brewer(palette = 'Paired') + 
    labs(x = 'SSE', y = 'Occurrence [%]', fill = 'Bulge type', title = title) +
    theme_bw() +
    customize_theme(stagger_labels = stagger_labels, turn_labels = turn_labels)
  if (show_confidence) { 
    g = g + geom_errorbar(aes(ymin = 100*freq_min, ymax = 100*freq_max), width = 0.4, position = position_dodge(.9)) }
  g
}

# Barchart with relative occurrence of each helix type (3_10/alpha/pi) for each SSE label
plot_contained_helix_types = function(df, ..., turn_labels = FALSE, stagger_labels = FALSE){
  changed_labels = c(contains_G = expression(paste('contains ', '3'[10])), contains_Η = expression(paste('contains ', alpha)), contains_Ι = expression(paste('contains ', pi)) )
  g = plot_criteria(df, contains_G = (df$longest_G >= 3), contains_Η = (df$longest_H >= 4), contains_Ι = (df$longest_I >= 5), ..., ignore_zero = TRUE)  +
    scale_fill_manual(values = HELIX_TYPE_PALETTE, labels = changed_labels) + 
    order_x_labels(include_strands = FALSE) + 
    theme_bw() +
    customize_theme(stagger_labels = stagger_labels, turn_labels = turn_labels)+
    theme(legend.text.align = 0)
  g
}

# Barchart with relative fraction of cases fulfulling specified criteria, for each SSE label
plot_criteria = function(df, ..., ignore_zero = FALSE, title = NULL, y_label = 'Fulfilling criterion'){
  table = table_criteria(df, ..., ignore_zero = ignore_zero) %>% group_by(label, set) %>% mutate(type = sse_category(label))
  g = ggplot(table) + 
    aes(x = label, y = 100*fulfilling, fill = set) +
    geom_col(position = 'dodge') + 
    geom_vline(xintercept = seq(1.5, length(unique(df$label))-0.5, 1), lwd = 0.4, colour = 'gray90') + 
    scale_fill_brewer(palette = 'Set1') + 
    coord_cartesian(ylim = c(0,105)) +
    scale_y_continuous(expand = expand_scale(add = 0), breaks = seq(0, 100, by = 20)) +
    labs(x = 'SSE', y = paste(y_label, ' [%]'), fill = 'Type', title = title) +
    order_x_labels()
  g
}

# Table with relative fraction of cases fulfulling specified criteria, for each SSE label
table_criteria = function(df, ..., by_label = TRUE, ignore_zero = FALSE) {
  tables = lapply(pairlist(...), function(crit) {table_criterion(df, crit, by_label = by_label, ignore_zero = ignore_zero)})
  t = bind_rows(tables, .id = 'set')
  t
}

# Table with relative fraction of cases fulfulling specified criterion, for each SSE label
table_criterion = function(df, criterion, by_label = TRUE, ignore_zero = FALSE) {
  df$ok = criterion
  if (ignore_zero){
    df = df %>% filter(length>0)
  }
  if (by_label){
    df = df %>% group_by(label)
  }
  df = df %>% summarise(fulfilling=mean(ok))
  df
}

# Table with relative occurrence for each beta-bulge label
table_bulge_occurrence = function(df, bulge_df, alpha = 0.05, include_parent_strands = TRUE) {
  bulge_df = mutate(bulge_df, label = paste(label, type, sep = '.'))
  if (include_parent_strands){
    bulge_df = bind_rows(mutate(df, type = 'strand'), bulge_df)
  }
  n_pdb = n_distinct(df$PDB)
  data = bulge_df %>% filter(length > 0) %>% group_by(label, type) %>% distinct(PDB) %>% count() %>% 
    mutate(freq = n/n_pdb, category = sse_category(label), freq_min = agresti_coull_confidence_interval(n, n_pdb, alpha)[1], freq_max = agresti_coull_confidence_interval(n, n_pdb, alpha)[2]) %>%
    filter(category == 'strand') %>% select(-'category')
  data
}

# Order SSE labels according to OUR_SSES, optionaly convert apostrophest to prime symbols or stagger labels (show every other label lower)
order_x_labels = function(apostrophe_to_prime = TRUE, stagger = FALSE, include_helices = TRUE, include_strands = TRUE) { 
  f = if (apostrophe_to_prime){
    function(x) { sub("'", "′", sub("''", "″", sub("'''", "‴", x))) }
  } else {
    function(x) { x }
  }
  label_func = if (stagger) {
    function(x) {paste(c('', '\n'), f(x))}
  } else {
    f
  }
  limits = if (include_helices & include_strands) {OUR_SSES} else if (include_helices) {HELIX_ORDER} else if (include_strands) {SHEET_ORDER}
  scale_x_discrete(limits = limits, labels = label_func)
}

# Customize theme_bw in ggplot
customize_theme = function(stagger_labels = FALSE, turn_labels = FALSE){
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = 'lightgray'),
    legend.box.margin = margin(l = -10, r = -8),
    axis.text.x = element_text(angle = if (turn_labels) {45} else {0}, margin = margin (t = if (turn_labels) {5} else {2}))
  )
}

# Boxplot with SSE length (or other specified column) distribution, for each SSE label
boxplot_sse = function(df, y_column = 'length', ignore_zero = FALSE, width = 0.9, group_col = 'type', title = NULL) {
  df = df %>% mutate(type = mapply(sse_category, label))
  if (ignore_zero) df = df %>% filter(df$length>0)
  ggplot(df) +
    aes(x = label, y = df[[y_column]], color = df[[group_col]]) +
    geom_boxplot(width = width) + 
    geom_vline(xintercept = seq(1.5, length(unique(df$label))-0.5, 1), lwd = 0.4, colour = 'gray90') +
    scale_color_brewer(palette = 'Set1') + 
    coord_cartesian(ylim = c(0, max(df[[y_column]]) + 1)) +
    scale_y_continuous(expand = expand_scale(add=0), minor_breaks = seq(0 , 100, 1), breaks = seq(0, 100, 5)) +
    labs(x = 'SSE', y = 'Length [#residues]', color = capitalize(group_col), title = title) +
    order_x_labels() + 
    theme_bw() + 
    customize_theme()
}

# Boxplot with SSE length (or other specified column) distribution, for each SSE label
violinplot_sse = function(df, y_column = 'length', ignore_zero = FALSE, scale = 'width', width = 0.9, group_col = 'type', title = NULL, stagger_labels = FALSE, turn_labels = FALSE) {
  #width = max. width wrt. bar width
  #scale = 'width' (violins have equal width), 'area' (violins have equal area)
  #point is mean, line is median
  df = df %>% mutate(type = mapply(sse_category, label))
  if (ignore_zero) df = df %>% filter(df$length > 0)
  ggplot(df) +
    aes(x = label, y = df[[y_column]], color = df[[group_col]]) +
    geom_violin(scale=scale, width = width, draw_quantiles = c(0.5), kernel = 'triangular', bw = 1/sqrt(6)) + #kernel='gaussian',bw=0.45 or kernel='triangular',bw=1/sqrt(6)
    geom_vline(xintercept = seq(1.5, length(unique(df$label))-0.5, 1), lwd = 0.4, colour = 'gray90') +
    stat_summary(fun.y = 'mean', geom = 'point') + 
    scale_color_brewer(palette = 'Set1') + 
    coord_cartesian(ylim = c(0, max(df[[y_column]]) + 1)) +
    scale_y_continuous(expand = expand_scale(add = 0), minor_breaks = seq(0 , 100, 1), breaks = seq(0, 100, 5)) +
    labs(x = 'SSE', y = 'Length [#residues]', color = capitalize(group_col), title = title) +
    order_x_labels(stagger = stagger_labels) +
    theme_bw()  +
    customize_theme(stagger_labels = stagger_labels, turn_labels = turn_labels)
}

# Boxplot with SSE length (or other specified column) distribution for each SSE label, comparing multiple datasets
violinplot_sse_multi = function(..., y_column = 'length', set_names = NULL, ignore_zero = FALSE, scale = 'width', width = 0.9, title = NULL, stagger_labels = FALSE, turn_labels = FALSE) {
  #width = max. width wrt. bar width
  #scale = 'width' (violins have equal width), 'area' (violins have equal area)
  #point is mean, line is median
  df = bind_rows(..., .id = 'set')
  if (ignore_zero) df = df %>% filter(length > 0)
  df = df %>% mutate(type = mapply(sse_category, label))
  ggplot(df) +
    aes(x = label, y = length, color = paste(set, ' (', type, ')', sep = '')) + 
    geom_violin(scale = scale, width = width, draw_quantiles = c(0.5), kernel = 'triangular', bw = 1/sqrt(6)) + #kernel='gaussian',bw=0.45 or kernel='triangular',bw=1/sqrt(6)
    geom_vline(xintercept = seq(1.5, length(unique(df$label))-0.5, 1), lwd = 0.4, colour = 'gray90') +
    stat_summary(fun.y = 'mean', geom = 'point') + 
    scale_color_manual(values = NONRETARDED_PAIRED_PALETTE_SEPARATED_DARK) + 
    coord_cartesian(ylim = c(0, max(df$length) + 1)) +
    scale_y_continuous(expand = expand_scale(add = 0), breaks = seq(0, 100, 5), minor_breaks = seq(0 , 100, 1)) +
    labs(x = 'SSE', y = 'Length [#residues]', color = 'Type', title = title) + 
    order_x_labels(stagger = stagger_labels) +
    theme_bw() + 
    customize_theme(stagger_labels = stagger_labels, turn_labels = turn_labels)
}

# Compare two binomial distributions by the test of equal proportions (prop.test), for each label
two_sample_occurrence_prop_test = function(df1, df2, p_limit = 0.05, print_all = FALSE, signif_digits = 2) {
  n_pdb1 = n_distinct(df1$PDB)
  data1 = df1 %>% filter(length > 0) %>% group_by(label) %>% count() %>% 
    mutate(present1 = n, absent1 = n_pdb1-n, freq1 = n/n_pdb1) %>% select(-n)
  n_pdb2 = n_distinct(df2$PDB)
  data2 = df2 %>% filter(length > 0) %>% group_by(label) %>% count() %>% 
    mutate(present2 = n, absent2 = n_pdb2-n, freq2 = n/n_pdb2) %>% select(-n)
  data = full_join(data1, data2, by = 'label') %>% mutate(p.value = 0, comparison = '', freq_diff = freq2-freq1)
  n_labels = length(data$label)
  for (the_label in unique(data$label)) {
    row = filter(data, label == the_label)
    test_result = prop.test(c(row$present1, row$present2), c(n_pdb1, n_pdb2))
    p = test_result$p.value
    freq_diff = row$freq_diff
    data[which(data$label == the_label), 'p.value'] = if (is.na(p)) '' else signif(p, digits=signif_digits) 
    data[which(data$label == the_label), 'compare'] = if (p>p_limit | freq_diff==0) '       ' else if (freq_diff<0) '  >    ' else if (freq_diff>0) '    <  '
    data[which(data$label == the_label), 'freq_diff'] = signif(freq_diff, digits=signif_digits)
  }
  data$freq1 = round(data$freq1, digits=signif_digits+1)
  data$freq2 = round(data$freq2, digits=signif_digits+1)
  
  if (!print_all) {
    data = data %>% filter(p.value<=p_limit)
  }
  data %>% select(label, freq1, compare, freq2, freq_diff, p.value) %>% data.frame()
}

# Compare two binomial distributions by the Fisher's test (fisher.test), for each label
two_sample_occurrence_fisher_test = function(df1, df2, p_limit = 0.05, print_all = FALSE, signif_digits = 2) {
  n_pdb1 = n_distinct(df1$PDB)
  data1 = df1 %>% filter(length>0) %>% group_by(label) %>% count() %>% 
    mutate(present1 = n, absent1 = n_pdb1-n, freq1 = n/n_pdb1) %>% select(-n)
  n_pdb2 = n_distinct(df2$PDB)
  data2 = df2 %>% filter(length>0) %>% group_by(label) %>% count() %>% 
    mutate(present2 = n, absent2 = n_pdb2-n, freq2 = n/n_pdb2) %>% select(-n)
  data = full_join(data1, data2, by='label') %>% mutate(p.value = 0, comparison = '', freq_diff = freq2-freq1)
  n_labels = length(data$label)
  for (the_label in data$label) {
    row = filter(data, label == the_label)
    # View(row)
    the_matrix = matrix(c(row$present1, row$absent1, row$present2, row$absent2), ncol=2)
    test_result = fisher.test(the_matrix)
    p = test_result$p.value
    freq_diff = row$freq_diff
    data[which(data$label == the_label), 'p.value'] = signif(p, digits=signif_digits) 
    data[which(data$label == the_label), 'compare'] = if (p>p_limit | freq_diff==0) '       ' else if (freq_diff<0) '  >    ' else if (freq_diff>0) '    <  '
    data[which(data$label == the_label), 'freq_diff'] = signif(freq_diff, digits=signif_digits) 
  }
  if (!print_all) {
    data = data %>% filter(p.value<=p_limit)
  }
  data %>% select(label, freq1, compare, freq2, freq_diff, p.value) %>% data.frame()
}

# Compare two distributions by the specified test, for each label
two_sample_test_by_labels_with_comparison = function(test, df1, df2, stat_col, label_col='label', p_limit=1, print_all=FALSE, reversed=FALSE, signif_digits=2) {
  result = unique(df1[label_col]) %>% mutate(mean1 = NaN, compare = '', mean2 = NaN, mean_diff = 0, median_diff = 0, p.value = NaN, p.g = NaN, p.l = NaN)
  result = result %>% arrange(result[[label_col]])
  n = length(result[[label_col]])
  for (label in result[[label_col]]) {
    sample1 = df1[which(df1[label_col]==label), stat_col][[stat_col]]
    sample2 = df2[which(df2[label_col]==label), stat_col][[stat_col]]
    p = test(sample1, sample2)$p.value
    pG = test(sample1, sample2, alternative='greater')$p.value
    pL = test(sample1, sample2, alternative='less')$p.value
    result[which(result[label_col]==label), 'p.value'] = signif(p, digits=signif_digits) 
    result[which(result[label_col]==label), 'p.g'] = signif(pG, digits=signif_digits) 
    result[which(result[label_col]==label), 'p.l'] = signif(pL, digits=signif_digits) 
    GT = '  >    '
    LT = '    <  '
    EQ = '       '
    comparison = if (p<p_limit & pG<p_limit & pL>=p_limit) GT else if (p<p_limit & pL<p_limit & pG>=p_limit) LT else EQ
    if (reversed) {
      comparison = if (comparison==GT) LT else if (comparison==LT) GT else EQ
    }
    result[which(result[label_col]==label), 'compare'] = comparison
    m1 = mean(sample1)
    m2 = mean(sample2)
    med1 = median(sample1)
    med2 = median(sample2)
    result[which(result[label_col]==label), 'mean1'] = signif(m1, digits=signif_digits)
    result[which(result[label_col]==label), 'mean2'] = signif(m2, digits=signif_digits)
    result[which(result[label_col]==label), 'mean_diff'] = signif(m2 - m1, digits=signif_digits)
    result[which(result[label_col]==label), 'median_diff'] = signif(med2 - med1, digits=signif_digits)
  }
  if (!print_all) {
    result = result %>% filter(p.value<=p_limit)
  }
  result %>% data.frame()
}

# Save as a PNG file
print_png = function(filename, width = 2000, ratio = 2/1, height = NULL, res = 250) { 
  if (is.null(height)) { height = width/ratio }
  dev.print(png, filename, width = width, height = height, res = res)
}

# Save as a TIFF file
print_tif = function(filename, width = 2000, ratio = 2/1, height = NULL, res = 250) { 
  if (is.null(height)) { height = width/ratio }
  dev.print(tiff, filename, width = width, height = height, res = res)
}
