# library(Hmisc)
library(dplyr)
library(tidyr)
library(ggplot2)
#library(readr)
#library(RColorBrewer)
#library(EMT)
#library(diptest)


fullname = function(filename) { 
  paste(DATADIR, filename, sep='/') 
}

read_sses_json = function(annot_file, uniprot_file, fields, vector_fields=NULL, add_missing_sses=TRUE) {
  syscall = paste('python3 to_csv.py', annot_file, uniprot_file, paste(cbind(fields,vector_fields),collapse=','), (if (add_missing_sses) 'True' else 'False'), sep=' ')
  csv_string = system(syscall, intern=TRUE)
  df = read.csv(text = csv_string, sep='\t')
  if (!is.null(vector_fields)) {
    df = split_vector_columns(df, vector_fields)
  }
  df
}

read_sses_json_new = function(annot_file, fields, vector_fields=NULL, add_missing_sses=TRUE, one_domain_per_pdb=FALSE, only_domain_lists=FALSE) {
  syscall = paste('python3 to_csv_new.py', annot_file, '--fields', paste(cbind(fields,vector_fields),collapse=','), (if (add_missing_sses) '--add_missing_sses'), (if (one_domain_per_pdb) '--one_domain_per_pdb'), (if (only_domain_lists) '--only_domain_lists'), sep=' ')
  csv_string = system(syscall, intern=TRUE)
  df = read.csv(text = csv_string, sep='\t')
  if (!is.null(vector_fields)) {
    df = split_vector_columns(df, vector_fields)
  }
  df
}

read_bulges_json = function(annot_file, uniprot_file, fields, vector_fields=NULL) {
  syscall = paste("python3 to_csv_bulges.py", annot_file, uniprot_file, paste(cbind(fields,vector_fields),collapse=','), sep=' ')
  csv_string = system(syscall, intern=TRUE)
  df = read.csv(text = csv_string, sep='\t')
  if (!is.null(vector_fields)) {
    df = split_vector_columns(df, vector_fields)
  }
  df
}

read_bulges_json_new = function(annot_file, fields, one_domain_per_pdb=FALSE) {
  syscall = paste('python3 to_csv_bulges_new.py', annot_file, '--fields', paste(fields, collapse=','), (if (one_domain_per_pdb) '--one_domain_per_pdb'), sep=' ')
  #syscall = paste("python3 to_csv_bulges_new.py", annot_file, paste(fields, collapse=','), sep=' ')
  csv_string = system(syscall, intern=TRUE)
  df = read.csv(text = csv_string, sep='\t')
  df
}

split_vector_columns = function(df, vector_cols) {
  for (vecf in vector_cols) {
    #df %>% mutate(vecf = 'a' ) %>% separate(vecf,c(paste(vecf,'_x'), paste(vecf,'_y'), paste(vecf,'_z')),sep=',') %>% View()
    new_col_names = c(paste(vecf,'_x',sep=''), paste(vecf,'_y',sep=''), paste(vecf,'_z',sep=''))
    df = df %>% mutate(vecf = gsub(' ', '', gsub('\\[', '', gsub(']', '', df[,vecf]))) ) %>% 
      separate(vecf, new_col_names,sep=',', remove=TRUE) %>%
      mutate_at(new_col_names, function(x){as.numeric(x)})
  }
  df[, !(colnames(df) %in% vector_cols)]
}

sse_length = function(start, end) {
  if (is.na(start))
    result = 0
  else
    result = end - start + 1
  result
}

add_lengths = function(df) {
  df$length = mapply(sse_length, df$start, df$end)
  df
}
is_letter <- function(x) grepl("[[:alpha:]]", x)
is_digit <- function(x) grepl("[[:digit:]]", x)

sse_category = function(name) {
  name = as.character(name)
  n=nchar(name)
  if (is_digit(substr(name,1,1)))
    result = 'strand'
  else if (substr(name,n,n)=='\'')
    result = 'minor h.' 
  else
    result = 'major h.'
  result
}

HELIX_ORDER = c("A'", "A", "B", "B'", "B''", "C", "D", "E", "F", "F'", "G'", "G", "H", "I", "J", "J'", "K", "K'", "K''", "L", "L'")
SHEET_ORDER = c("1-0", "1-1", "1-2", "1-3", "1-4", "1-5", "2-1", "2-2", "3-1", "3-2", "3-3", "4-1", "4-2", "5-1", "5-2", "6-1", "6-2")
sse_order = c("A'", "A", "1a", "1b", "B", "1e", "B''", "B'", "B'''", "C", "D", "3a", "E", "F", "F'", "G'", "G", "H", "I", "J", "J'", "K", "4c", "1d", "2a", "2b", "1c", "K'", "K'''", "K''", "L", "3c", "L'", "4a", "4b", "3b", "h", "e")
OUR_SSES = c(HELIX_ORDER, SHEET_ORDER)

nonretarded_paired_palette = c("#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")  # Colors from palette Paired, ordered consistently with Set1
nonretarded_paired_palette_separated      = c("#FB9A99", "#A6CEE3", "#B2DF8A", "#E31A1C", "#1F78B4", "#33A02C")  # Colors from palette Paired, ordered consistently with Set1, first light then dark
nonretarded_paired_palette_separated_dark = c("#FB6A69", "#769EE3", "#82CF5A", "#C30A0C", "#0F68A4", "#23901C")  # Colors from palette Paired, ordered consistently with Set1, first light then dark, darkened

only_helices = function(df){
  df %>% filter(label %in% HELIX_ORDER)
}
only_sheets = function(df){
  df %>% filter(label %in% SHEET_ORDER)
}

order_x_labels_old = function(apostrophe_to_prime=TRUE) { 
  if (apostrophe_to_prime){
    scale_x_discrete(limits=c(HELIX_ORDER, SHEET_ORDER), 
                     labels=function(x){ paste(c('', '\n'), sub("'", "′", sub("''", "″", sub("'''", "‴", x)))) }
    )
  } else {
    scale_x_discrete(limits=c(HELIX_ORDER, SHEET_ORDER))
  }
}

order_x_labels = function(apostrophe_to_prime=TRUE, stagger=FALSE) { 
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
  scale_x_discrete(limits = c(HELIX_ORDER, SHEET_ORDER), labels = label_func)
}

order_x_labels_helices = function() { scale_x_discrete(limits=HELIX_ORDER) }
order_fills = function() { scale_fill_discrete(limits=c('mainly H', 'mainly G', 'other')) }

table_sse_lengths = function(df, ignore_zero=FALSE) {
  if (ignore_zero) df = df %>% filter(length>0)
  data = df %>% group_by(label) %>% summarise(avg_len=mean(length),sd_len=sd(length),type=sse_category(label[1]))
}

plot_sse_lengths = function(df, ignore_zero=FALSE) {
  ggplot(table_sse_lengths(df,ignore_zero)) +
    aes(x=label,y=avg_len,ymin=avg_len-sd_len,ymax=avg_len+sd_len,fill=type) +
    geom_col() + 
    scale_fill_brewer(palette='Set1') + 
    #geom_point() + 
    geom_errorbar(width=0.4) +
    labs(x='SSE', y='Mean length [#residues]', fill='Type') + 
    order_x_labels()
}

jitterplot_sse_lengths = function(df, ignore_zero=FALSE) {
  df = df %>% mutate(type=mapply(sse_category,label))
  if (ignore_zero) df = df %>% filter(length>0)
  ggplot(df) +
    aes(x=label,y=length,color=type) +
    geom_jitter(size=0.4) + 
    scale_color_brewer(palette='Set1') + 
    labs(x='SSE', y='Length [#residues]', color='Type') + 
    order_x_labels()
}

jitterplot_xy = function(df, x, y, ignore_zero=FALSE) {
  x = enquo(x)
  y = enquo(y)
  df = df %>% mutate(type=mapply(sse_category,label))
  if (ignore_zero) df = df %>% filter(length>0)
  ggplot(df) +
    #aes(x=df[[x_column]], y=df[[y_column]], color=type) +
    aes(x=!!x, y=!!y, color=type) +
    geom_jitter(size=0.5, height=0.25, width=0.25) + 
    scale_color_brewer(palette='Set1') + 
    labs(x=quo_name(x), y=quo_name(y), color='Type')
}

customize_theme = function(stagger_labels = FALSE, turn_labels = FALSE){
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = 'lightgray'),
    legend.box.margin = margin(l = -10, r = -8),
    axis.text.x = element_text(angle = if (turn_labels) {45} else {0}, margin = margin (t = if (turn_labels) {5} else {2}))
  )
}

violinplot_sse_lengths = function(df, ignore_zero=FALSE, scale='width', width=0.9, group_col='type', title=NULL, stagger_labels=FALSE, turn_labels=FALSE) {
  violinplot_sse(df, 'length', ignore_zero=ignore_zero, scale=scale, width=width, group_col=group_col, title=title, stagger_labels=stagger_labels, turn_labels=turn_labels) 
}

violinplot_sse = function(df, y_column, ignore_zero = FALSE, scale = 'width', width = 0.9, group_col = 'type', title = NULL, stagger_labels = FALSE, turn_labels = FALSE) {
  #width = max. width wrt. bar width
  #scale = 'width' (violins have equal width), 'area' (violins have equal area)
  #point is mean, line is median
  df = df %>% mutate(type = mapply(sse_category, label))
  if (ignore_zero) df = df %>% filter(df$length > 0)
  ggplot(df) +
    aes(x = label, y = df[[y_column]], color = df[[group_col]])
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

boxplot_sse = function(df, y_column, ignore_zero = FALSE, scale = 'width', width = 0.9, group_col = 'type', title = NULL) {
  df = df %>% mutate(type = mapply(sse_category, label))
  if (ignore_zero) df = df %>% filter(df$length>0)
  ggplot(df) +
    aes(x = label, y = df[[y_column]], color = df[[group_col]]) +
    geom_boxplot(scale = scale, width = width) + 
    scale_color_brewer(palette = 'Set1') + 
    scale_y_continuous(minor_breaks = seq(0 , 100, 1), breaks = seq(0, 100, 5)) +
    labs(x = 'SSE', y = 'Length [#residues]', color = capitalize(group_col), title = title) +
    order_x_labels()
}

capitalize = function(x) { 
  paste(toupper(substring(x, 1, 1)), tolower(substring(x, 2)), sep='')
}

violinplot_sse_lengths_multiset = function(..., ignore_zero=FALSE, scale='width', width=0.9, title='') {
  #width = max. width wrt. bar width
  #scale = 'width' (violins have equal width), 'area' (violins have equal area)
  #point is mean, line is median
  
  #df = df %>% mutate(type=mapply(sse_category,label))
  df = bind_rows(..., .id='Dataset')
  if (ignore_zero) df = df %>% filter(length>0)
  ggplot(df) +
    aes(x=label, y=length, color=Dataset) + 
    geom_violin(scale=scale, width=width, draw_quantiles = c(0.5), kernel='triangular', bw=1/sqrt(6)) + #kernel='gaussian',bw=0.45 or kernel='triangular',bw=1/sqrt(6)
    stat_summary(fun.y='mean', geom='point') + 
    scale_color_brewer(palette='Set1') + 
    scale_y_continuous(minor_breaks = seq(0 , 100, 1), breaks = seq(0, 100, 5)) +
    labs(x='SSE', y='Length [#residues]', title=title) + 
    order_x_labels()
}

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
    scale_color_manual(values = nonretarded_paired_palette_separated_dark) + 
    coord_cartesian(ylim = c(0, max(df$length) + 1)) +
    scale_y_continuous(expand = expand_scale(add = 0), breaks = seq(0, 100, 5), minor_breaks = seq(0 , 100, 1)) +
    labs(x = 'SSE', y = 'Length [#residues]', color = 'Type', title = title) + 
    order_x_labels(stagger = stagger_labels) +
    theme_bw() + 
    customize_theme(stagger_labels = stagger_labels, turn_labels = turn_labels)
}

find_outlier_thresholds = function(df, y_column, group_column) {
  y_column = enquo(y_column)  # This is dark magic! https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
  group_column = enquo(group_column)  # This is dark magic! https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
  df %>% group_by(!!group_column) %>% summarise(ymin = boxplot.stats(!!y_column)$stats[1], ymax = boxplot.stats(!!y_column)$stats[5])
}

remove_outliers = function(df, y_column, ignore_zero=FALSE, extra_tolerance=0){
  y_column = enquo(y_column)  # This is dark magic! https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
  if (ignore_zero) df = df %>% filter(!!y_column > 0)
  thresholds = df %>% 
    group_by(label) %>% 
    summarise(ymin = boxplot.stats(!!y_column)$stats[1] - extra_tolerance, ymax = boxplot.stats(!!y_column)$stats[5] + extra_tolerance)
  df = left_join(df, thresholds, by='label')
  df = df %>% filter(!!y_column >= ymin & !!y_column <= ymax) 
  df = df %>% select(-c('ymin','ymax'))
  df
}

table_sse_occurrence = function(df, alpha=0.05) {
  if (!is.null(df$Domain)) {
    n_pdb = n_distinct(df$Domain)
  } else {
    warning('Column Domain not found, using column PDB to count domains')
    n_pdb = n_distinct(df$PDB)
  }
  data = df %>% filter(length>0) %>% group_by(label) %>% count() %>% 
    mutate(freq=n/n_pdb, type=sse_category(label), freq_min=agresti_coull_confidence_interval(n,n_pdb,alpha)[1], freq_max=agresti_coull_confidence_interval(n,n_pdb,alpha)[2])
  data
  }

table_bulge_occurrence = function(df, bulge_df, alpha=0.05, include_parent_strands=TRUE) {
  bulge_df = mutate(bulge_df, label=paste(label, bulge, sep='.'))
  if (include_parent_strands){
    bulge_df = bind_rows (mutate(df, bulge='Strand'), bulge_df)
  }
  n_pdb = n_distinct(df$PDB)
  data = bulge_df %>% filter(length>0) %>% group_by(label, bulge) %>% distinct(PDB) %>% count() %>% 
    mutate(freq=n/n_pdb, type=sse_category(label), freq_min=agresti_coull_confidence_interval(n,n_pdb,alpha)[1], freq_max=agresti_coull_confidence_interval(n,n_pdb,alpha)[2]) %>%
    filter(type=='strand')
  data
}

table_criterion = function(df, criterion, by_label=TRUE, ignore_zero=FALSE) {
  #df = df %>% mutate(ok=criterion) 
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

table_criteria = function(df, ..., by_label=TRUE, ignore_zero=FALSE) {
  tables = lapply(pairlist(...), function(crit){table_criterion(df, crit, by_label=by_label, ignore_zero=ignore_zero)})
  t = bind_rows(tables, .id='set')
  #t = spread(t, set, fulfilling)
  t
}

plot_criterion = function(df, criterion, ignore_zero=FALSE, title='', y_label='Fulfilling criterion'){
  table = table_criterion(df, criterion, ignore_zero = ignore_zero) %>% group_by(label) %>% mutate(type=sse_category(label))
  g = ggplot(table) + 
    aes(x=label, y=100*fulfilling, fill=type) +
    geom_col() + 
    scale_fill_brewer(palette='Set1') + 
    coord_cartesian(ylim=c(0,100)) +
    labs(x='SSE', y=paste(y_label, ' [%]'), fill='Type', title=title) +
    order_x_labels()
  g
}

plot_criteria = function(df, ..., ignore_zero=FALSE, title='', y_label='Fulfilling criterion'){
  table = table_criteria(df, ..., ignore_zero = ignore_zero) %>% group_by(label, set) %>% mutate(type=sse_category(label))
  g = ggplot(table) + 
    aes(x=label, y=100*fulfilling, fill=set) +
    geom_col(position='dodge') + 
    scale_fill_brewer(palette='Set1') + 
    coord_cartesian(ylim=c(0,100)) +
    labs(x='SSE', y=paste(y_label, ' [%]'), fill='Type', title=title) +
    order_x_labels()
  g
}

plot_contained_helix_types = function(df, ...){
  #type_palette = c("#F89850", "#43B03C", "#9048A0", "#CCCCCC") # brighter green, darker violet
  #type_palette = c("#F89850", "#33A02C", "#B86EC3", "#CCCCCC") # darker green, brighter violet
  type_palette = c("#B080C8", "#30A028", "#F8A870", "#CCCCCC") # darker green, brighter violet
  changed_labels = c(contains_G = expression(paste('contains ', '3'[10])), contains_Η = expression(paste('contains ', alpha)), contains_Ι = expression(paste('contains ', pi)) )
  #changed_labels = c(contains_G = 'Contain Geee', contains_H = 'Contain α', contains_I = 'Contain I' )
  g = plot_criteria(df, contains_G=(df$longest_G>=3), contains_Η=(df$longest_H>=4), contains_Ι=(df$longest_I>=5), ..., ignore_zero=TRUE)  + 
    scale_fill_manual(values=type_palette, labels=changed_labels) + 
    order_x_labels_helices() + 
    theme(legend.text.align = 0)
  g
}

plot_sse_occurrence = function(df, show_confidence=TRUE, alpha=0.05, title=NULL, stagger_labels=FALSE, turn_labels=FALSE) {
  g = ggplot(table_sse_occurrence(df, alpha)) +
    aes(x=label, y=100*freq, fill=type) +
    geom_col() + 
    geom_vline(xintercept=seq(1.5, length(unique(df$label))-0.5, 1), lwd=0.4, colour='gray90') +
    scale_fill_brewer(palette='Set1') + 
    coord_cartesian(ylim=c(0, 105)) +
    scale_y_continuous(expand=expand_scale(add=0), breaks=seq(0,100,by=20)) +
    labs(x='SSE', y='Occurrence [%]', fill='Type', title=title) +
    order_x_labels(stagger=stagger_labels) +
    theme_bw() +
    customize_theme(stagger_labels=stagger_labels, turn_labels=turn_labels)
  if (show_confidence) { g = g + geom_errorbar(aes(ymin=100*freq_min, ymax=100*freq_max), width=0.4, position=position_dodge(.9)) }
  g
}

plot_bulge_occurrence = function(df, bulge_df, show_confidence=TRUE, alpha=0.05, title='', include_parent_strands=TRUE) {
  g = ggplot(table_bulge_occurrence(df, bulge_df, alpha=alpha, include_parent_strands=include_parent_strands)) +
    aes(x=label, y=100*freq, fill=bulge) +
    geom_col() + 
    scale_fill_brewer(palette='Paired') + 
    #coord_cartesian(ylim=c(0,100)) +
    labs(x='SSE', y='Occurrence [%]', fill='Bulge type', title=title)
  if (show_confidence) { g = g + geom_errorbar(aes(ymin=100*freq_min, ymax=100*freq_max), width=0.4, position=position_dodge(.9)) }
  g
}

plot_sse_occurrence_multi = function(..., show_confidence = TRUE, alpha = 0.05, title = NULL, stagger_labels = FALSE, turn_labels = FALSE) {
  tables = lapply(pairlist(...), table_sse_occurrence, alpha)
  table = bind_rows(tables, .id = 'set')
  g = ggplot(table) +
    aes(x = label, y = 100*freq, fill = paste(set, ' (', type, ')', sep = '')) +
    geom_col(position = 'dodge') + 
    geom_vline(xintercept = seq(1.5, length(unique(table$label))-0.5, 1), lwd = 0.4, colour = 'gray90') +
    scale_fill_manual(values = nonretarded_paired_palette_separated) + 
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

plot_compare_sse_lengths = function(..., ignore_zero=FALSE, set_names=NULL) {
  ts = lapply(list(...), table_sse_lengths, ignore_zero=ignore_zero)
  t = bind_rows(ts, .id='set')
  if (!is.null(set_names)) t$set = factor(t$set, labels=set_names)
  t$group = paste(t$type,' (',t$set,')',sep='')
  ggplot(t) +
    aes(x=label,y=avg_len,ymin=avg_len-sd_len,ymax=avg_len+sd_len,fill=group) +
    geom_col(position='dodge') + 
    #scale_fill_brewer(palette='Paired') + 
    scale_fill_manual(values=nonretarded_paired_palette) + 
    #geom_point() + 
    geom_errorbar(width=0.9, position='dodge') +
    labs(x='SSE', y='Mean length [#residues]', fill='Type') + 
    order_x_labels()
}

add_predominant_type1 = function(df){
  df$predom_type = ifelse (df$longest_G >= 3 & df$bonds_G > df$bonds_H + df$bonds_I, "mainly G",
                   ifelse (df$longest_H >= 4 & df$bonds_H > df$bonds_G + df$bonds_I, "mainly H",
                   ifelse (df$longest_I >= 5 & df$bonds_I > df$bonds_G + df$bonds_H, "mainly I",
                   "other")))
  df
}

add_predominant_type3 = function(df){
  df$predom_type = ifelse (df$longest_G >= 3 & df$bonds_G > df$bonds_H & df$bonds_G > df$bonds_I, "mainly G",
                   ifelse (df$longest_H >= 4 & df$bonds_H > df$bonds_G & df$bonds_H > df$bonds_I, "mainly H",
                   ifelse (df$longest_I >= 5 & df$bonds_I > df$bonds_G & df$bonds_I > df$bonds_H, "mainly I",
                   ifelse (df$longest_G < 3 & df$longest_H < 4 & df$longest_I < 5, "shorty",
                   ifelse (df$bonds_G == df$bonds_H, "mongrel",
                   "other")))))
  df
}

add_predominant_type4 = function(df){
  df$predom_type = ifelse (df$bonds_G > df$bonds_H & df$bonds_G > df$bonds_I, "mainly G",
                   ifelse (df$bonds_H > df$bonds_G & df$bonds_H > df$bonds_I, "mainly H",
                   ifelse (df$bonds_I > df$bonds_G & df$bonds_I > df$bonds_H, "mainly I",
                   "other")))
  df
}

add_predominant_type = function(df) {add_predominant_type4(df)}

table_predominant_type = function(df){
  relevant = df %>% add_predominant_type() %>% filter(label %in% HELIX_ORDER & length > 0)
  total = relevant  %>% group_by(label) %>% count() %>% transmute(total=n)
  by_type = relevant  %>% group_by(label, predom_type) %>% count()
  result = full_join(by_type, total, by='label') %>% transmute(freq=n/total)
  result
}

plot_predominant_type = function(df, title=''){
  table = table_predominant_type(df)
  table$predom_type = factor(table$predom_type, levels=c('mainly H', 'other', 'mainly G', 'mainly I'))
  #type_palette = c("#E88840", "#CCCCCC", "#33A02C", "#A6CEE3", "#1F78B4", "#B2DF8A")
  type_palette = c("#30A028", "#C8C8C8", "#B080C8", "#F8A870", "#602000", "#A6CEE3", "#1F78B4", "#B2DF8A")
  # red from Set1: "#E41A1C", blue from Set1: "#377EB8", blue from Set2: "#8DA0CB"
  # orange from Set1: "#FF7F00", orange from Set2: "#FC8D62", orange from Set3: "#FDB462"
  changed_labels = c('mainly H' = expression(paste('mainly ', alpha)), 'mainly I' = expression(paste('mainly ', pi)), 'mainly G' = expression(paste('mainly ', 3[10])) )
  g = ggplot(table) + 
    aes(x=label, y=100*freq, fill=predom_type) + 
    geom_col() +
    #scale_fill_brewer(palette='Set2') +
    scale_fill_manual(values=type_palette, labels=changed_labels) +
    coord_cartesian(ylim=c(0,100)) +
    labs(x='SSE', y='Fraction [%]', fill='Type', title=title) +
    order_x_labels_helices() + 
    theme(legend.text.align = 0)
  g
}

hist_sse_length = function(df, sse_label){
  data = df %>% filter(label==sse_label)
  ggplot(data) + aes(x=length) + geom_bar() + labs(title = sse_label) + xlim(0,NA)
}

generate_hists_sse_length = function(df, outfile_prefix='plots/length-histogram-'){
  labels = (df %>% distinct(label)) $ label
  for (l in labels) {
    filename = paste(outfile_prefix, l, '.png', sep='')
    hist_sse_length(df, l)
    ggsave(filename, width=8, height=4)
  }
}

spread_lengths = function(df) {
  df[,c("UniProt","PDB","chain_id","label","length")]%>%spread(label,length)
}

spread_coords = function(df) {
  df[,c("UniProt","PDB","chain_id","label","x1","x2","y1","y2","z1","z2")] %>% 
    gather(x1,x2,y1,y2,z1,z2,key='coord',value='coordval') %>% 
    mutate(coord=paste(label,coord,sep='_')) %>% 
    select(-label) %>%
    spread(coord, coordval)
}

replace_nas = function(df, new_value) {
  df %>% mutate_all( funs(replace(., is.na(.), new_value)) )
}

pca_labels = function(spread_df) {
  pca=prcomp(spread_df[,c(-1,-2,-3)], center=TRUE, scale. = FALSE)
  pca_df = data.frame(pca$rotation)
  ggplot(pca_df) + 
    aes(x=PC1,y=PC2,color=PC3,label=rownames(pca_df)) + 
    scale_color_gradient2(low='green',mid='black',high='red',midpoint=0.5*(min(pca_df$PC3)+max(pca_df$PC3))) + 
    geom_label()
  #ggsave('R_plot.png',width=12,height=8)
}

pca_pdbs = function(spread_df) {
  pca=prcomp(spread_df[,c(-1,-2,-3)], center=TRUE, scale. = FALSE)
  pca_df = data.frame(pca$x)
  rownames(pca_df) = spread_df$PDB
  ggplot(pca_df) + 
    aes(x=PC1,y=PC2,color=PC3,label=rownames(pca_df)) + 
    scale_color_gradient2(low='green',mid='black',high='red',midpoint=0.5*(min(pca_df$PC3)+max(pca_df$PC3))) + 
    geom_label()
  ggsave('R_plot.png',width=12,height=8)
  pca
}

pca_sse_coords = function(df, pc_x='PC1',pc_y='PC2', show_lines=FALSE, use_ggsave=TRUE) {
  df = df %>% filter(!is.na(x1))
  pca=prcomp(df[,c('x1','y1','z1','x2','y2','z2')], center=TRUE, scale. = FALSE)
  pca_df = data.frame(pca$x)
  rownames(pca_df) = make.names(paste(df$label, df$PDB, sep='@'), unique=TRUE)
  pca_df$label = df$label
  pca_df$PDB = df$PDB
  #pca_df_distinct = pca_df %>% distinct(label, .keep_all=TRUE)
  pca_df_distinct = pca_df %>% group_by(label) %>% summarise_at(.vars=c(pc_x,pc_y), .funs=mean) %>% ungroup() %>% arrange(match(label,sse_order))
  #View(pca_df_distinct)
  fucking_colors = setNames(rep(brewer.pal(12,"Paired"),times=100)[1:length(sse_order)], sse_order)
  fucking_colors['h'] = 'black'
  fucking_colors['e'] = 'black'
  
  ggplot(pca_df) + 
    aes(x=pca_df[,pc_x], y=pca_df[,pc_y], color=label,label=rownames(pca_df)) + 
    #scale_color_distiller(palette='Paired') + 
    #scale_colour_manual(values=rep(brewer.pal(12,"Paired"),times=100)) + #brewer.pal(8,"Set2")
    scale_color_manual(values=fucking_colors) + #brewer.pal(8,"Set2")
    geom_point() + 
    (if (show_lines) geom_path(data=pca_df_distinct, aes(x=pca_df_distinct[,pc_x], y=pca_df_distinct[,pc_y], label=label, group=1))) + 
    geom_label(data=pca_df_distinct, aes(x=pca_df_distinct[,pc_x], y=pca_df_distinct[,pc_y], label=label)) + 
    labs(x=pc_x, y=pc_y)
  if (use_ggsave) ggsave('R_plot.png',width=12,height=8) 
  pca
}

two_sample_test_by_labels = function(test, df1, df2, stat_col, label_col='label', p_limit=1) {
  result = unique(df1[label_col]) %>% mutate(p.value = NaN)
  for (label in result[[label_col]]) {
    sample1 = df1[which(df1[label_col]==label), stat_col]
    sample2 = df2[which(df2[label_col]==label), stat_col]
    p = test(sample1, sample2)$p.value
    result[which(result[label_col]==label), 'p.value'] = p 
  }
  result[which(result$p.value<=p_limit),]
}

two_sample_test_by_labels_with_comparison_bonf = function(test, df1, df2, stat_col, label_col='label', p_limit=1, print_all=FALSE, reversed=FALSE, signif_digits=2) {
  if (reversed){
    result = two_sample_test_by_labels_with_comparison_bonf(test, df2, df1, stat_col, label_col=label_col, p_limit=p_limit, print_all=print_all)
    result = result %>% mutate(mean_diff = -mean_diff)
  } else {
    result = unique(df1[label_col]) %>% mutate(p.value = NaN, comparison = '', mean_diff = 0)
    result = result %>% arrange(result[[label_col]])
    n = length(result[[label_col]])
    p_limit_bonf = p_limit / n
    for (label in result[[label_col]]) {
      sample1 = df1[which(df1[label_col]==label), stat_col]
      sample2 = df2[which(df2[label_col]==label), stat_col]
      p = test(sample1, sample2)$p.value
      pG = test(sample1, sample2, alternative='greater')$p.value
      pL = test(sample1, sample2, alternative='less')$p.value
      result[which(result[label_col]==label), 'p.value'] = signif(p, digits=signif_digits) 
      result[which(result[label_col]==label), 'comparison'] = if (p>p_limit | pG==pL) '   =   ' else if (pG<pL & p<=p_limit_bonf) '>>     ' else if (pG<pL) ' >     ' else if (pL<pG & p<=p_limit_bonf) '     <<' else if (pL<pG) '     < '
      m1 = mean(sample1)
      m2 = mean(sample2)
      result[which(result[label_col]==label), 'mean_diff'] = signif(m2 - m1, digits=signif_digits)
    }
    if (!print_all) {
      result = result %>% filter(p.value<=p_limit)
    }
  }
  result
  # explanation (> left set is stochastically greater, >> detto even after Bonferroni correction, == statistically insignificant)
}

two_sample_test_by_labels_with_comparison = function(test, df1, df2, stat_col, label_col='label', p_limit=1, print_all=FALSE, reversed=FALSE, signif_digits=2) {
  if (reversed){
    result = two_sample_test_by_labels_with_comparison(test, df2, df1, stat_col, label_col=label_col, p_limit=p_limit, print_all=print_all, signif_digits=signif_digits)
    result = result %>% mutate(mean_diff = -mean_diff, median_diff = -median_diff)
  } else {
    result = unique(df1[label_col]) %>% mutate(p.value = NaN, p.g = NaN, p.l = NaN, compare = '', mean_diff = 0)
    result = result %>% arrange(result[[label_col]])
    n = length(result[[label_col]])
    for (label in result[[label_col]]) {
      sample1 = df1[which(df1[label_col]==label), stat_col]
      sample2 = df2[which(df2[label_col]==label), stat_col]
      p = test(sample1, sample2)$p.value
      pG = test(sample1, sample2, alternative='greater')$p.value
      pL = test(sample1, sample2, alternative='less')$p.value
      result[which(result[label_col]==label), 'p.value'] = signif(p, digits=signif_digits) 
      result[which(result[label_col]==label), 'p.g'] = signif(pG, digits=signif_digits) 
      result[which(result[label_col]==label), 'p.l'] = signif(pL, digits=signif_digits) 
      result[which(result[label_col]==label), 'compare'] = if (p>p_limit | pG==pL) '       ' else if (pG<pL) ' >     ' else if (pL<pG) '     < '
      m1 = mean(sample1)
      m2 = mean(sample2)
      med1 = median(sample1)
      med2 = median(sample2)
      result[which(result[label_col]==label), 'mean_diff'] = signif(m2 - m1, digits=signif_digits)
      result[which(result[label_col]==label), 'median_diff'] = signif(med2 - med1, digits=signif_digits)
    }
    if (!print_all) {
      result = result %>% filter(p.value<=p_limit)
    }
  }
  result
}

two_sample_occurrence_prop.test_bonf = function(df1, df2, p_limit=0.05, print_all=FALSE, signif_digits=2) {
  n_pdb1 = n_distinct(df1$PDB)
  data1 = df1 %>% filter(length>0) %>% group_by(label) %>% count() %>% 
    mutate(present1 = n, absent1 = n_pdb1-n, freq1 = n/n_pdb1) %>% select(-n)
  n_pdb2 = n_distinct(df2$PDB)
  data2 = df2 %>% filter(length>0) %>% group_by(label) %>% count() %>% 
    mutate(present2 = n, absent2 = n_pdb2-n, freq2 = n/n_pdb2) %>% select(-n)
  data = full_join(data1, data2, by='label') %>% mutate(p.value = 0, comparison = '', freq_diff = freq2-freq1)
  
  n_labels = length(data$label)
  p_limit_bonf = p_limit / n_labels
  #View(data)
  
  for (the_label in data$label) {
    row = filter(data, label == the_label)
    # View(row)
    test_result = prop.test(c(row$present1, row$present2), c(n_pdb1, n_pdb2))
    p = test_result$p.value
    freq_diff = row$freq_diff
    data[which(data$label == the_label), 'p.value'] = signif(p, digits=signif_digits) 
    data[which(data$label == the_label), 'comparison'] = if (p>p_limit | freq_diff==0) '       ' else if (freq_diff<0 & p<=p_limit_bonf) '>>     ' else if (freq_diff<0) ' >     ' else if (freq_diff>0 & p<=p_limit_bonf) '     <<' else if (freq_diff>0) '     < '
    data[which(data$label == the_label), 'freq_diff'] = signif(freq_diff, digits=signif_digits) 
  }
  if (!print_all) {
    data = data %>% filter(p.value<=p_limit)
  }
  data %>% select(label, p.value, comparison, freq_diff) %>% data.frame()
  # explanation (> left set is stochastically greater, >> detto even after Bonferroni correction, == statistically insignificant)
}

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

two_sample_occurrence_fisher.test = function(df1, df2, p_limit=0.05, print_all=FALSE, signif_digits=2) {
  n_pdb1 = n_distinct(df1$PDB)
  data1 = df1 %>% filter(length>0) %>% group_by(label) %>% count() %>% 
    mutate(present1 = n, absent1 = n_pdb1-n, freq1 = n/n_pdb1) %>% select(-n)
  n_pdb2 = n_distinct(df2$PDB)
  data2 = df2 %>% filter(length>0) %>% group_by(label) %>% count() %>% 
    mutate(present2 = n, absent2 = n_pdb2-n, freq2 = n/n_pdb2) %>% select(-n)
  data = full_join(data1, data2, by='label') %>% mutate(p.value = 0, comparison = '', freq_diff = freq2-freq1)
  
  n_labels = length(data$label)
  p_limit_bonf = p_limit / n_labels
  #View(data)
  
  for (the_label in data$label) {
    row = filter(data, label == the_label)
    # View(row)
    the_matrix = matrix(c(row$present1, row$absent1, row$present2, row$absent2), ncol=2)
    test_result = fisher.test(the_matrix)
    p = test_result$p.value
    freq_diff = row$freq_diff
    data[which(data$label == the_label), 'p.value'] = signif(p, digits=signif_digits) 
    data[which(data$label == the_label), 'comparison'] = if (p>p_limit | freq_diff==0) '   =   ' else if (freq_diff<0 & p<=p_limit_bonf) '>>     ' else if (freq_diff<0) ' >     ' else if (freq_diff>0 & p<=p_limit_bonf) '     <<' else if (freq_diff>0) '     < '
    data[which(data$label == the_label), 'freq_diff'] = signif(freq_diff, digits=signif_digits) 
  }
  if (!print_all) {
    data = data %>% filter(p.value<=p_limit)
  }
  data %>% select(label, p.value, comparison, freq_diff) %>% data.frame()
  # explanation (> left set is stochastically greater, >> detto even after Bonferroni correction, == statistically insignificant)
}

agresti_coull_confidence_interval = function(successes, total, alpha=0.05) {
  z = qnorm(1-alpha/2)
  n_ = total + z^2
  p_ = 1/n_ * (successes + z^2/2)
  d = z * sqrt(p_ * (1-p_) / n_)
  c(p_-d, p_+d)
}

table_bulges_with_missing_with_ligands = function(df, bulge_df, ligads_df, only_UniProts_with_more_PDBs=FALSE) {
  all_PDBs_UniProts = df[c('UniProt', 'PDB')] %>% distinct(PDB, .keep_all=TRUE) 
  bulge_df = bulge_df %>% select(- UniProt)
  if (only_UniProts_with_more_PDBs){
    selected_UniProts = (all_PDBs_UniProts %>% group_by(UniProt) %>% count() %>% filter(n > 1))['UniProt']
    all_PDBs_UniProts = left_join(selected_UniProts, all_PDBs_UniProts, by='UniProt')
  }
  bulges_with_missing = left_join(all_PDBs_UniProts, bulge_df, by='PDB')
  bulges_with_missing_with_ligands = left_join(bulges_with_missing, ligands, by='PDB')
  bulges_with_missing_with_ligands
}

compare_annotations = function(dfX, dfY){
  oldies = anti_join(dfX %>% select(Domain, label, start, end, type, length), dfY %>% select(Domain, label, start, end, type, length))
  newbies = anti_join(dfY %>% select(Domain, label, start, end, type, length), dfX %>% select(Domain, label, start, end, type, length))
  diff = full_join(oldies, newbies, by=c('Domain', 'label'))
  diff
}

print_png = function(filename, width = 2000, ratio = 2/1, height = NULL, res = 250) { 
  if (is.null(height)) { height = width/ratio }
  dev.print(png, filename, width = width, height = height, res = res)
}

print_tif = function(filename, width = 2000, ratio = 2/1, height = NULL, res = 250) { 
  if (is.null(height)) { height = width/ratio }
  dev.print(tiff, filename, width = width, height = height, res = res)
}

