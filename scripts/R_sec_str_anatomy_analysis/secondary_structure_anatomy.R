# SCRIPT FOR ANALYSIS OF SECONDARY STRUCTURE ANATOMY OF A PROTEIN FAMILY (CYTOCHROMES P450)

source('secondary_structure_anatomy_lib.R')  # Contains a few CytochromeP450-specific constants!
DATADIR = '/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/CytochromesP450-20200323'


# READ DATASETS
taxons = read.csv(full_path('domain_taxons_groups.tsv'), sep='\t', header = FALSE, col.names=c('Domain','TaxID','Group'))
setNR = read_tsv(full_path('annotations_with_reference_residues_NR.tsv')) %>%
  left_join(taxons, by = 'Domain') %>% select(-starts_with('longest_'), -starts_with('bonds_')) %>% filter(label %in% OUR_SSES)
setALL = read_tsv(full_path('annotations_with_reference_residues_ALL.tsv')) %>% 
  left_join(taxons, by = 'Domain') %>% select(-starts_with('longest_'), -starts_with('bonds_')) %>% filter(label %in% OUR_SSES) %>% keep_one_domain_per_pdb()
dir.create(full_path('plots'), showWarnings = FALSE)

setNR_Bact = filter(setNR, Group=='Bact')
setNR_Euka = filter(setNR, Group=='Euka')

domainsNR = get_domains(setNR)
domainsALL = get_domains(setALL, summarize_by_UniProt = TRUE)
domains_combined = full_join(domainsNR, domainsALL, by = 'UniProt') %>% left_join(taxons) %>% 
  transmute(Group = Group, UniProt = UniProt, SetNR = Domain, SetALL_count = Count, SetALL = Domains) %>% arrange(Group)
write_tsv(domains_combined, full_path('domain_lists_table.tsv'))

table((setNR %>% distinct(PDB, Group))$Group)
barplot(table((setNR %>% distinct(PDB, Group))$Group), main = 'Number of PDB entries in superkingdoms (Set-NR)')


# PLOTS FOR OCCURRENCE
occurrence_table_NR = table_sse_occurrence(setNR, alpha = 0.05)
plot_sse_occurrence(setNR, show_confidence = TRUE, alpha = 0.05, turn_labels = TRUE)
print_png(full_path('plots/occurrence-setNR-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/occurrence-setNR-500t.tif'), width = 4000, ratio = 2/1, res = 500)

plot_sse_occurrence_multi(Bact = setNR_Bact, Euka = setNR_Euka, turn_labels = TRUE)
print_png(full_path('plots/occurrence-setNR-Bact-Euka-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/occurrence-setNR-Bact-Euka-500t.tif'), width = 4000, ratio = 2/1, res = 500)


# PLOTS FOR LENGTH DISTRIBUTION
boxplot_sse(setNR, ignore_zero = TRUE, title = 'Set-NR')

violinplot_sse(setNR, ignore_zero = TRUE, turn_labels = TRUE)
print_png(full_path('plots/length-setNR-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/length-setNR-500t.tif'), width = 4000, ratio = 2/1, res = 500)

violinplot_sse_multi(Bact = setNR_Bact, Euka = setNR_Euka, ignore_zero = TRUE, turn_labels = TRUE) 
print_png(full_path('plots/length-setNR-Bact-Euka-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/length-setNR-Bact-Euka-500t.tif'), width = 4000, ratio = 2/1, res = 500)
# point = mean, horizontal line = median


# # PLOTS FOR ONE ONE HELIX (EXAMPLE FOR GRAPHICAL ABSTRACT)
# plot_sse_occurrence(setNR, show_confidence = TRUE, alpha = 0.05) + scale_x_discrete(limits="J") + theme(legend.position = 'none')
# print_png(full_path('plots/occurrence-setNR-J.png'), width = 450, height = 1200, res = 500)
# plot_sse_occurrence_multi(Bact = setNR_Bact, Euka = setNR_Euka) + scale_x_discrete(limits="J") + theme(legend.position = 'none')
# print_png(full_path('plots/occurrence-setNR-Bact-Euka-J.png'), width = 450, height = 1200, res = 500)
# violinplot_sse(setNR, ignore_zero = TRUE) + scale_x_discrete(limits="J") + theme(legend.position = 'none') + coord_cartesian(ylim = c(0,25))
# print_png(full_path('plots/length-setNR-J.png'), width = 450, height = 1200, res = 500)
# violinplot_sse_multi(Bact = setNR_Bact, Euka = setNR_Euka, ignore_zero = TRUE) + scale_x_discrete(limits="J") + theme(legend.position = 'none') +
#   scale_color_manual(values = NONRETARDED_PAIRED_PALETTE_DARK) + coord_cartesian(ylim = c(0,25))
# print_png(full_path('plots/length-setNR-Bact-Euka-J.png'), width = 450, height = 1200, res = 500)


# STATISTICAL COMPARISON Set-NR-Bact vs. Set-NR-Euka
setNR_Bact_nonzero = filter(setNR_Bact, length > 0)
setNR_Euka_nonzero = filter(setNR_Euka, length > 0)

two_sample_occurrence_prop_test(setNR_Bact, setNR_Euka, p_limit = 0.05, print_all = TRUE)
two_sample_occurrence_fisher_test(setNR_Bact, setNR_Euka, p_limit = 0.05, print_all = TRUE)

two_sample_test_by_labels_with_comparison(ks.test, setNR_Bact_nonzero, setNR_Euka_nonzero, 'length', label_col='label', p_limit=0.05, print_all=TRUE, reversed=TRUE)
two_sample_test_by_labels_with_comparison(wilcox.test, setNR_Bact_nonzero, setNR_Euka_nonzero, 'length', label_col='label', p_limit=0.05, print_all=TRUE)
# alternative hypothesis for KS test is reversed ("greater" mean stochastically smaller)!


# BETA-BULGES - will only work if SecStrAnnotator was run with option  --verbose
bulgesNR = read_tsv(full_path('beta_bulges_NR.tsv')) %>%
  left_join(taxons, by = 'Domain') %>% filter(label %in% OUR_SSES) %>% mutate(bulge_label = paste(label, type, sep = '.'))
bulgesALL = read_tsv(full_path('beta_bulges_ALL.tsv')) %>%
  left_join(taxons, by = 'Domain') %>% filter(label %in% OUR_SSES) %>% mutate(bulge_label = paste(label, type, sep = '.')) %>% keep_one_domain_per_pdb()

View(table_bulge_occurrence(setNR, bulgesNR))
View(table_bulge_occurrence(setNR %>% filter(Group=='Bact'), bulgesNR %>% filter(Group=='Bact')))
View(table_bulge_occurrence(setNR %>% filter(Group=='Euka'), bulgesNR %>% filter(Group=='Euka')))
plot_bulge_occurrence(setNR, bulgesNR, show_confidence = TRUE, alpha = 0.05, include_parent_strands = FALSE, turn_labels = TRUE)
print_png(full_path('plots/bulge_occurrence-setNR.png'), width=4000, ratio=2/1, res=500)


# HELIX TYPES - will only work if SecStrAnnotator was run with option  --verbose
helicesNR = read_tsv(full_path('annotations_with_reference_residues_NR.tsv')) %>%
  left_join(taxons, by = 'Domain') %>% filter(label %in% HELIX_ORDER & length > 0)

plot_contained_helix_types(helicesNR, y_label = 'Fraction')
print_png(full_path('plots/contained_types-setNR.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/contained_types-setNR.tif'), width = 4000, ratio = 2/1, res = 500)


# Comment: Ctrl+Shift+C
# Go to function definition: F2 / Ctrl+Click
# Fold all: Alt+O
