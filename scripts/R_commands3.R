# Comment: Ctrl+Shift+C
# Go to function definition: F2 / Ctrl+Click
# Fold all: Alt+O

DATADIR = '/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/testing_20200128-verbose'
source('R_lib.R')
dir.create(full_path('plots'), showWarnings = FALSE)


# READ DATASETS
taxons = read.csv(full_path('domain_taxons_groups.tsv'), sep='\t', header = FALSE, col.names=c('Domain','TaxID','Group'))
setNR = read_tsv(full_path('annotations_with_reference_residues_NR.tsv')) %>%
  left_join(taxons, by = 'Domain') %>% select(-starts_with('longest_'), -starts_with('bonds_')) %>% filter(label %in% OUR_SSES)
setALL = read_tsv(full_path('annotations_with_reference_residues_ALL.tsv')) %>% 
  left_join(taxons, by = 'Domain') %>% select(-starts_with('longest_'), -starts_with('bonds_')) %>% filter(label %in% OUR_SSES) %>% keep_one_domain_per_pdb()
# ligands = read.csv(full_path('ligands.tsv'), sep='\t', header=FALSE, col.names=c('PDB','Ligands'))

setNR_Bact = filter(setNR, Group=='Bact')
setNR_Euka = filter(setNR, Group=='Euka')
# setCAM = setALL %>% filter(UniProt=='P00183')
# setBM3 = setALL %>% filter(UniProt=='P14779')
# set3A4 = setALL %>% filter(UniProt=='P08684')

domainsNR = get_domains(setNR)
domainsALL = get_domains(setALL, summarize_by_UniProt = TRUE)
domains_combined = full_join(domainsNR, domainsALL, by = 'UniProt') %>% left_join(taxons) %>% 
  transmute(Group = Group, UniProt = UniProt, SetNR = Domain, SetALL_count = Count, SetALL = Domains) %>% arrange(Group)
write_tsv(domains_combined, full_path('plots/domain_lists_table.tsv'))

table((setNR %>% distinct(PDB, Group))$Group)
barplot(table((setNR %>% distinct(PDB, Group))$Group), main = 'Number of PDB entries in superkingdoms (Set-NR)')


# PLOTS FOR OCCURRENCE
occurrence_table_NR = table_sse_occurrence(setNR, alpha = 0.05)
# plot_sse_occurrence(setNR, show_confidence = TRUE, alpha = 0.05)
# print_png(full_path('plots/occurrence-setNR.png'), width = 4000, ratio = 2/1, res = 400)
plot_sse_occurrence(setNR, show_confidence = TRUE, alpha = 0.05, turn_labels = TRUE)
print_png(full_path('plots/occurrence-setNR-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/occurrence-setNR-500t.tif'), width = 4000, ratio = 2/1, res = 500)
# plot_sse_occurrence(setNR, show_confidence = TRUE, alpha = 0.05, stagger_labels = TRUE)
# print_png(full_path('plots/occurrence-setNR-500s.png'), width = 4000, ratio = 2/1, res = 500)

# plot_sse_occurrence_multi(Bact = setNR_Bact, Euka = setNR_Euka)
# print_png(full_path('plots/occurrence-setNR-Bact-Euka.png'), width = 4000, ratio = 2/1, res = 400)
plot_sse_occurrence_multi(Bact = setNR_Bact, Euka = setNR_Euka, turn_labels = TRUE)
print_png(full_path('plots/occurrence-setNR-Bact-Euka-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/occurrence-setNR-Bact-Euka-500t.tif'), width = 4000, ratio = 2/1, res = 500)
# plot_sse_occurrence_multi(Bact = setNR_Bact, Euka = setNR_Euka, stagger_labels = TRUE)
# print_png(full_path('plots/occurrence-setNR-Bact-Euka-500s.png'), width = 4000, ratio = 2/1, res = 500)


# PLOTS FOR LENGTH DISTRIBUTION
boxplot_sse(setNR, ignore_zero = TRUE, title = 'Set-NR')

# violinplot_sse(setNR, ignore_zero = TRUE)
# print_png(full_path('plots/length-setNR.png'), width = 4000, ratio = 2/1, res = 400)
violinplot_sse(setNR, ignore_zero = TRUE, turn_labels = TRUE)
print_png(full_path('plots/length-setNR-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/length-setNR-500t.tif'), width = 4000, ratio = 2/1, res = 500)
# violinplot_sse(setNR, ignore_zero = TRUE, stagger_labels = TRUE)
# print_png(full_path('plots/length-setNR-500s.png'), width = 4000, ratio = 2/1, res = 500)

# violinplot_sse_multi(Bact = setNR_Bact, Euka = setNR_Euka, ignore_zero = TRUE)
# print_png(full_path('plots/length-setNR-Bact-Euka.png'), width = 4000, ratio = 2/1, res = 400)
violinplot_sse_multi(Bact = setNR_Bact, Euka = setNR_Euka, ignore_zero = TRUE, turn_labels = TRUE) 
print_png(full_path('plots/length-setNR-Bact-Euka-500t.png'), width = 4000, ratio = 2/1, res = 500)
print_tif(full_path('plots/length-setNR-Bact-Euka-500t.tif'), width = 4000, ratio = 2/1, res = 500)
# violinplot_sse_multi(Bact = setNR_Bact, Euka = setNR_Euka, ignore_zero = TRUE, stagger_labels = TRUE)
# print_png(full_path('plots/length-setNR-Bact-Euka-500s.png'), width = 4000, ratio = 2/1, res = 500)
# point = mean, horizontal line = median


# STATISTICAL COMPARISON Set-NR-Bact vs. Set-NR-Euka
setNR_Bact_nonzero = filter(setNR_Bact, length > 0)
setNR_Euka_nonzero = filter(setNR_Euka, length > 0)

two_sample_occurrence_prop_test(setNR_Bact, setNR_Euka, p_limit = 0.05, print_all = TRUE)
two_sample_occurrence_fisher_test(setNR_Bact, setNR_Euka, p_limit = 0.05, print_all = TRUE)

two_sample_test_by_labels_with_comparison(ks.test, setNR_Bact_nonzero, setNR_Euka_nonzero, 'length', label_col='label', p_limit=0.05, print_all=TRUE, reversed=TRUE)
two_sample_test_by_labels_with_comparison(wilcox.test, setNR_Bact_nonzero, setNR_Euka_nonzero, 'length', label_col='label', p_limit=0.05, print_all=TRUE)
# alternative hypothesis for KS test is reversed ("greater" mean stochastically smaller)!



# # Removing outliers
# sNR_noout = remove_outliers(setNR, length, ignore_zero=TRUE)
# sNR_noout1 = remove_outliers(setNR, length, ignore_zero=TRUE, extra_tolerance=1)
# 
# two_sample_test_by_labels_with_comparison(ks.test, filter(setNR, length>0), sNR_noout, 'length', label_col='label', p_limit=0.05, print_all=TRUE, reversed=TRUE)
# two_sample_test_by_labels_with_comparison(wilcox.test, filter(setNR, length>0), sNR_noout, 'length', label_col='label', p_limit=0.05, print_all=TRUE)
# 
# two_sample_test_by_labels_with_comparison(ks.test, filter(sNR_noout, Group=='Bact'), filter(sNR_noout, Group=='Euka'), 'length', label_col='label', p_limit=0.05, print_all=TRUE, reversed=TRUE)
# two_sample_test_by_labels_with_comparison(wilcox.test, filter(sNR_noout, Group=='Bact'), filter(sNR_noout, Group=='Euka'), 'length', label_col='label', p_limit=0.05, print_all=TRUE)

# setNR - 4940 obs., setNR nonzero - 4161, noout1 - 4021 (140 outliers), noout - 3838 (323 outliers)


# # Bulges

# bNR = read_bulges_json_new(full_path('annotations_with_reference_residues_NR.json'), c('chain_id', 'start','end'))
# bNR = left_join(add_sse_lengths(bNR), taxons, by='Domain') %>% filter(label %in% OUR_SSES)
# bALL = read_bulges_json_new(full_path('annotations_with_reference_residues_ALL.json'), c('chain_id', 'start','end'), one_domain_per_pdb=TRUE)
# bALL = left_join(add_sse_lengths(bALL), taxons, by='Domain') %>% filter(label %in% OUR_SSES)

# View(table_bulge_occurrence(setNR, bNR))
# View(table_bulge_occurrence(setNR %>% filter(Group=='Bact'), bNR %>% filter(Group=='Bact')))
# View(table_bulge_occurrence(setNR %>% filter(Group=='Euka'), bNR %>% filter(Group=='Euka')))
# plot_bulge_occurrence(setNR, bNR, show_confidence = TRUE, alpha=0.05, title='Set-NR', include_parent_strands=FALSE)
# print_png(full_path('plots/bulge_occurrence-setNR.png'), width=4000, ratio=2/1, res=410)
# 
# violinplot_sse(bNR, 'length', ignore_zero = TRUE, title='Set-NR')
# violinplot_sse(setNR, 'longest_G', ignore_zero = TRUE, title='Set-NR')
# print_png('plots/violinplot.png', width=4000, height=2000)
# print_png('plots/violinplot-small.png', width=1800, height=900)
# print_png('plots/violinplot-small.png', width=1800, height=1800)
# 
# table_bulges_with_missing_with_ligands(setALL, bALL %>% filter(label=='4.2' & bulge=='N'), ligands, only_UniProts_with_more_PDBs=TRUE) %>% View()
# table_bulges_with_missing_with_ligands(setALL, bALL %>% filter(label=='3.3' & bulge=='N'), ligands, only_UniProts_with_more_PDBs=TRUE) %>% View()


# # Helix types
# helices = setNR %>% filter(label %in% helix_order & length > 0)
# 
# helices %>% add_predominant_type() %>% group_by(predom_type) %>% count() %>% View()
# pt = table_predominant_type(setNR)  # Requires nested_sses in helices (SecStrAnnotator --verbose)
# pt$predom_type = factor(pt$predom_type)
# pt %>% spread(predom_type, freq) %>% View()
# plot_predominant_type(setNR, title='Predominant helix type (Set-NR)')  # Requires nested_sses in helices (SecStrAnnotator --verbose)
# print_png(full_path('plots/predominant_type-setNR.png'), width=4000, ratio=2/1, res=450)
# 
# table_criteria(helices, have_G=(helices$longest_G>=3), have_H=(helices$longest_H>=4), have_I=(helices$longest_I>=5), by_label=FALSE) %>% spread(set, fulfilling) %>% View()
# table_criteria(helices, have_G=(helices$longest_G>=3), have_H=(helices$longest_H>=4), have_I=(helices$longest_I>=5)) %>% spread(set, fulfilling) %>% View()
# table_criteria(setNR_Bact, have_G=(setNR_Bact$longest_G>=3), have_H=(setNR_Bact$longest_H>=4), have_I=(setNR_Bact$longest_I>=5), ignore_zero=TRUE) %>% spread(set, fulfilling) %>% filter(label %in% helix_order) %>% View()
# table_criteria(setNR_Euka, have_G=(setNR_Euka$longest_G>=3), have_H=(setNR_Euka$longest_H>=4), have_I=(setNR_Euka$longest_I>=5), ignore_zero=TRUE) %>% spread(set, fulfilling) %>% filter(label %in% helix_order) %>% View()
# table_criteria(setNR, have_I=(setNR$longest_I>=5), start_with_I=(setNR$longest_I_start>=5), end_with_I=(setNR$longest_I_end>=5), ignore_zero=TRUE) %>% spread(set, fulfilling) %>% filter(label %in% helix_order) %>% View()
# 
# plot_contained_helix_types(setNR, title='Contained helix types (Set-NR)', y_label='Fraction')
# print_png(full_path('plots/contained_types-setNR.png'), width=4000, ratio=2/1, res=450)
# 
# 
# plot_sse_occurrence(setNR, group_col='predType', show_confidence = TRUE, alpha=0.05, title='Set-NR')
# plot_sse_occurrence(setNR %>% filter(longest_I>=5), show_confidence = TRUE, alpha=0.05, title='Set-NR, longest_I >= 5')
# jitterplot_xy(filter(setNR_Euka, label=="J'"), longest_G, longest_H)
