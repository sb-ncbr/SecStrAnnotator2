# Comment: Ctrl+Shift+C
# Go to function definition: F2 / Ctrl+Click
# Fold all: Alt+O

DATADIR = '/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/testing_20200128-verbose'
source('R_lib.R')
dir.create(fullname('plots'), showWarnings=FALSE)

# Read  SSEs
distinction = read.csv(fullname('domain_taxons_groups.tsv'), sep='\t', header=FALSE, col.names=c('Domain','TaxID','Group'))
# ligands = read.csv(fullname('ligands.tsv'), sep='\t', header=FALSE, col.names=c('PDB','Ligands'))
sNR = read_sses_json_new(fullname('annotations_with_reference_residues_NR.json'), c('start','end','type'))
sNR = left_join(add_lengths(sNR), distinction, by='Domain') %>% filter(label %in% our_sses)
sNR = sNR %>% select(-starts_with("longest_"), -starts_with("bonds_"))
write.csv(sNR, fullname('SetNR.csv'), row.names=FALSE)
sALL = read_sses_json_new(fullname('annotations_with_reference_residues_ALL.json'), c('start','end','type'), one_domain_per_pdb=TRUE)
sALL = left_join(add_lengths(sALL), distinction, by='Domain') %>% filter(label %in% our_sses)

# sNR = read.csv(fullname('SetNR.csv'))

dALL = read_sses_json_new(fullname('annotations_with_reference_residues_ALL.json'), c('start'), one_domain_per_pdb=TRUE, only_domain_lists=TRUE)
dNR = read_sses_json_new(fullname('annotations_with_reference_residues_NR.json'), c('start'), one_domain_per_pdb=TRUE, only_domain_lists=TRUE)
dCombined = dNR %>%
  transmute(UniProt=UniProt, Domain=Domains) %>%
  left_join(distinction,by="Domain") %>%
  left_join(dALL, by="UniProt") %>%
  transmute(Group=Group, UniProt=UniProt, SetNR=Domain, CountInSetALL=Count, SetALL=Domains) %>%
  arrange(Group, UniProt)
write.table(dCombined, file=fullname('plots/domain_lists_table.tsv'), quote=FALSE, sep='\t', row.names=FALSE)

# barplot(table((sALL %>% distinct(PDB, Group))$Group), main='Number of PDB entries in superkingdoms (Set-ALL)')
barplot(table((sNR %>% distinct(PDB, Group))$Group), main='Number of PDB entries in superkingdoms (Set-NR)')
table((sNR %>% distinct(PDB, Group))$Group)


# sCAM = sALL %>% filter(UniProt=='P00183')
# sBM3 = sALL %>% filter(UniProt=='P14779')
# s3A4 = sALL %>% filter(UniProt=='P08684')

# # Compare results obtained by different methods
# diff = compare_annotations(sNR2, sNR2s)
# diff = compare_annotations(sNR2s, sNR2se)
# diff = compare_annotations(sNRsw, sNRs)
# diff = compare_annotations(sNRo, sNR)
# diff$label = factor(diff$label)

# the_label = "B'"
# View(diff_    %>% filter(label==the_label))
# diff_ = diff_ %>% filter(label!=the_label)


# # Read bulges
# bNR = read_bulges_json_new(fullname('annotations_with_reference_residues_NR.json'), c('chain_id', 'start','end'))
# bNR = left_join(add_lengths(bNR), distinction, by='Domain') %>% filter(label %in% our_sses)
# bALL = read_bulges_json_new(fullname('annotations_with_reference_residues_ALL.json'), c('chain_id', 'start','end'), one_domain_per_pdb=TRUE)
# bALL = left_join(add_lengths(bALL), distinction, by='Domain') %>% filter(label %in% our_sses)

# # Statistical comparison Set-NR vs. Set-ALL
# p_limit = 0.05
# my_filter = function(data){data} # function(data){ filter(data, Group=='Bact') }
# two_sample_test_by_labels_with_comparison(ks.test, my_filter(sNR), my_filter(sALL), 'length', label_col='label', p_limit=p_limit, print_all=TRUE, reversed=TRUE)
# two_sample_test_by_labels_with_comparison(wilcox.test, my_filter(sNR), my_filter(sALL), 'length', label_col='label', p_limit=p_limit, print_all=TRUE, signif_digits = 4)
# # alternative hypothesis for KS test is reversed ("greater" mean stochastically smaller)!


# Statistical comparison Set-NR-Bact vs. Set-NR-Euka
sNR_Bact = filter(sNR, Group=='Bact')
sNR_Euka = filter(sNR, Group=='Euka')
sNR_Bact_nz = filter(sNR_Bact, length>0)
sNR_Euka_nz = filter(sNR_Euka, length>0)

two_sample_occurrence_prop.test(sNR_Bact, sNR_Euka, p_limit=0.05, print_all=TRUE)
two_sample_occurrence_fisher.test(sNR_Bact, sNR_Euka, p_limit=0.05, print_all=TRUE)

two_sample_test_by_labels_with_comparison(ks.test, sNR_Bact_nz, sNR_Euka_nz, 'length', label_col='label', p_limit=0.05, print_all=TRUE, reversed=TRUE)
two_sample_test_by_labels_with_comparison(wilcox.test, sNR_Bact_nz, sNR_Euka_nz, 'length', label_col='label', p_limit=0.05, print_all=TRUE)
# alternative hypothesis for KS test is reversed ("greater" mean stochastically smaller)!


# Plots for occurrence
occurrence_table_NR = table_sse_occurrence(sNR, alpha=0.05)
# plot_sse_occurrence(sNR, show_confidence = TRUE, alpha=0.05)
# print_png(fullname('plots/occurrence-sNR.png'), width=4000, ratio=2/1, res=400)
plot_sse_occurrence(sNR, show_confidence = TRUE, alpha=0.05, turn_labels=TRUE)
print_png(fullname('plots/occurrence-sNR-500t.png'), width=4000, ratio=2/1, res=500)
print_tif(fullname('plots/occurrence-sNR-500t.tif'), width=4000, ratio=2/1, res=500)
# plot_sse_occurrence(sNR, show_confidence = TRUE, alpha=0.05, stagger_labels=TRUE)
# print_png(fullname('plots/occurrence-sNR-500s.png'), width=4000, ratio=2/1, res=500)

# plot_sse_occurrence(sCAM, show_confidence = TRUE, alpha=0.05, title='Set-CAM')
# print_png(fullname('plots/occurrence-sCAM.png'), width=4000, ratio=2/1, res=400)
# plot_sse_occurrence(sBM3, show_confidence = TRUE, alpha=0.05, title='Set-BM3')
# print_png(fullname('plots/occurrence-sBM3.png'), width=4000, ratio=2/1, res=400)
# plot_sse_occurrence(s3A4, show_confidence = TRUE, alpha=0.05, title='Set-3A4')
# print_png(fullname('plots/occurrence-s3A4.png'), width=4000, ratio=2/1, res=400)

# plot_compare_sse_occurrence('Set-NR'=sNR, 'Set-ALL'=sALL, show_confidence = TRUE, alpha=0.05, title='Set-ALL vs Set-NR')
# plot_compare_sse_occurrence(Bact=sNR_Bact, Euka=sNR_Euka)
# print_png(fullname('plots/occurrence-sNR-Bact-Euka.png'), width=4000, ratio=2/1, res=400)
plot_compare_sse_occurrence(Bact=sNR_Bact, Euka=sNR_Euka, turn_labels=TRUE)
print_png(fullname('plots/occurrence-sNR-Bact-Euka-500t.png'), width=4000, ratio=2/1, res=500)
print_tif(fullname('plots/occurrence-sNR-Bact-Euka-500t.tif'), width=4000, ratio=2/1, res=500)
# plot_compare_sse_occurrence(Bact=sNR_Bact, Euka=sNR_Euka, stagger_labels=TRUE)
# print_png(fullname('plots/occurrence-sNR-Bact-Euka-500s.png'), width=4000, ratio=2/1, res=500)


# # Removing outliers
# sNR_noout = remove_outliers(sNR, length, ignore_zero=TRUE)
# sNR_noout1 = remove_outliers(sNR, length, ignore_zero=TRUE, extra_tolerance=1)
# 
# two_sample_test_by_labels_with_comparison(ks.test, filter(sNR, length>0), sNR_noout, 'length', label_col='label', p_limit=0.05, print_all=TRUE, reversed=TRUE)
# two_sample_test_by_labels_with_comparison(wilcox.test, filter(sNR, length>0), sNR_noout, 'length', label_col='label', p_limit=0.05, print_all=TRUE)
# 
# two_sample_test_by_labels_with_comparison(ks.test, filter(sNR_noout, Group=='Bact'), filter(sNR_noout, Group=='Euka'), 'length', label_col='label', p_limit=0.05, print_all=TRUE, reversed=TRUE)
# two_sample_test_by_labels_with_comparison(wilcox.test, filter(sNR_noout, Group=='Bact'), filter(sNR_noout, Group=='Euka'), 'length', label_col='label', p_limit=0.05, print_all=TRUE)

# sNR - 4940 obs., sNR nonzero - 4161, noout1 - 4021 (140 outliers), noout - 3838 (323 outliers)

# Plots for length distribution
boxplot_sse(sNR, y_column='length', ignore_zero=TRUE, title='Set-NR')

# violinplot_sse_lengths(sNR, ignore_zero=TRUE)
# print_png(fullname('plots/length-sNR.png'), width=4000, ratio=2/1, res=400)
violinplot_sse_lengths(sNR, ignore_zero=TRUE, turn_labels=TRUE)
print_png(fullname('plots/length-sNR-500t.png'), width=4000, ratio=2/1, res=500)
print_tif(fullname('plots/length-sNR-500t.tif'), width=4000, ratio=2/1, res=500)
# violinplot_sse_lengths(sNR, ignore_zero=TRUE, stagger_labels=TRUE)
# print_png(fullname('plots/length-sNR-500s.png'), width=4000, ratio=2/1, res=500)

# violinplot_sse_lengths_multiset2('Bact'=sNR_Bact, 'Euka'=sNR_Euka, ignore_zero=TRUE) 
# print_png(fullname('plots/length-sNR-Bact-Euka.png'), width=4000, ratio=2/1, res=400)
violinplot_sse_lengths_multiset2('Bact'=sNR_Bact, 'Euka'=sNR_Euka, ignore_zero=TRUE, turn_labels=TRUE) 
print_png(fullname('plots/length-sNR-Bact-Euka-500t.png'), width=4000, ratio=2/1, res=500)
print_tif(fullname('plots/length-sNR-Bact-Euka-500t.tif'), width=4000, ratio=2/1, res=500)
# violinplot_sse_lengths_multiset2('Bact'=sNR_Bact, 'Euka'=sNR_Euka, ignore_zero=TRUE, stagger_labels=TRUE) 
# print_png(fullname('plots/length-sNR-Bact-Euka-500s.png'), width=4000, ratio=2/1, res=500)
# point = mean, horizontal line = median

# violinplot_sse_lengths(sCAM, ignore_zero=TRUE, title='Set-CAM')
# print_png(fullname('plots/length-sCAM.png'), width=4000, ratio=2/1, res=400)
# violinplot_sse_lengths(sBM3, ignore_zero=TRUE, title='Set-BM3')
# print_png(fullname('plots/length-sBM3.png'), width=4000, ratio=2/1, res=400)
# violinplot_sse_lengths(s3A4, ignore_zero=TRUE, title='Set-3A4')
# print_png(fullname('plots/length-s3A4.png'), width=4000, ratio=2/1, res=400)

# # Bulges
# View(table_bulge_occurrence(sNR, bNR))
# View(table_bulge_occurrence(sNR %>% filter(Group=='Bact'), bNR %>% filter(Group=='Bact')))
# View(table_bulge_occurrence(sNR %>% filter(Group=='Euka'), bNR %>% filter(Group=='Euka')))
# plot_bulge_occurrence(sNR, bNR, show_confidence = TRUE, alpha=0.05, title='Set-NR', include_parent_strands=FALSE)
# print_png(fullname('plots/bulge_occurrence-sNR.png'), width=4000, ratio=2/1, res=410)
# 
# violinplot_sse(bNR, 'length', ignore_zero = TRUE, title='Set-NR')
# violinplot_sse(sNR, 'longest_G', ignore_zero = TRUE, title='Set-NR')
# print_png('plots/violinplot.png', width=4000, height=2000)
# print_png('plots/violinplot-small.png', width=1800, height=900)
# print_png('plots/violinplot-small.png', width=1800, height=1800)
# 
# table_bulges_with_missing_with_ligands(sALL, bALL %>% filter(label=='4.2' & bulge=='N'), ligands, only_UniProts_with_more_PDBs=TRUE) %>% View()
# table_bulges_with_missing_with_ligands(sALL, bALL %>% filter(label=='3.3' & bulge=='N'), ligands, only_UniProts_with_more_PDBs=TRUE) %>% View()


# # Helix types
# helices = sNR %>% filter(label %in% helix_order & length > 0)
# 
# helices %>% add_predominant_type() %>% group_by(predom_type) %>% count() %>% View()
# pt = table_predominant_type(sNR)  # Requires nested_sses in helices (SecStrAnnotator --verbose)
# pt$predom_type = factor(pt$predom_type)
# pt %>% spread(predom_type, freq) %>% View()
# plot_predominant_type(sNR, title='Predominant helix type (Set-NR)')  # Requires nested_sses in helices (SecStrAnnotator --verbose)
# print_png(fullname('plots/predominant_type-sNR.png'), width=4000, ratio=2/1, res=450)
# 
# table_criteria(helices, have_G=(helices$longest_G>=3), have_H=(helices$longest_H>=4), have_I=(helices$longest_I>=5), by_label=FALSE) %>% spread(set, fulfilling) %>% View()
# table_criteria(helices, have_G=(helices$longest_G>=3), have_H=(helices$longest_H>=4), have_I=(helices$longest_I>=5)) %>% spread(set, fulfilling) %>% View()
# table_criteria(sNR_Bact, have_G=(sNR_Bact$longest_G>=3), have_H=(sNR_Bact$longest_H>=4), have_I=(sNR_Bact$longest_I>=5), ignore_zero=TRUE) %>% spread(set, fulfilling) %>% filter(label %in% helix_order) %>% View()
# table_criteria(sNR_Euka, have_G=(sNR_Euka$longest_G>=3), have_H=(sNR_Euka$longest_H>=4), have_I=(sNR_Euka$longest_I>=5), ignore_zero=TRUE) %>% spread(set, fulfilling) %>% filter(label %in% helix_order) %>% View()
# table_criteria(sNR, have_I=(sNR$longest_I>=5), start_with_I=(sNR$longest_I_start>=5), end_with_I=(sNR$longest_I_end>=5), ignore_zero=TRUE) %>% spread(set, fulfilling) %>% filter(label %in% helix_order) %>% View()
# 
# plot_contained_helix_types(sNR, title='Contained helix types (Set-NR)', y_label='Fraction')
# print_png(fullname('plots/contained_types-sNR.png'), width=4000, ratio=2/1, res=450)
# 
# 
# plot_sse_occurrence(sNR, group_col='predType', show_confidence = TRUE, alpha=0.05, title='Set-NR')
# plot_sse_occurrence(sNR %>% filter(longest_I>=5), show_confidence = TRUE, alpha=0.05, title='Set-NR, longest_I >= 5')
# jitterplot_xy(filter(sNR_Euka, label=="J'"), longest_G, longest_H)
