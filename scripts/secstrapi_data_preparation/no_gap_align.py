
import numpy as np
import argparse
import os
from os import path
import shutil
import sys
import json
import heapq
import itertools
from typing import List, Tuple, Dict
from collections import defaultdict

from PIL import Image
from Bio import AlignIO, SeqIO
from Bio.SubsMat import MatrixInfo
import logomaker
from matplotlib import pyplot
# pip3 install Biopython logomaker matplotlib

import lib

################################################################################

GAP_CHAR = '-'
UNKNOWN_CHAR = '?'
GAP_PENALTY = 10
WEBLOGO2 = './weblogo/seqlogo'
WEBLOGO3 = 'weblogo'

def print_aln(seqs, names=None, tree=None, output_file=None):
    columns = [] 
    if tree is not None:
        columns.append(draw_reordered_tree(tree))
    if names is not None:
        max_name_length = max( len(name) for name in names )
        columns.append([ name.ljust(max_name_length) for name in names ])
        columns.append('|' * len(seqs))
    columns.append(seqs)
    columns.append('|' * len(seqs))
    if output_file is None:
        for row in zip(*columns):
            print(*row)
    else:
        with open(output_file, 'w', encoding=lib.DEFAULT_ENCODING) as w:
            for row in zip(*columns):
                w.write(' '.join(row) + '\n')
    #     for seq, name in zip(seqs, names):
    #         print(f'{name.ljust(max_name_length)}', '|', seq, '|')
    # else:
    #     for seq in seqs:
    #         print('|', seq, '|')

def write_fasta(names: List[str], sequences: List[str], filename: str):
    with open(filename, 'w', encoding=lib.DEFAULT_ENCODING) as w:
        for name, seq in zip(names, sequences):
            w.write(f'>{name}\n{seq}\n')

def substitution_matrix(scores: Dict[Tuple[str, str], float], alphabet=None, gap_penalty=0):
    if alphabet is None:
        alphabet = [GAP_CHAR] + sorted(set( letter for substitution, score in scores for letter in substitution ) - {GAP_CHAR, UNKNOWN_CHAR}) + [UNKNOWN_CHAR]
    letter2index = { letter: i for i, letter in enumerate(alphabet) }
    n = len(alphabet)
    matrix = np.zeros((n, n))
    for (x, y), score in scores.items():
        xi = letter2index.get(x, None)
        yi = letter2index.get(y, None)
        if xi is not None and yi is not None:
            matrix[xi, yi] = score
            matrix[yi, xi] = score
    matrix[1:, 1:] += gap_penalty
    return matrix, alphabet, letter2index

def letter2vector(letter, letter2index):
    vec = np.zeros((len(letter2index),))
    vec[letter2index[letter]] = 1
    return vec

def sequence2matrix(sequence, letter2index):
    n = len(sequence)
    m = len(letter2index)
    matrix = np.zeros((n, m))
    for i, letter in enumerate(sequence):
        matrix[i, letter2index[letter]] = 1
    return matrix

def calculate_score(letter1, letter2, subst_matrix):
    return np.matmul(np.matmul(letter1, subst_matrix), letter2)

def calculate_score_mat(seq_mat_1, seq_mat_2, subst_matrix):
    scores = np.matmul(seq_mat_1, subst_matrix) * seq_mat_2
    return scores.sum()

def optimal_shift_and_score(seq_mat_1, seq_mat_2, subst_matrix, return_bestness=False):
    n1 = seq_mat_1.shape[0]
    n2 = seq_mat_2.shape[0]
    seq_mat_1_M = np.matmul(seq_mat_1, subst_matrix)
    seq_mat_1_M_2 = np.matmul(seq_mat_1_M, seq_mat_2.transpose())
    MIN_SHIFT = -n2+1  # shift of sequence 2 rightwards wrt. sequence 1
    MAX_SHIFT = n1-1
    # shifts_scores = []
    # for shift in range(MIN_SHIFT, MAX_SHIFT+1):
        # start = max(-shift, 0)  # wrt. sequence 2
        # stop = min(n1-shift, n2)  # wrt. sequence 2
        # score = np.sum(seq_mat_1_M[start+shift:stop+shift, :] * seq_mat_2[start:stop, :])
        # shifts_scores.append((shift, score))
    shifts_scores = [ (shift, seq_mat_1_M_2.diagonal(offset=-shift).sum()) for shift in range(MIN_SHIFT, MAX_SHIFT+1) ]
    (best_shift, best_score), second_best = two_max(shifts_scores, key = lambda t: t[1])
    # print(best_score, second_best_score, bestness)
    if return_bestness:
        if second_best is not None:
            (second_best_shift, second_best_score) = second_best
            bestness = (best_score - second_best_score) / best_score
            # print(bestness, best_shift, best_score, second_best_shift, second_best_score)
        else:
            bestness = 1.0
        return best_shift, best_score, bestness
    else:
        return best_shift, best_score

def two_max(iterable, key=lambda x: x):
    first_max = None
    second_max = None
    for x in iterable:
        if first_max is None or key(x) > key(first_max):
            first_max, second_max = x, first_max
        elif second_max is None or key(x) > key(second_max):
            second_max = x
    return first_max, second_max

def combine_sequences(seq_mat_1, seq_mat_2, weights=(1, 1), shift=None, subst_matrix=None):
    if shift is None:
        if subst_matrix is None:
            raise Exception('Either shift or substitution matrix must be provided')
        else:
            shift, score = optimal_shift_and_score(seq_mat_1, seq_mat_2, subst_matrix)
    if shift < 0:
        return combine_sequences(seq_mat_2, seq_mat_1, weights=tuple(reversed(weights)), shift=-shift)
    n1, m = seq_mat_1.shape
    n2 = seq_mat_2.shape[0]
    n = max(n1, n2+shift)
    w1 = weights[0] / sum(weights)
    w2 = weights[1] / sum(weights)
    result = np.zeros((n, m))
    result[0:n1, :] += w1 * seq_mat_1
    result[shift:n2+shift, :] += w2 * seq_mat_2
    result[:, 0] = 1.0 - result[:, 1:].sum(axis=1)  # calculate probabilities for GAP_CHAR
    return result

def multialign(sequence_matrices, subst_matrix):
    sequence_matrices = sequence_matrices[:]
    n = len(sequence_matrices)
    N = 2*n - 1
    scores = np.zeros((N, N))
    active_nodes = set(range(n))
    weights = [1] * n
    tree = np.full((N, 3), -1, dtype=int)  # (left child, right child, shift of right wrt left) for each node, -1 = leaf/uninitialized
    i_j_shift_scores = [ (i, j, *optimal_shift_and_score(sequence_matrices[i], sequence_matrices[j], subst_matrix)) 
        for (i, j) in itertools.combinations(active_nodes, 2) ]
    queue = PriorityQueue( ((i, j), (-score, shift)) for i, j, shift, score in i_j_shift_scores )
    while len(active_nodes) > 1:
        best_pair = queue.pop_min_which(lambda ij: ij[0] in active_nodes and ij[1] in active_nodes)
        if best_pair is None:
            break  # no joinable pairs, algorithm has converged
        (i, j), (neg_score, shift) = best_pair
        new = len(sequence_matrices)
        new_matrix = combine_sequences(sequence_matrices[i], sequence_matrices[j], weights=(weights[i], weights[j]), shift=shift)
        new_weight = weights[i] + weights[j]
        sequence_matrices.append(new_matrix)
        weights.append(new_weight)
        active_nodes.remove(i)
        active_nodes.remove(j)
        tree[new, :] = (i, j, shift)
        for node in active_nodes:
            shift, score = optimal_shift_and_score(sequence_matrices[node], sequence_matrices[new], subst_matrix)
            queue.add((node, new), (-score, shift))
        active_nodes.add(new)
    alignment_matrix = sequence_matrices[-1]
    shifts = shifts_from_tree(tree)
    return alignment_matrix, shifts, tree

def multirealign(reference_sequence_matrix, sequence_matrices, subst_matrix):
    shifts, scores, bestnesses = zip(*( optimal_shift_and_score(reference_sequence_matrix, mat, subst_matrix, return_bestness=True) for mat in sequence_matrices ))
    min_shift = min(shifts)
    shifts = [ shift - min_shift for shift in shifts ]
    max_length = max( shift + seq_mat.shape[0] for shift, seq_mat in zip(shifts, sequence_matrices) )
    alignment_matrix = np.zeros((max_length, reference_sequence_matrix.shape[1]))
    for shift, seq_mat in zip(shifts, sequence_matrices):
        # print(shift)
        # print(seq_mat.shape)
        # print(alignment_matrix[shift:seq_mat.shape[0]+shift, :].shape)
        alignment_matrix[shift:seq_mat.shape[0]+shift, :] += seq_mat
    alignment_matrix /= len(sequence_matrices)
    # print(*( f'{be:.4f} {sc:.2f} {sh},' for be, sc, sh in sorted(zip(bestnesses, scores, shifts)) ))
    # print(sorted(scores))
    alignment_matrix[:, 0] = 1.0 - alignment_matrix[:, 1:].sum(axis=1)  # calculate probabilities for GAP_CHAR
    return alignment_matrix, shifts, bestnesses

def shifts_from_tree_aux(tree, root, current_shift, result_array):
    left, right, shift = tree[root, :]
    if left < 0:  # leaf
        result_array[root] = current_shift
    else:  # internal node
        if shift >= 0:
            shifts_from_tree_aux(tree, left, current_shift, result_array)
            shifts_from_tree_aux(tree, right, current_shift + shift, result_array)
        else:
            shifts_from_tree_aux(tree, left, current_shift - shift, result_array)
            shifts_from_tree_aux(tree, right, current_shift, result_array)

def shifts_from_tree(tree):
    N, m = tree.shape
    if m != 3:
        raise
    n_leaves = (N+1) // 2
    root = N-1
    current_shift = 0
    result_array = np.zeros(N, dtype=int)
    shifts_from_tree_aux(tree, root, current_shift, result_array)
    return result_array[:n_leaves]

def reordering_from_tree(tree, root=None):
    if root is None:
        root = tree.shape[0] - 1
    left, right, shift = tree[root, :]
    if left < 0:  # leaf
        yield root
    else:  # internal node
        if shift >= 0:
            yield from reordering_from_tree(tree, left)
            yield from reordering_from_tree(tree, right)
        else:
            yield from reordering_from_tree(tree, right)
            yield from reordering_from_tree(tree, left)

def apply_shifts(sequences: List[str], shifts: List[int]) -> List[str]:
    max_length = max( len(seq) + shift for seq, shift in zip(sequences, shifts) )
    return [ GAP_CHAR*shift + seq + GAP_CHAR*(max_length-len(seq)-shift) for seq, shift in zip(sequences, shifts) ]

def draw_reordered_tree(tree):
    root = tree.shape[0] - 1
    lines, start = draw_reordered_tree_aux(tree, root)
    return lines

def draw_reordered_tree_aux(tree, root, put_root_down=False):
    left, right, shift = tree[root, :]
    if left < 0:  # leaf
        return ['─'], 0
    else:  # internal node
        top, bottom = (left, right) if shift >= 0 else (right, left)
        fig1, start1 = draw_reordered_tree_aux(tree, top, put_root_down=True)
        fig2, start2 = draw_reordered_tree_aux(tree, bottom)
        height1 = len(fig1)
        height2 = len(fig2)
        height = height1 + height2
        width1 = len(fig1[0])
        width2 = len(fig2[0])
        width = max(width1, width2)
        for i in range(height1):
            if i == start1:
                fig1[i] = fig1[i].rjust(width, '─')
            else:
                fig1[i] = fig1[i].rjust(width)
        for i in range(height2):
            if i == start2:
                fig2[i] = fig2[i].rjust(width, '─')
            else:
                fig2[i] = fig2[i].rjust(width)
        fig = fig1 + fig2
        start2 += height1
        # start = (start1 + start2 + 1) // 2 if put_root_down else (start1 + start2) // 2
        if start2 - start1 < 2:
            start = start2 if put_root_down else start1
        else:
            start = start2-1 if put_root_down else start1+1
        for i in range(height):
            if i == start and i == start1:
                fig[i] = '┬' + fig[i]
            elif i == start and i == start2:
                fig[i] = '┴' + fig[i]
            elif i == start:
                fig[i] = '┤' + fig[i]
            elif i == start1:
                fig[i] = '┌' + fig[i]
            elif start1 < i < start2:
                fig[i] = '│' + fig[i]
            elif i == start2:
                fig[i] = '└' + fig[i]
            else:
                fig[i] = ' ' + fig[i]
        return fig, start

def logo_heights_widths_areas(sequence_matrix, first_column_is_gap=True):
    if first_column_is_gap:
        gap_prob = sequence_matrix[:, 0]
        probs = sequence_matrix[:, 1:].copy()
    else:
        gap_prob = 1 - sequence_matrix.sum(axis=1)
        probs = sequence_matrix.copy()
    probs /= probs.sum(axis=1, keepdims=True)  # normalize probabilities
    probs[probs<=0] = 1  # so that np.log2 throws no error, does not affect the result
    background_entropy = -np.log2(1/20)
    entropies = -np.sum(probs * np.log2(probs), axis=1)
    heights = background_entropy - entropies
    widths = 1 - gap_prob
    areas = heights * widths
    return heights, widths, areas

def logo_maxprobs(sequence_matrix, first_column_is_gap=True):
    if first_column_is_gap:
        probs = sequence_matrix[:, 1:]
    else:
        probs = sequence_matrix
    maxprobs = probs.max(axis=1)
    return maxprobs

def get_widest_and_highest_column_index(sequence_matrix):
    heights, widths, areas = logo_heights_widths_areas(sequence_matrix)
    max_width, max_height, index = max( (width, height, i) for i, (width, height) in enumerate(zip(widths, heights)) )
    return index

def get_highest_column_index(sequence_matrix):
    heights, widths, areas = logo_heights_widths_areas(sequence_matrix)
    max_height, max_height_index = max( (height, i) for i, height in enumerate(heights) )
    return max_height_index

def get_max_area_column_index(sequence_matrix):
    heights, widths, areas = logo_heights_widths_areas(sequence_matrix)
    # max_area, max_area_index = max( (area, i) for i, area in enumerate(areas) )  # Takes the last in case of ties
    max_area = max(areas)
    max_area_index = next(i for i, area in enumerate(areas) if area == max_area)  # Takes the first in case of ties
    # print(f'Max. area column: {max_area_index}, area: {max_area}', file=sys.stderr)
    return max_area_index
    
def get_max_maxprob_column_index(sequence_matrix):
    maxprobs = logo_maxprobs(sequence_matrix)
    max_mp, max_mp_index = max( (mp, i) for i, mp in enumerate(maxprobs) )
    return max_mp_index

def get_pivot_column_index(sequence_matrix, probability=False):
    if probability:
        return get_max_maxprob_column_index(sequence_matrix)
    else:
        return get_max_area_column_index(sequence_matrix)


def run_weblogo2(alignment_file, logo_file, first_index=0):
    with open(alignment_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        for line in r:
            if line[0] != '>':
                n_residues = len(line.strip())
                break
    height = 8
    width_per_residue = 0.8
    title = path.split(logo_file)[1]
    title = path.splitext(title)[0]
    title = 'Helix ' + title if title[0].isalpha() else 'Strand ' + title
    command = f'{WEBLOGO2} -caMnY -F PNG -h {height} -w {width_per_residue * n_residues} -s {first_index} -t "{title}" -f "{alignment_file}" > "{logo_file}"'
    # print(command)
    os.system(command)

def run_weblogo3(alignment_file, logo_file, first_index=0):
    # generate sequence logos using WebLogo
    # WebLogo documentation: http://weblogo.threeplusone.com/manual.html#CLI
    title = path.split(logo_file)[1]
    title = path.splitext(title)[0]
    title = 'Helix ' + title if title[0].isalpha() else 'Strand ' + title
    composition = 'equiprobable' # 'equiprobable' 'none' 'auto'
    units = 'bits' # 'bits' 'probability'
    command = f'''{WEBLOGO3}  --format png  --resolution 600  --stacks-per-line 60  --fineprint ""  --errorbars NO
        --rotate-numbers YES  --number-interval 1  --aspect-ratio 6  --logo-font ArialBold  --title-font TimesNewRomanBold  --scale-width YES
        --sequence-type protein  --units {units}  --composition {composition}
        --first-index {first_index}  --title "{title}"  --fin "{alignment_file}"  --fout "{logo_file}"
        '''.replace('\n', ' ')
    # print(command)
    os.system(command)

def run_logomaker_from_matrix(matrix, widths, logo_file, first_index=0, dpi=600, units='bits', color_scheme='weblogo_protein', title=None):
    '''Generate sequence logo using Logomaker.
    matrix: pandas dataframe (columns correspond to amino acid letters)
    widths: array of widths of individual columns
    units: 'bits' or 'probability', color_scheme: 'weblogo_protein', 'hydrophobicity' ...
    # Logomaker documentation: https://logomaker.readthedocs.io/en/latest/'''
    scale = 1.0
    height = 3.0
    width_per_residue = 0.3  # 0.35
    width_extra = 1.0  # 0.5
    width_minimum = 2.0
    vpad = 0.05
    hpad = 0.05
    font_name = 'DejaVu Sans Mono'
    y_tick_spacing = 0.5
    y_label_spacing = 1
    y_label_format = '{:.0f}'

    if title is None:
        title = path.split(logo_file)[1]
        title = path.splitext(title)[0]
        title = 'Helix ' + title if title[0].isalpha() else 'Strand ' + title

    matrix.index += first_index
    n_residues = matrix.shape[0]

    width = max(n_residues*width_per_residue + width_extra, width_minimum)
    figsize = (width * scale, height * scale)
    logo = logomaker.Logo(matrix, figsize=figsize, vpad=vpad, color_scheme=color_scheme, font_name=font_name)
    logo.fig.set_dpi(100)
    for glyph in logo.glyph_list:
        p = int(glyph.p)  # position
        c = glyph.c  # character
        logo.style_single_glyph(p, c, width=(1 - hpad) * widths[p-first_index])

    max_y = np.log2(20) if units == 'bits' else 1.0

    logo.ax.yaxis.set_view_interval(0, max_y)
    logo.ax.title.set_text(title)
    logo.ax.title.set_fontsize(20*scale)
    logo.ax.title.set_fontweight('bold')
    # logo.ax.set_xlabel('position', labelpad=0, fontsize=14*scale)
    y_label = 'Information content [bits]' if units=='bits' else units
    logo.ax.set_ylabel(y_label, labelpad=3.0, fontsize=12*scale)
    logo.style_xticks(spacing=1, anchor=0, rotation=90, fmt='%d', fontsize=14*scale)
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom', 'top', 'right'], visible=True, linewidth=1)

    logo.ax.set_yticks(np.arange(0, max_y, y_tick_spacing))
    logo.fig.tight_layout()
    ytl = logo.ax.get_yticklabels()
    for lab in ytl:
        y = lab.get_position()[1]
        if y % y_label_spacing == 0:
            lab.set_text(y_label_format.format(y))
        else:
            lab.set_text('')
    logo.ax.set_yticklabels(ytl, fontsize=14*scale)
    logo.fig.tight_layout()
    logo.fig.savefig(logo_file, dpi=dpi)
    pyplot.close(logo.fig)
    if str(logo_file).endswith(('.tif', '.tiff')):
        compress_tiff(logo_file)

def run_logomaker(alignment_file, logo_file, first_index=0, dpi=600, units='bits', color_scheme='weblogo_protein', title=None):
    '''Generate sequence logo using Logomaker.
    units: 'bits' or 'probability', color_scheme: 'weblogo_protein', 'hydrophobicity' ...
    # Logomaker documentation: https://logomaker.readthedocs.io/en/latest/'''
    inp = SeqIO.parse(alignment_file, 'fasta')
    names, seqs = zip(*( (x.id, str(x.seq)) for x in inp ))
    counts_mat = logomaker.alignment_to_matrix(seqs, characters_to_ignore=GAP_CHAR+UNKNOWN_CHAR)
    prob_mat = counts_mat / len(seqs)
    heights, widths, areas = logo_heights_widths_areas(prob_mat.values, first_column_is_gap=False)
    if units == 'bits':
        matrix = heights.reshape((-1, 1)) * prob_mat / prob_mat.values.sum(axis=1, keepdims=True)
    elif units == 'probability':
        matrix = prob_mat / prob_mat.values.sum(axis=1, keepdims=True)
    else:
        raise ValueError("units must be 'bits' or 'probability'")
    run_logomaker_from_matrix(matrix, widths, logo_file, first_index=first_index, dpi=dpi, units=units, color_scheme=color_scheme, title=title)
    
def run_logomaker_backup(alignment_file, logo_file, first_index=0, dpi=600, units='bits', color_scheme='weblogo_protein'):
    '''Generate sequence logo using Logomaker.
    units: 'bits' or 'probability', color_scheme: 'weblogo_protein', 'hydrophobicity' ...
    # Logomaker documentation: https://logomaker.readthedocs.io/en/latest/'''
    scale = 1.0
    height = 3.0
    width_per_residue = 0.3  # 0.35
    width_extra = 1.0  # 0.5
    width_minimum = 2.0
    vpad = 0.05
    hpad = 0.05
    font_name = 'DejaVu Sans Mono'
    y_tick_spacing = 0.5
    y_label_spacing = 1
    y_label_format = '{:.0f}'

    title = path.split(logo_file)[1]
    title = path.splitext(title)[0]
    title = 'Helix ' + title if title[0].isalpha() else 'Strand ' + title

    inp = SeqIO.parse(alignment_file, 'fasta')
    names, seqs = zip(*( (x.id, str(x.seq)) for x in inp ))
    counts_mat = logomaker.alignment_to_matrix(seqs, characters_to_ignore=GAP_CHAR+UNKNOWN_CHAR)
    prob_mat = counts_mat / len(seqs)
    
    heights, widths, areas = logo_heights_widths_areas(prob_mat.values, first_column_is_gap=False)
    if units == 'bits':
        matrix = heights.reshape((-1, 1)) * prob_mat / prob_mat.values.sum(axis=1, keepdims=True)
    elif units == 'probability':
        matrix = prob_mat / prob_mat.values.sum(axis=1, keepdims=True)
    else:
        raise ValueError("units must be 'bits' or 'probability'")

    matrix.index += first_index
    n_residues = matrix.shape[0]

    width = max(n_residues*width_per_residue + width_extra, width_minimum)
    figsize = (width * scale, height * scale)
    logo = logomaker.Logo(matrix, figsize=figsize, vpad=vpad, color_scheme=color_scheme, font_name=font_name)
    logo.fig.set_dpi(100)
    for glyph in logo.glyph_list:
        p = int(glyph.p)  # position
        c = glyph.c  # character
        logo.style_single_glyph(p, c, width=(1 - hpad) * widths[p-first_index])

    max_y = np.log2(20) if units == 'bits' else 1.0

    logo.ax.yaxis.set_view_interval(0, max_y)
    logo.ax.title.set_text(title)
    logo.ax.title.set_fontsize(20*scale)
    logo.ax.title.set_fontweight('bold')
    # logo.ax.set_xlabel('position', labelpad=0, fontsize=14*scale)
    y_label = 'Information content [bits]' if units=='bits' else units
    logo.ax.set_ylabel(y_label, labelpad=3.0, fontsize=12*scale)
    logo.style_xticks(spacing=1, anchor=0, rotation=90, fmt='%d', fontsize=14*scale)
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom', 'top', 'right'], visible=True, linewidth=1)

    logo.ax.set_yticks(np.arange(0, max_y, y_tick_spacing))
    logo.fig.tight_layout()
    ytl = logo.ax.get_yticklabels()
    for lab in ytl:
        y = lab.get_position()[1]
        if y % y_label_spacing == 0:
            lab.set_text(y_label_format.format(y))
        else:
            lab.set_text('')
    logo.ax.set_yticklabels(ytl, fontsize=14*scale)
    logo.fig.tight_layout()
    logo.fig.savefig(logo_file, dpi=dpi)
    pyplot.close(logo.fig)
    if logo_file.endswith(('.tif', '.tiff')):
        compress_tiff(logo_file)

def compress_tiff(filename: str):
    name, ext = path.splitext(filename)
    temp_name = f'{name}-compressed{ext}'
    im = Image.open(filename)
    im.save(temp_name, compression='tiff_lzw')
    im.close()  # `with Image.open(filename)...` fails to close the file on Windows!
    os.remove(filename)
    shutil.move(temp_name, filename)


class PriorityQueue:
	def __init__(self, elements_keys):
		self.heap = list(( (key, i, elem) for (i, (elem, key)) in enumerate(elements_keys) ))
		heapq.heapify(self.heap)
		self.seq = len(self.heap)
	def is_empty(self):
		return len(self.heap) == 0
	def add(self, elem, key):
		heapq.heappush(self.heap, (key, self.seq, elem))
		self.seq += 1
	def pop_min(self):
		if not self.is_empty():
			key, i, elem = heapq.heappop(self.heap)
			return elem, key
		else:
			return None
	def pop_min_which(self, predicate):
		while not self.is_empty():
			key, i, elem = heapq.heappop(self.heap)
			if predicate(elem):
				return elem, key
		return None


################################################################################

class NoGapAligner:
    def __init__(self, subst_matrix_info=MatrixInfo.blosum62, gap_penalty=10, realign=True):
        self.subst_matrix, self.alphabet, self.letter2index = substitution_matrix(subst_matrix_info, gap_penalty=gap_penalty)
        self.realign = realign

    def align(self, sequences, names=None):
        if isinstance(sequences, str):
            inp = SeqIO.parse(sequences, 'fasta')
            self.names, self.seqs = zip(*( (x.id, str(x.seq)) for x in inp ))
        elif names is None:
            self.seqs = list(sequences)
            self.names = [ str(i) for i in range(len(self.seqs)) ]
        else:
            self.seqs = list(sequences)
            self.names = list(names)
            if len(self.seqs) != len(self.names):
                raise Exception('There must be the same number of sequences and names')
        sequence_matrices = [ sequence2matrix(seq, self.letter2index) for seq in self.seqs ]
        self.alignment_matrix, self.shifts, self.tree = multialign(sequence_matrices, self.subst_matrix)
        if self.realign:
            last_shifts = None
            while last_shifts is None or any( last != curr for last, curr in zip(last_shifts, self.shifts) ):
                last_shifts = self.shifts
                self.alignment_matrix, self.shifts, bestnesses = multirealign(self.alignment_matrix, sequence_matrices, self.subst_matrix)
                # print(f'Bestness: min {min(bestnesses):.4f}, max {max(bestnesses):.4f}, mean {np.mean(bestnesses):.4f}, median {np.median(bestnesses):.4f}')
        self.aln_seqs = apply_shifts(self.seqs, self.shifts)

    def output_alignment(self, output_file, keep_order=False):
        if not keep_order:
            aln_seqs, names = zip(*( (self.aln_seqs[i], self.names[i]) for i in reordering_from_tree(self.tree) ))
        write_fasta(names, aln_seqs, output_file)

    def print_tree(self, output_file=None):
        aln_seqs, names = zip(*( (self.aln_seqs[i], self.names[i])  for i in reordering_from_tree(self.tree) ))
        print_aln(aln_seqs, names=names, tree=self.tree, output_file=output_file)

    def output_logo(self, output_file, tool='weblogo3', pivot_as=0, units='bits'):
        ''' tool in ['weblogo2', 'weblogo3', 'logomaker'] '''
        height = 8
        if pivot_as is not None:
            first_index = pivot_as - get_pivot_column_index(self.alignment_matrix, probability=(units=='probability'))
        else:
            first_index = 1
        alignment_file = output_file + '.fasta.tmp'
        self.output_alignment(alignment_file)
        if tool == 'weblogo2':
            run_weblogo2(alignment_file, output_file, first_index=first_index)
        elif tool == 'weblogo3':
            run_weblogo3(alignment_file, output_file, first_index=first_index)
        elif tool == 'logomaker':
            # run_logomaker(alignment_file, output_file, first_index=first_index, color_scheme='weblogo_protein')
            run_logomaker(alignment_file, output_file, first_index=first_index, color_scheme='weblogo_protein', units=units)
        else:
            raise NotImplementedError(f'Unknown tool: {tool}')
        os.remove(alignment_file)
    
    def print_column_statistics(self, label='', include_header=False, only_header=False, file=None):
        file = file or sys.stdout
        sep = '\t'
        rows = []
        if include_header or only_header:
            rows.append(('Label', 'Weighted-average height', 'Max. area', 'Max. area index', 'Max. height', 'Max. height index'))
        if not only_header:
            heights, widths, areas = logo_heights_widths_areas(self.alignment_matrix)
            max_area_index = get_max_area_column_index(self.alignment_matrix)
            max_area = areas[max_area_index]
            max_height_index = get_highest_column_index(self.alignment_matrix)
            max_height = heights[max_height_index]
            avg_height = sum(areas) / sum(widths)
            rows.append((label, avg_height, max_area, max_area_index, max_height, max_height_index))
        for row in rows:
            print(*row, sep=sep, file=file)


class Realigner:
    def __init__(self, reference_alignment_file, subst_matrix_info=MatrixInfo.blosum62, gap_penalty=10):
        self.subst_matrix, self.alphabet, self.letter2index = substitution_matrix(subst_matrix_info, gap_penalty=gap_penalty)
        inp = SeqIO.parse(reference_alignment_file, 'fasta')
        names, seqs = zip(*( (x.id, str(x.seq)) for x in inp ))
        sequence_matrices = [ sequence2matrix(seq, self.letter2index) for seq in seqs ]
        self.reference_alignment_matrix = sum(sequence_matrices) / len(sequence_matrices)
        self.pivot_index = get_pivot_column_index(self.reference_alignment_matrix)
        
    def aligning_shift_and_pivot(self, sequence):
        seq_matrix = sequence2matrix(sequence, self.letter2index)
        shift, score = optimal_shift_and_score(self.reference_alignment_matrix, seq_matrix, self.subst_matrix)
        pivot = self.pivot_index - shift
        return shift, pivot

# def main(input_file, output_file, logo_output=None, keep_order=False, print_tree=False, realign=True):
#     inp = SeqIO.parse(input_file, 'fasta')
#     names, seqs = zip(*( (x.id, str(x.seq)) for x in inp ))

#     subst_matrix, alphabet, letter2index = substitution_matrix(MatrixInfo.blosum62, gap_penalty=GAP_PENALTY)

#     # print('Matrix:', subst_matrix)
#     # print('Alphabet:', alphabet)
#     # print('Letter2index:', letter2index)

#     sequence_matrices = [ sequence2matrix(seq, letter2index) for seq in seqs ]
#     matAln, shifts, tree = multialign(sequence_matrices, subst_matrix)

#     if realign:
#         last_shifts = None
#         while last_shifts is None or any( last != curr for last, curr in zip(last_shifts, shifts) ):
#             last_shifts = shifts
#             matAln, shifts, bestnesses = multirealign(matAln, sequence_matrices, subst_matrix)
#             # print(f'Bestness: min {min(bestnesses):.4f}, max {max(bestnesses):.4f}, mean {np.mean(bestnesses):.4f}, median{np.median(bestnesses):.4f}')

#     aln_seqs = apply_shifts(seqs, shifts)

#     # print('Aln\n', matAln)
#     # print('Shifts\n', shifts)
#     # print('Tree\n', tree)

#     # print_aln(seqs, names=names)
#     # print('#')
#     if not keep_order:
#         aln_seqs, names = zip(*( (aln_seqs[i], names[i])  for i in reordering_from_tree(tree) ))
#     write_fasta(names, aln_seqs, output_file)
#     if print_tree:
#         print_aln(aln_seqs, names=None, tree=tree)

#     if logo_output is not None:
#         height = 8
#         width_per_residue = 0.8
#         n_residues = matAln.shape[0]
#         max_height_index = get_pivot_column_index(matAln)
#         # run_weblogo2(output_file, logo_output, first_index=-max_height_index)
#         run_weblogo3(output_file, logo_output, first_index=-max_height_index)
            
# ################################################################################

# if __name__ == "__main__":
#     # parser = argparse.ArgumentParser()
#     # parser.add_argument('input', help='Input file with sequences in multi-FASTA format', type=str)
#     # parser.add_argument('output', help='Output file for aligned sequences in multi-FASTA format', type=str)
#     # parser.add_argument('--logo_output', help='Output file for aligned sequences in multi-FASTA format', type=str, default=None)
#     # parser.add_argument('--keep_order', help='Print aligned sequence in the same order as they appeared in the input (otherwise will order them by the tree)', action='store_true')
#     # args = parser.parse_args()

#     # main(input_file = args.input, output_file = args.output, logo_output = args.logo_output, keep_order = args.keep_order)

#     directory = "/home/adam/Workspace/Python/Ubertemplate/alignment_data/"
#     labels = ["1a", "1b", "1c", "1d", "1e", "2a", "2b", "3a", "3b", "3c", "4a", "4b", "4c", "A", "A'", "B", "B'", "B''", "B'''", "C", "D", "E", "F", "F'", "G", "G'", "H", "I", "J", "J'", "K", "K'", "K''", "K'''", "L", "L'"]
#     # labels = ["F"]
#     os.makedirs(path.join(directory, 'aligned'), exist_ok=True)
#     os.makedirs(path.join(directory, 'logos'), exist_ok=True)
#     for label in labels:
#         print(label)
#         main(
#             path.join(directory, 'sequences', f'extracted_sequences_{label}.fasta'),
#             path.join(directory, 'aligned', f'{label}.fasta'),
#             logo_output=path.join(directory, 'logos', f'{label}.png'),
#             keep_order=False, 
#             print_tree=True)
