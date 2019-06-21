#!/usr/bin/env python3

"""
Copyright 2018 EMBL - European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
either express or implied. See the License for the specific
language governing permissions and limitations under the
License.
"""

import json
from unittest import TestCase

from validator.residue_index import ResidueIndexes

with open("data/test_data.json", "r") as mock_data_file:
    mock_data = json.load(mock_data_file)

with open("data/test_data_multichain.json", "r") as mock_data_file_multichain:
    mock_data_multichain = json.load(mock_data_file_multichain)

mock_data_no_pdb_id = {"foo": "bar"}

mock_data_bad_numbering = {"pdb_id": "2aqa",
                           "chains": [{"chain_label": "A",
                                       "residues": [{"pdb_res_label": "2",
                                                     "aa_type": "ALA"}]}]}


def mock_get_residue_numbering_false(self):
    return False


def mock_get_residue_numbering_true(self):
    return True


def mock_compare_residue_number(self, foo, bar, asd):
    return False


class TestCheckResidueIndices(TestCase):

    def setUp(self):
        self.cri = ResidueIndexes(mock_data)

    def test_loop_chains(self):
        self.cri._get_residue_numbering = mock_get_residue_numbering_false
        result = self.cri.check_every_residue()
        self.assertFalse(result)
        self.cri._get_residue_numbering = mock_get_residue_numbering_true
        result = self.cri.check_every_residue()
        self.assertTrue(result)
        self.cri.pdb_id = None
        self.assertFalse(self.cri.check_every_residue())

    def test_set_pdb_id(self):
        self.assertIsNotNone(self.cri._set_pdb_id())
        bad_cri = ResidueIndexes(mock_data_no_pdb_id)
        self.assertIsNone(bad_cri._set_pdb_id())

    def test_check_numbering(self):
        result = self.cri._check_numbering({}, {}, "A")
        self.assertFalse(result)
        self.cri._compare_residue_number = mock_compare_residue_number
        result = self.cri._check_numbering({}, {"residues": [{"pdb_res_label": 0, "aa_type": "ALA"}]}, "A")
        self.assertFalse(result)

    def test_get_residue_numbering(self):
        mock_data = {"chain_label": "A"}
        self.cri.pdb_id = "1CBS"
        self.cri._check_numbering = lambda x, y, z : True
        result = self.cri._get_residue_numbering(mock_data)
        self.assertTrue(result)
        self.cri.pdb_id = "2H58"
        result = self.cri._get_residue_numbering(mock_data)
        self.assertFalse(result)

    def test_recursive_loop(self):
        result = self.cri._recursive_loop([{"foo": "bar"}], "foo", None, None, "A")
        self.assertFalse(result)

    def test_with_bad_numbering(self):
        cri_with_bad_numbering = ResidueIndexes(mock_data_bad_numbering)
        result = cri_with_bad_numbering.check_every_residue()
        self.assertFalse(result)

    def test_process_residues(self):
        result = self.cri._process_residues(
            [{"author_residue_number": 1, "residue_name": "ALA", "author_insertion_code": ""}], "1", "ALA", "A")
        self.assertTrue(result)
        result = self.cri._process_residues(
            [{"author_residue_number": 1, "residue_name": "ALA", "author_insertion_code": "C"}], "1C", "ALA", "A")
        self.assertTrue(result)
        result = self.cri._process_residues(
            [{"author_residue_number": 2, "residue_name": "ALA", "author_insertion_code": ""}], "1", "ALA", "A")
        self.assertFalse(result)
        result = self.cri._process_residues(
            [{"author_residue_number": 1, "residue_name": "ALA", "author_insertion_code": ""}], "1", "HIS", "A")
        self.assertFalse(result)

    def test_with_multichain(self):
        self.cri = ResidueIndexes(mock_data_multichain)
        result = self.cri.check_every_residue()
        self.assertTrue(result)