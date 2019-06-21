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

import unittest
from validator.validator import Validator


class TestValidator(unittest.TestCase):

    def test_no_file_at_path(self):
        validator = Validator("resource_name")
        validator._parse_json("invalid_path")
        self.assertIn("File error", validator.error_log)

    def test_json_parser(self):
        validator = Validator("resource_name")
        data = validator._parse_json("data/test_data.json")
        self.assertIsNotNone(data)
        validator._parse_json("data/test_data_malformed.json")
        self.assertIn("JSON error", validator.error_log)

    def test_basic_checks(self):
        validator = Validator("test")
        validator.json_data = {"data_resource": "test", "pdb_id": "1abc"}
        self.assertTrue(validator.basic_checks())
        validator.json_data = {"data_resource": "test"}
        self.assertFalse(validator.basic_checks())
        validator.json_data = {"pdb_id": "1abc"}
        self.assertFalse(validator.basic_checks())

    def test_no_resource_name(self):
        validator = Validator("test")
        validator.json_data = {"pdb_id": "1abc"}
        self.assertFalse(validator._test_resource())

    def test_resource_name_mismatch(self):
        validator = Validator("test")
        validator.json_data = {"data_resource": "test2", "pdb_id": "1abc"}
        self.assertFalse(validator._test_resource())

    def test_no_pdb_id(self):
        validator = Validator("test")
        validator.json_data = {"data_resource": "test"}
        self.assertFalse(validator._test_pdb_id())

    def test_invalid_pdb_id(self):
        validator = Validator("test")
        validator.json_data = {"data_resource": "test", "pdb_id": "invalid"}
        self.assertFalse(validator._test_pdb_id())

    def test_json_validation(self):
        validator = Validator("ProKinO")
        validator.load_json("data/test_data.json")
        validator.load_schema("data/funpdbe_schema.json")
        validation = validator.validate_against_schema()
        self.assertTrue(validation)
        validator.load_json("data/test_data_invalid.json")
        validation = validator.validate_against_schema()
        self.assertFalse(validation)