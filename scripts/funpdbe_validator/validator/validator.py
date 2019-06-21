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

import jsonschema
import json
import re


class Validator(object):
    """
    Validates FunPDBe JSON files

    Example usage:
    validator = Validator()
    validator.load_schema("path/to/schema")
    validator.load_json("path/to/json")
    if validator.basic_checks() and validator.validate_against_schema():
        # the input JSON is a valid FunPDBe file
    """

    def __init__(self, resource):
        self.resource = resource
        self.schema = None
        self.json_data = None
        self.error_log = None

    def load_json(self, path_to_file):
        """
        Loads the user JSON
        :param path_to_file: String, path to JSON file
        :return: None
        """
        self.json_data = self._parse_json(path_to_file)

    def load_schema(self, path_to_schema):
        """
        Loads the FunPDBe schema
        :param path_to_schema: String, path to FunPDBe schema
        :return: None
        """
        self.schema = self._parse_json(path_to_schema)

    def basic_checks(self):
        """
        Performs basic data checks
        :return: Bool, True if valid, False if invalid
        """
        if self._test_resource() and self._test_pdb_id():
            return True
        return False

    def validate_against_schema(self):
        """
        Calls jsonschema.validate() to compare the JSON with
        the FunPDBe JSON schema
        :return: Bool, True is valid, False if invalid
        """
        try:
            jsonschema.validate(self.json_data, self.schema)
            return True
        except jsonschema.exceptions.ValidationError as err:
            self.error_log = "JSON does not comply with schema: %s" % err
            return False

    def _parse_json(self, path):
        """
        Parses a FunPDBe JSON file and in case of file error
        or JSON error, the error message is saved to self.error_log

        :return: Bool, True is parsed, False if failed
        """
        try:
            with open(path, "r") as json_file:
                try:
                    return json.load(json_file)
                except json.decoder.JSONDecodeError as err:
                    self.error_log = "JSON error: %s" % err
                    return None
        except IOError as ioerr:
            self.error_log = "File error: %s" % ioerr
            return None

    def _test_resource(self):
        """
        Check if data_resource field exists in the JSON,
        and if it is the same as the provided resource name
        :return: Bool, True if valid, False if invalid
        """
        if "data_resource" not in self.json_data.keys():
            self.error_log = "No data resource name found"
            return False
        elif self.json_data["data_resource"] != self.resource:
            self.error_log = "Data resource name mismatch"
            return False
        return True

    def _test_pdb_id(self):
        """
        Check if PDB id exists in the JSON, and if it follows
        the PDB id pattern
        :return: Bool, True if valid, False if invalid
        """
        if "pdb_id" not in self.json_data.keys():
            self.error_log = "No PDB id found"
            return False
        elif not re.match("^[1-9][a-zA-Z0-9]{3}$", self.json_data["pdb_id"]):
            self.error_log = "Invalid PDB id found"
            return False
        return True
