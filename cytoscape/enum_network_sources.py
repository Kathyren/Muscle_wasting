
from enum import Enum
import json

import yaml


class NetworkSource(Enum):
    GENE_MANIA = 'gene_mania'
    STRING = 'stringdb'
    INTACT = 'intact'

    def __init__(self, source_name):
        self.source_name = source_name
        self._node_keys = None
        self._edge_keys = None
        self._main_name = None

    def get_node_keys(self):
        if self._node_keys is None:
            self._load_keys()
        return self._node_keys

    def get_edge_keys(self):
        if self._edge_keys is None:
            self._load_keys()
        return self._edge_keys

    def get_main_name(self):
        if self._main_name is None:
            self._load_keys()
        return self._main_name

    def _load_keys(self):
        with open(f'Data/Apps_details.yaml', 'r') as f:
            keys = yaml.safe_load(f)
            keys = keys[self.source_name]
        self._node_keys = keys.get('node_keys', {})
        self._edge_keys = keys.get('edge_keys', {})
        self._main_name = keys.get('main_name', {})

