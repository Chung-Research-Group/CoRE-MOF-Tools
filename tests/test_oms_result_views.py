"""Verify OMS table views without rerunning chemistry or loading CIFs."""
import ast
import copy
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock


SOURCE = Path(__file__).resolve().parents[1] / 'CoREMOF/calculation/mof_collection.py'
tree = ast.parse(SOURCE.read_text())
original_class = next(item for item in tree.body if isinstance(item, ast.ClassDef) and item.name == 'MofCollection')
properties = [item for item in original_class.body if getattr(item, 'name', None) in {'mof_oms_df', 'metal_site_df'}]
view_class = ast.ClassDef(name='ResultViews', bases=[], keywords=[], body=properties, decorator_list=[])
code = compile(ast.fix_missing_locations(ast.Module(body=[view_class], type_ignores=[])), str(SOURCE), 'exec')


class Frame:
    def __init__(self, data):
        self.data = copy.deepcopy(data)

    @classmethod
    def from_dict(cls, data, orient):
        assert orient == 'index'
        return cls(data)


class OmsResultViewTests(unittest.TestCase):
    def make_views(self, dataframe=Frame):
        namespace = {'pd': SimpleNamespace(DataFrame=dataframe)}
        exec(code, namespace)
        result = namespace['ResultViews']()
        result._metal_site_df = result._mof_oms_df = None
        result._validate_properties = Mock(return_value=([], True))
        result.mof_coll = [{'checksum': 'fixture'}]
        result.properties = {'fixture': {
            'name': 'MOF_A', 'metal_species': ['Cu'], 'has_oms': True,
            'metal_sites': [{'metal': 'Cu', 'is_open': True, 'unique': True,
                            'all_dihedrals': [0.0, 15.0], 'min_dihedral': 0.0}],
        }}
        return result

    def test_mof_and_site_views_have_separate_caches_in_either_order(self):
        for order in (('mof_oms_df', 'metal_site_df'), ('metal_site_df', 'mof_oms_df')):
            with self.subTest(order=order):
                result = self.make_views()
                first, second = (getattr(result, name) for name in order)
                self.assertIsNot(first, second)
                self.assertIs(getattr(result, order[0]), first)
                self.assertIs(getattr(result, order[1]), second)
                self.assertEqual(result.mof_oms_df.data, {
                    'MOF_A': {'Metal Types': 'Cu', 'Has OMS': 'Yes', 'OMS Types': 'Cu'}})
                self.assertEqual(result.metal_site_df.data, {
                    'MOF_A_0': {'metal': 'Cu', 'is_open': True, 'unique': True, 'mof_name': 'MOF_A'}})

    def test_generating_views_does_not_mutate_recorded_site_evidence(self):
        result = self.make_views()
        before = copy.deepcopy(result.properties)
        result.metal_site_df
        result.mof_oms_df
        self.assertEqual(result.properties, before)

    def test_incomplete_results_never_become_a_cached_complete_view(self):
        result = self.make_views()
        result._validate_properties.return_value = (['fixture'], False)
        self.assertIs(result.mof_oms_df, False)
        self.assertIs(result.metal_site_df, False)
        self.assertIsNone(result._mof_oms_df)
        self.assertIsNone(result._metal_site_df)

    def test_real_pandas_view_columns_remain_distinct(self):
        try:
            from pandas import DataFrame
        except ImportError:
            self.skipTest('pandas is optional in dependency-minimal tests')
        result = self.make_views(DataFrame)
        mof = result.mof_oms_df
        sites = result.metal_site_df
        self.assertEqual(list(mof.columns), ['Metal Types', 'Has OMS', 'OMS Types'])
        self.assertEqual(list(sites.index), ['MOF_A_0'])
        self.assertEqual(sites.loc['MOF_A_0', 'mof_name'], 'MOF_A')
        self.assertIs(result.mof_oms_df, mof)
        self.assertIs(result.metal_site_df, sites)


if __name__ == '__main__':
    unittest.main()
