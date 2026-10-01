"""Crossref date semantics without HTTP requests or mandatory dependencies."""
import contextlib
import io
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from CoREMOF.calculation import get_info as api


class DateTests(unittest.TestCase):
    def test_online_has_priority_over_print(self):
        self.assertEqual(api._publication_date({'published-online': {'date-parts': [[2024, 2, 29]]},
                                               'published-print': {'date-parts': [[2025, 1, 1]]}}), '2024-02-29')

    def test_partial_dates_retain_their_precision(self):
        for parts, expected in (([2024], '2024'), ([2024, 6], '2024-06'), ([2024, 6, 2], '2024-06-02')):
            with self.subTest(parts=parts):
                self.assertEqual(api._publication_date({'published-print': {'date-parts': [parts]}}), expected)

    def test_licence_deposit_creation_and_index_are_not_publication(self):
        metadata = {'license': [{'start': {'date-parts': [[2026, 1, 1]]}}],
                    'created': {'date-parts': [[2020, 1, 1]]}, 'deposited': {'date-parts': [[2021]]},
                    'indexed': {'date-parts': [[2025]]}}
        self.assertEqual(api._publication_date(metadata), 'unknown')

    def test_invalid_online_date_can_use_valid_print(self):
        self.assertEqual(api._publication_date({'published-online': {'date-parts': [[2023, 2, 29]]},
                                               'published-print': {'date-parts': [[2024, 3]]}}), '2024-03')

    def test_malformed_or_invalid_dates_never_become_years(self):
        for value in (None, '2024', {}, {'date-parts': []}, {'date-parts': [[True]]},
                      {'date-parts': [[2024, 13]]}, {'date-parts': [[0]]},
                      {'date-parts': [[2024, 1, 1, 1]]}, {'date-parts': [[2024], [2025]]},
                      {'date-parts': [['2024']]}):
            with self.subTest(value=value):
                self.assertEqual(api._publication_date({'published-online': value}), 'unknown')


class RequestTests(unittest.TestCase):
    def setUp(self):
        class Timeout(Exception):
            pass
        class ConnectionError(Exception):
            pass
        self.requests = SimpleNamespace(get=Mock(), exceptions=SimpleNamespace(Timeout=Timeout, ConnectionError=ConnectionError))
        self.context = patch.dict('sys.modules', {'requests': self.requests})
        self.context.start()
        self.addCleanup(self.context.stop)

    def response(self, message, status=200):
        return SimpleNamespace(status_code=status, json=lambda: {'message': message})

    def test_doi_is_quoted_as_path_not_query(self):
        self.requests.get.return_value = self.response({'published-online': {'date-parts': [[2024]]}})
        self.assertEqual(api.get_publication_date(' 10.1234/a?x=2#frag '), '2024')
        self.requests.get.assert_called_once_with('https://api.crossref.org/works/10.1234/a%3Fx%3D2%23frag', timeout=15)

    def test_no_sleep_after_last_attempt(self):
        self.requests.get.side_effect = self.requests.exceptions.Timeout()
        with patch.object(api.time, 'sleep') as sleep:
            self.assertEqual(api.get_publication_date('10.1234/test', max_retries=3, delay=2), 'unknown')
            self.assertEqual(self.requests.get.call_count, 3)
            self.assertEqual(sleep.call_count, 2)

    def test_connection_failure_can_recover(self):
        self.requests.get.side_effect = [self.requests.exceptions.ConnectionError(),
                                        self.response({'published-print': {'date-parts': [[2025]]}})]
        with patch.object(api.time, 'sleep'):
            self.assertEqual(api.get_publication_date('10.1234/test'), '2025')

    def test_licence_only_response_is_unknown(self):
        self.requests.get.return_value = self.response({'license': [{'start': {'date-parts': [[2025]]}}]})
        self.assertEqual(api.get_publication_date('10.1234/test'), 'unknown')

    def test_non_success_http_does_not_supply_date(self):
        self.requests.get.return_value = self.response({'published-online': {'date-parts': [[2024]]}}, 404)
        self.assertEqual(api.get_publication_date('10.1234/test'), 'unknown')

    def test_bad_json_is_unknown(self):
        self.requests.get.return_value = SimpleNamespace(status_code=200, json=Mock(side_effect=ValueError('bad JSON')))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(api.get_publication_date('10.1234/test'), 'unknown')

    def test_invalid_limits_do_not_make_requests(self):
        for options in ({'max_retries': 0}, {'max_retries': True}, {'delay': float('nan')}, {'delay': -1}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                api.get_publication_date('10.1234/test', **options)
        self.requests.get.assert_not_called()

    def test_empty_doi_does_not_make_requests(self):
        self.assertEqual(api.get_publication_date(' '), 'unknown')
        self.requests.get.assert_not_called()


if __name__ == '__main__':
    unittest.main()
