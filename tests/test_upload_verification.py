"""Protect code-only contributions when copying into an existing fork."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from verify_upload import inspect


ROOT = Path(__file__).resolve().parents[1]


class UploadVerificationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='coremof-upload-test-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write(self, relative, text):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def track(self, relative):
        subprocess.run(['git', 'init', '--quiet', str(self.root)], check=True,
                       capture_output=True)
        subprocess.run(['git', '-C', str(self.root), 'add', '-f', '--', relative],
                       check=True, capture_output=True)

    def cli(self, *args):
        return subprocess.run([sys.executable, '-B', '-S', str(ROOT / 'verify_upload.py'),
                               '--root', str(self.root), *args],
                              capture_output=True, text=True)

    def test_untracked_ignored_local_assets_stay_excluded(self):
        self.write('.gitignore', '/local/\n')
        self.write('local/authorized.cif', 'data_local\n')
        result = inspect(self.root)
        self.assertEqual(result['status'], 'PASS_CODE_SURFACE')
        self.assertNotIn('local/authorized.cif', [r['path'] for r in result['files']])

    def test_tracked_ignored_database_table_is_rejected(self):
        self.write('.gitignore', '/CoREMOF/data/CR.json\n')
        self.write('CoREMOF/data/CR.json', '{}\n')
        self.track('CoREMOF/data/CR.json')
        self.assertTrue(any('local data/model/archive' in e
                            for e in inspect(self.root)['errors']))

    def test_tracked_extensionless_model_is_rejected(self):
        self.write('.gitignore', '/CoREMOF/models/cp_app/ensemble_models_*/\n')
        self.write('CoREMOF/models/cp_app/ensemble_models_legacy/300/model_0', 'weights')
        self.track('CoREMOF/models/cp_app/ensemble_models_legacy/300/model_0')
        self.assertTrue(any('local data/model/archive' in e
                            for e in inspect(self.root)['errors']))

    def test_parent_repository_check_is_scoped_to_requested_contribution(self):
        self.write('other_project/private.h5', 'weights')
        self.write('contribution/README.md', '# Contribution\n')
        self.track('other_project/private.h5')
        self.track('contribution/README.md')
        result = inspect(self.root / 'contribution')
        self.assertEqual(result['status'], 'PASS_CODE_SURFACE')
        self.assertEqual([r['path'] for r in result['files']], ['README.md'])

    def test_removed_worker_is_rejected(self):
        self.write('CoREMOF/_release_setc_worker.py', 'pass\n')
        self.assertTrue(any('removed checker implementation' in e
                            for e in inspect(self.root)['errors']))

    def test_unlisted_agent_documents_are_rejected(self):
        for name in ('skills/internal/SKILL.md', '.agents/skills/private/notes.md',
                     'docs/AGENTS.md', 'docs/copied/SKILL.md', '.codex/prompts/task.md'):
            with self.subTest(name=name):
                path = self.write(name, '# Internal instructions\n')
                self.assertTrue(any('consumer-document allowlist' in e
                                    for e in inspect(self.root)['errors']))
                path.unlink()

    def test_private_instructions_in_renamed_document_are_rejected(self):
        for name in ('docs/usage.md', 'README.rst', 'example.ipynb'):
            with self.subTest(name=name):
                path = self.write(name, 'Read /home/yuc/private/INSTRUCTIONS.md\n')
                self.assertTrue(any('private development instructions' in e
                                    for e in inspect(self.root)['errors']))
                path.unlink()

    def test_consumer_skill_is_allowed(self):
        self.write('.agents/skills/coremof-dataset-use/SKILL.md',
                   '# Read released metadata using the existing API\n')
        self.assertEqual(inspect(self.root)['status'], 'PASS_CODE_SURFACE')

    def test_allowed_skill_path_does_not_allow_private_content(self):
        self.write('.agents/skills/coremof-dataset-use/SKILL.md',
                   'Read latest_evidence_registry.json before release construction.\n')
        self.assertTrue(any('private development instructions' in e
                            for e in inspect(self.root)['errors']))

    def test_credentials_are_rejected(self):
        self.write('token.txt', 'ghp_' + 'a' * 40)
        self.assertTrue(any('credential' in e for e in inspect(self.root)['errors']))

    def test_symlinks_are_rejected(self):
        target = self.write('safe.txt', 'safe')
        (self.root / 'link.txt').symlink_to(target)
        self.assertTrue(any('symlink' in e for e in inspect(self.root)['errors']))

    def test_saved_notebook_results_are_rejected(self):
        self.write('example.ipynb', json.dumps({'cells': [{
            'cell_type': 'code', 'outputs': [{'text': 'result'}], 'execution_count': 1}]}))
        self.assertTrue(any('notebook output' in e for e in inspect(self.root)['errors']))

    def test_manifest_mismatch_fails_without_repair(self):
        self.write('example.py', 'answer = 42\n')
        self.assertEqual(self.cli('--write-manifest').returncode, 0)
        manifest = (self.root / 'GIT_UPLOAD_MANIFEST.json').read_bytes()
        self.write('example.py', 'answer = 43\n')
        result = self.cli()
        self.assertEqual(result.returncode, 1)
        self.assertIn('differ from GIT_UPLOAD_MANIFEST', result.stdout)
        self.assertEqual((self.root / 'GIT_UPLOAD_MANIFEST.json').read_bytes(), manifest)

    def test_manifest_refresh_then_read_only_check_pass(self):
        self.write('example.py', 'answer = 42\n')
        self.assertEqual(self.cli('--write-manifest').returncode, 0)
        self.assertEqual(self.cli().returncode, 0)


if __name__ == '__main__':
    unittest.main()
