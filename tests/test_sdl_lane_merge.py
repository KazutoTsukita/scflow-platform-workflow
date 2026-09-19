import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_scope_regressions import load_legacy_module


def member(text: str) -> bytes:
    return gzip.compress(f"@{text}\nACGT\n+\nIIII\n".encode())


class SdlLaneMergeTests(unittest.TestCase):
    """GSE254185 / SRR26934827: SRA never normalized the run, so SDL only offers the submitter's bcl2fastq lane files
    (L001/L002 x I1/R1/R2/R3). Strictly named lane/chunk files are concatenated per role in (lane, chunk) order;
    every other same-role duplication stays refused."""

    def setUp(self):
        self.sdl = load_legacy_module('download_ncbi_sdl_sources')
        self.contents = {
            'Sample2ATAC_S2_L001_R1_001.fastq.gz': member('L1R1'), 'Sample2ATAC_S2_L002_R1_001.fastq.gz': member('L2R1'),
            'Sample2ATAC_S2_L001_R2_001.fastq.gz': member('L1R2'), 'Sample2ATAC_S2_L002_R2_001.fastq.gz': member('L2R2'),
        }

    def payload(self, names):
        return {'result': [{'files': [dict(type='source', name=n, size=str(len(self.contents.get(n, b'x'))),
                                           locations=[dict(link=f'https://example.org/{n}')]) for n in names]}]}

    def run_one(self, root, names, download):
        with patch.object(self.sdl, 'fetch_sdl', return_value=self.payload(names)), \
                patch.object(self.sdl, 'download', side_effect=download):
            return self.sdl.run_one(dict(run_accession='SRR1'), root, 100, False, 'gzip', 0, 'full', 0)

    def fetch(self, url, output):
        output.write_bytes(self.contents[url.rsplit('/', 1)[1]])

    def test_lane_split_files_are_merged_per_role_in_lane_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = self.run_one(root, list(self.contents), self.fetch)
            finals = {r['path'].rsplit('/', 1)[1]: r for r in records if r['mode'] == 'fastq'}
            self.assertEqual(set(finals), {'SRR1_1.fastq.gz', 'SRR1_2.fastq.gz'})
            for r in finals.values():
                self.assertEqual(r['status'], 'downloaded', r)
                self.assertTrue(r['integrity'].startswith('lane_merged:2 parts'), r)
            self.assertEqual(gzip.decompress((root / 'SRR1_1.fastq.gz').read_bytes()), b'@L1R1\nACGT\n+\nIIII\n@L2R1\nACGT\n+\nIIII\n')
            self.assertEqual(gzip.decompress((root / 'SRR1_2.fastq.gz').read_bytes()), b'@L1R2\nACGT\n+\nIIII\n@L2R2\nACGT\n+\nIIII\n')
            self.assertEqual(sorted(p.name for p in root.iterdir()), ['SRR1_1.fastq.gz', 'SRR1_2.fastq.gz'])
            # part rows never look like inputs to the downstream coverage check
            self.assertTrue(all(r['status'].startswith('part_') and r['mode'] == 'fastq_part' for r in records if r not in finals.values()))

    def test_merged_output_is_reused_on_the_next_invocation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.run_one(root, list(self.contents), self.fetch)
            def unexpected(*args):
                self.fail('merged input was downloaded again')
            records = self.run_one(root, list(self.contents), unexpected)
            self.assertEqual({r['status'] for r in records}, {'skipped_existing'})

    def test_failed_part_fails_the_merged_role_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def flaky(url, output):
                if url.endswith('L002_R2_001.fastq.gz'):
                    raise subprocess.CalledProcessError(4, ['wget'])
                self.fetch(url, output)
            records = self.run_one(root, list(self.contents), flaky)
            finals = {r['path'].rsplit('/', 1)[1]: r for r in records if r['mode'] == 'fastq'}
            self.assertEqual(finals['SRR1_1.fastq.gz']['status'], 'downloaded')
            self.assertEqual(finals['SRR1_2.fastq.gz']['status'], 'failed')
            self.assertTrue(finals['SRR1_2.fastq.gz']['reason'].startswith('lane_merge_incomplete:'))
            self.assertFalse((root / 'SRR1_2.fastq.gz').exists())
            self.assertTrue((root / 'SRR1_2.lane01.part').exists())   # kept for resume

    def test_other_duplications_stay_refused(self):
        # lane set differs between roles; non-bcl2fastq names; two sample numbers
        for names in (['Sample2ATAC_S2_L001_R1_001.fastq.gz', 'Sample2ATAC_S2_L002_R1_001.fastq.gz', 'Sample2ATAC_S2_L001_R2_001.fastq.gz'],
                      ['read_R1.fastq.gz', 'read2_R1.fastq.gz'],
                      ['A_S1_L001_R1_001.fastq.gz', 'A_S2_L002_R1_001.fastq.gz', 'A_S1_L001_R2_001.fastq.gz', 'A_S2_L002_R2_001.fastq.gz']):
            with tempfile.TemporaryDirectory() as tmp:
                records = self.run_one(Path(tmp), names, self.fetch)
                self.assertEqual(len(records), 1, names)
                self.assertEqual(records[0]['status'], 'failed')
                self.assertTrue(records[0]['reason'].startswith('unsafe_sdl_source_roles:'), records[0])

    def test_single_file_per_role_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.contents['Sample2ATAC_S2_L001_I1_001.fastq.gz'] = member('I1')
            records = self.run_one(root, ['Sample2ATAC_S2_L001_R1_001.fastq.gz', 'Sample2ATAC_S2_L001_R2_001.fastq.gz', 'Sample2ATAC_S2_L001_I1_001.fastq.gz'], self.fetch)
            self.assertEqual(sorted(r['path'].rsplit('/', 1)[1] for r in records), ['SRR1_1.fastq.gz', 'SRR1_2.fastq.gz', 'SRR1_3.fastq.gz'])
            self.assertEqual({r['status'] for r in records}, {'downloaded'})


if __name__ == '__main__':
    unittest.main()
