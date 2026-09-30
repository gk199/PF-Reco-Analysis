"""Merge the exact successful input jobs; publish only validated results."""
import json
import os
from pathlib import Path
import subprocess
import time
from root_counts import tree_counts


def command(args):
    print('Running:', ' '.join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), check=True)


def copy_file(source, target):
    for attempt in range(3):
        try:
            command(['xrdcp', '--force', '--cksum', 'adler32', '--posc', source, target])
            return
        except subprocess.CalledProcessError:
            if attempt == 2:
                raise
            time.sleep(10)


def run():
    host = os.environ['EOS_HOST'].rstrip('/')
    base = os.environ['EOS_RUN'].rstrip('/')
    save_reco = os.environ['SAVE_RECO'] == '1'
    diagnostics = Path('diagnostics')
    jobs = [row.split() for row in Path('jobs.txt').read_text().splitlines() if row.strip()]
    inputs = Path('selected_inputs.txt').read_text().splitlines()
    if len(jobs) != len(inputs) or not jobs:
        raise RuntimeError('Job/input manifest mismatch')
    records = []
    for index, tag in jobs:
        local = diagnostics / f'{tag}.done.json'
        copy_file(f'{host}/{base}/completed/{tag}.done.json', local)
        data = json.loads(local.read_text())
        if (data['tag'] != tag or data['input'] != inputs[int(index)] or
                data['max_events'] != int(os.environ['MAX_EVENTS']) or
                data['threshold'] != float(os.environ['TIMING_THRESHOLD']) or
                data['save_reco'] != save_reco):
            raise RuntimeError('Completion record does not match this submission: ' + tag)
        records.append(data)

    parts_to_remove, results = [], {}
    for variant in ('standardPF', 'seedTimingPF'):
        folder = variant if variant == 'standardPF' else variant + '_' + os.environ['TIMING_TAG']
        part_paths = [f'{base}/parts/{folder}/pfObjectsNtuple_{tag}_{variant}.root' for _, tag in jobs]
        part_urls = [host + '/' + path for path in part_paths]
        expected = {}
        schema = set(records[0]['trees'][variant])
        for record in records:
            counts = record['trees'][variant]
            if set(counts) != schema:
                raise RuntimeError('Inconsistent tree names across input ntuples')
            for name, count in counts.items():
                expected[name] = expected.get(name, 0) + count
        list_file = diagnostics / f'ntuples_{variant}.txt'
        list_file.write_text('\n'.join(part_urls) + '\n')
        merged = Path(f'pfObjectsNtuple_{variant}.root')
        # Only TFileService ntuples go through hadd. Never hadd CMS EDM files.
        command(['hadd', '-f', '-n', '32', merged, '@' + str(list_file)])
        actual = tree_counts(merged)
        if actual != expected:
            raise RuntimeError(f'Merged entry counts differ: {actual} != {expected}')
        final_url = f'{host}/{base}/{merged.name}'
        copy_file(merged, final_url)
        if tree_counts(final_url) != expected:
            raise RuntimeError('Remote merged ntuple verification failed')
        results[variant] = {'ntuple': final_url, 'trees': actual}
        parts_to_remove.extend(part_paths)
        merged.unlink()

        if save_reco:
            edm_paths = [f'{base}/parts/{folder}/pf_reReco_{tag}_{variant}.root' for _, tag in jobs]
            edm_urls = [host + '/' + path for path in edm_paths]
            merged_edm = Path(f'pf_reReco_{variant}.root')
            config = diagnostics / f'merge_edm_{variant}_cfg.py'
            config.write_text(
                'import FWCore.ParameterSet.Config as cms\n'
                'process = cms.Process("MERGE")\n'
                'process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))\n'
                'process.source = cms.Source("PoolSource",\n'
                f'    fileNames=cms.untracked.vstring({edm_urls!r}),\n'
                '    duplicateCheckMode=cms.untracked.string("checkAllFilesOpened"))\n'
                'process.out = cms.OutputModule("PoolOutputModule",\n'
                f'    fileName=cms.untracked.string({str(merged_edm)!r}),\n'
                '    outputCommands=cms.untracked.vstring("keep *"))\n'
                'process.end = cms.EndPath(process.out)\n'
            )
            # Use CMSSW-aware merging to preserve event provenance/references.
            command(['cmsRun', '-j', diagnostics / f'merge_{variant}.xml', config])
            command(['edmFileUtil', merged_edm])
            expected_events = sum(record['edm_events'][variant] for record in records)
            if tree_counts(merged_edm).get('Events') != expected_events:
                raise RuntimeError('EDM event count differs from the input sum; check for duplicate event IDs.')
            final_edm = f'{host}/{base}/{merged_edm.name}'
            copy_file(merged_edm, final_edm)
            if tree_counts(final_edm).get('Events') != expected_events:
                raise RuntimeError('Remote merged EDM event count differs from the input sum.')
            results[variant]['reco'] = final_edm
            parts_to_remove.extend(edm_paths)
            merged_edm.unlink()

    # No parts are removed until BOTH algorithms' final files are validated.
    summary = diagnostics / 'merge_complete.json'
    summary.write_text(json.dumps({'input_files': inputs, 'outputs': results}, indent=2) + '\n')
    copy_file(summary, f'{host}/{base}/merge_complete.json')
    if os.environ['KEEP_PARTS'] == '0':
        leftovers = []
        for path in parts_to_remove:
            # Only exact paths constructed for this submission are removed.
            try:
                command(['xrdfs', host, 'rm', path])
            except subprocess.CalledProcessError:
                leftovers.append(path)
        if leftovers:
            (diagnostics / 'cleanup_remaining.txt').write_text('\n'.join(leftovers) + '\n')
            print('Final files are complete, but some temporary parts could not be removed; see cleanup_remaining.txt.', flush=True)
    print('All selected inputs merged successfully.', flush=True)


if __name__ == '__main__':
    run()
