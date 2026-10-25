"""Offline end-to-end test of phase111_hardware_run.ipynb and p111_runner.py: IBM service and sampler replaced by fakes
that run the circuits on Aer with the FakeFez noise model.

Usage:
    python offline_notebook_test.py <phase111 folder> <empty work folder>

Checks: manifest, timestamp record and file verification; target check; one submission per step; no resubmission on
rerun or after lost local records; a stale job of another experiment ignored; a duplicate job, a foreign local record,
submitted circuits that differ from the frozen ones, a job adopted from its tag whose inputs cannot be checked, a
changed frozen file and a wrong manifest hash each stop the run.
The development build is marked as live here (with a stand-in calibration snapshot) only so the runner accepts it."""
import builtins, datetime, hashlib, json, os, shutil, sys, uuid

SRC = sys.argv[1]
WORK = sys.argv[2]
nb = json.load(open(os.path.join(SRC, 'phase111_hardware_run.ipynb'), encoding='utf8'))
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()

# ---- frozen files into the work folder, with a test manifest and timestamp record
D = os.path.join(WORK, 'phase111')
os.makedirs(D, exist_ok=True)
for f in ('DESIGN.md', 'p111_build.py', 'p111_verify.py', 'p111_power.py', 'p111_analyse.py', 'p111_runner.py',
          'p111_candidates.json', 'p111_circuits.qpy', 'p111_target.pkl', 'p111_power.json',
          'phase111_hardware_run.ipynb'):
    shutil.copy(os.path.join(SRC, f), D)
for f in ('testB_transpiled.qpy', 'testB_manifest.json'):
    shutil.copy(os.path.join(SRC, '..', 'phase107', 'IBM run data', f), D)
shutil.copy(os.path.join(SRC, '..', 'phase107', 'testB_plan.json'), D)
plan = json.load(open(os.path.join(SRC, 'p111_plan.json')))
plan['live_target'] = True
plan['backend'] = 'ibm_fez'
json.dump(plan, open(os.path.join(D, 'p111_plan.json'), 'w'), indent=1)
rep = json.load(open(os.path.join(SRC, 'p111_verify_report.json')))
rep['plan_sha256'] = sha(os.path.join(D, 'p111_plan.json'))
json.dump(rep, open(os.path.join(D, 'p111_verify_report.json'), 'w'), indent=1)
json.dump({'stand-in': True}, open(os.path.join(D, 'p111_backend_snapshot.json'), 'w'))
with open(os.path.join(D, 'FREEZE_MANIFEST.txt'), 'w') as m:
    m.write('# test manifest\n')
    for f in sorted(os.listdir(D)):
        if f != 'FREEZE_MANIFEST.txt':
            m.write(f"{sha(os.path.join(D, f))}  {f}\n")
MAN_SHA = sha(os.path.join(D, 'FREEZE_MANIFEST.txt'))
json.dump({'manifest_sha256': MAN_SHA, 'attested_utc': '2026-10-24T12:00:00Z', 'bitcoin_block': 0,
           'verified_with': 'test', 'verified_utc': '2026-10-24T13:00:00Z'},
          open(os.path.join(D, 'TIMESTAMP_CHECK.json'), 'w'))
os.chdir(WORK)

# ---- fakes
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeFez

BACKEND = FakeFez()
SIM = AerSimulator.from_backend(BACKEND)
JOBS = {}
SUBMISSIONS = []
HIDE_INPUTS = set()


class _Meas:
    def __init__(self, counts):
        self._c = counts
        self.num_shots = sum(counts.values())

    def get_counts(self):
        return self._c


class FakeJob:
    def __init__(self, circs, shots, tags, run=True):
        self._id = 'fake-' + uuid.uuid4().hex[:8]
        self.tags = list(tags)
        self.creation_date = datetime.datetime.utcnow()
        self._circs = list(circs)
        if run:
            res = SIM.run(circs, shots=shots, seed_simulator=len(JOBS) + 11).result()
            self._res = [type('R', (), {'data': type('D', (), {'meas': _Meas(res.get_counts(i))})()})()
                         for i in range(len(circs))]

    @property
    def inputs(self):
        if self._id in HIDE_INPUTS:
            raise RuntimeError('inputs not available')
        return {'pubs': [(c,) for c in self._circs]}

    def job_id(self):
        return self._id

    def backend(self):
        return type('B', (), {'name': 'ibm_fez'})()

    def result(self):
        return self._res

    def usage(self):
        return 1.0

    def status(self):
        return 'DONE'


class _NS:
    def __init__(self):
        self.__dict__['_d'] = {}

    def __getattr__(self, k):
        return self._d.setdefault(k, _NS())

    def __setattr__(self, k, v):
        self._d[k] = v


class FakeSampler:
    def __init__(self, mode):
        assert mode is BACKEND
        self.options = _NS()

    def run(self, circs, shots):
        tags = self.options.environment._d.get('job_tags', [])
        o = self.options
        assert o.dynamical_decoupling._d.get('enable') is False and o.twirling._d.get('enable_gates') is False
        job = FakeJob(circs, shots, tags)
        JOBS[job.job_id()] = job
        SUBMISSIONS.append((tags, len(circs), shots))
        return job


class FakeService:
    def __init__(self, channel, token, instance):
        pass

    def backend(self, name):
        BACKEND.status = lambda: type('S', (), {'pending_jobs': 0})()
        return BACKEND

    def jobs(self, job_tags=None, limit=20, created_after=None):       # IBM: a job must carry all the given tags
        return [j for j in JOBS.values() if job_tags is None or set(job_tags) <= set(j.tags)]

    def job(self, jid):
        return JOBS[jid]

    def usage(self):
        return {'usage_remaining_seconds': 600}


stale = FakeJob([], 0, ['phase111-B-full'], run=False)       # an earlier design's generic tag: must be ignored
JOBS[stale.job_id()] = stale

builtins.input = lambda prompt='': (print(prompt + '<answer>'), MAN_SHA if 'SHA-256' in prompt else 'YES')[1]
os.environ['IBM_QUANTUM_TOKEN'] = 'test'
ns = {'__name__': '__main__'}
code_cells = [''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code']
cell = lambda text: next(i for i, c in enumerate(code_cells) if text in c)


def run(i):
    src = '\n'.join(l for l in code_cells[i].splitlines() if not l.lstrip().startswith('!'))
    src = src.replace('from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2', '')
    ns.update({'QiskitRuntimeService': FakeService, 'SamplerV2': FakeSampler})
    print(f'--- cell {i}')
    exec(compile(src, f'cell{i}', 'exec'), ns)


def expect_stop(i, what):
    try:
        run(i)
    except (AssertionError, RuntimeError) as e:
        print(f'stopped as expected ({what}): {str(e)[:150]}')
        return
    raise SystemExit(f'FAILED: no stop for {what}')


for i in range(len(code_cells)):
    run(i)
print('\nsubmissions:', SUBMISSIONS)
assert len(SUBMISSIONS) == 3, 'expected exactly three submissions'
submit = {s: cell(f"submit_once('{s}')") for s in ('A2', 'pilot', 'B')}
analyse_B = cell("analyse('B')")
EXP = ns['EXP']
# ---- reruns, and lost local records: nothing resubmitted, records rebuilt from the experiment tags
for i in submit.values():
    run(i)
for s in ('A2', 'pilot', 'B'):
    os.remove(os.path.join(D, f'{s}_jobid.json'))
for i in submit.values():
    run(i)
assert all(json.load(open(os.path.join(D, f'{s}_jobid.json')))['recovered_from_tag'] for s in ('A2', 'pilot', 'B'))
run(analyse_B)                                              # recovered record, inputs returned: analysed
assert len(SUBMISSIONS) == 3, f'resubmission happened: {SUBMISSIONS}'
jid_B = json.load(open(os.path.join(D, 'B_jobid.json')))['job_id']
# ---- recovered record and the service does not return the inputs: not analysed
HIDE_INPUTS.add(jid_B)
expect_stop(analyse_B, 'job adopted from its tag, inputs not returned')
HIDE_INPUTS.discard(jid_B)
# ---- the job's submitted circuits differ from the frozen ones
JOBS[jid_B]._circs = list(reversed(JOBS[jid_B]._circs))
expect_stop(analyse_B, 'submitted circuits differ')
JOBS[jid_B]._circs = list(reversed(JOBS[jid_B]._circs))
# ---- a second job with this experiment's B tag
dup = FakeJob([], 0, [EXP.tags['B']], run=False)
JOBS[dup.job_id()] = dup
expect_stop(submit['B'], 'two jobs carry the experiment tag')
del JOBS[dup.job_id()]
# ---- a local record from another experiment
rec = json.load(open(os.path.join(D, 'B_jobid.json')))
json.dump(dict(rec, experiment_id='p111-000000000000'), open(os.path.join(D, 'B_jobid.json'), 'w'))
expect_stop(submit['B'], 'foreign local record')
expect_stop(analyse_B, 'foreign local record at retrieval')
json.dump(rec, open(os.path.join(D, 'B_jobid.json'), 'w'))
# ---- a frozen file changed after loading: caught just before a submission
os.remove(os.path.join(D, 'B_jobid.json'))
del JOBS[jid_B]
with open(os.path.join(D, 'DESIGN.md'), 'a') as f:
    f.write('\nchanged\n')
expect_stop(submit['B'], 'frozen file changed before submission')
# ---- a wrong manifest hash typed in
builtins.input = lambda prompt='': '0' * 64
expect_stop(cell('RUN.Experiment('), 'manifest hash not the timestamped one')
assert len(SUBMISSIONS) == 3
print('\nNOTEBOOK OFFLINE TEST PASS: 3 submissions; no resubmission after reruns or lost records; stale job ignored; '
      'unverifiable or mismatched payload, duplicate job, foreign record, changed file and wrong manifest hash all stop '
      'the run')
