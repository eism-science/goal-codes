"""Phase 111 run logic: every check, submission, retrieval and analysis step of the hardware run. The run notebook only
calls this module, after checking this file's own hash against FREEZE_MANIFEST.txt, so the procedure that touches the
hardware is frozen and hashed, not notebook cells that Colab rewrites.

Safeguards (DESIGN.md, "Job identity" and Part C step 7):
- every file in FREEZE_MANIFEST.txt (except the notebook) matches; the operator's manifest SHA-256 (from the
  timestamped commit) matches; TIMESTAMP_CHECK.json records the verified OpenTimestamps attestation of that manifest,
  attested before any submission;
- the verification report passes for these exact files; QPY names match the plan in order; the power check is for these
  circuits; software versions match the build; a development build is refused;
- immediately before each submission: all hashes again, every instruction supported by the backend's target now, the
  step's shot count as frozen, the earlier step done, one job at most per experiment tag;
- a job record (job ID, experiment ID, tag, ordered names, circuit fingerprints, shots) is saved before waiting;
- retrieval checks tag, backend, record, per-circuit shots, and the submitted circuits themselves (from the service's
  record of the job inputs) against the frozen fingerprints. A job adopted from its tag alone, whose inputs the service
  does not return, is not analysed.
"""
import datetime
import hashlib
import json
import os

import p111_analyse as AN
import p111_verify as PV

NEEDED = ['DESIGN.md', 'p111_build.py', 'p111_verify.py', 'p111_power.py', 'p111_analyse.py', 'p111_runner.py',
          'p111_plan.json', 'p111_candidates.json', 'p111_circuits.qpy', 'p111_target.pkl', 'p111_verify_report.json',
          'p111_power.json', 'p111_backend_snapshot.json', 'testB_transpiled.qpy', 'testB_manifest.json',
          'testB_plan.json']
NOTEBOOK = 'phase111_hardware_run.ipynb'
STEPS = ('A2', 'pilot', 'B')
PILOT_SHOTS = 1000
PILOT_PAIRS = ['main_B_not_A/B', 'both/A']
BUDGET = {'A2': 120, 'pilot': 120, 'B': 300}          # max_execution_time, seconds


def now():
    return datetime.datetime.now(datetime.timezone.utc)


def read_manifest(folder):
    manifest = {}
    for line in open(os.path.join(folder, 'FREEZE_MANIFEST.txt')):
        if line.strip() and not line.startswith('#'):
            h, n = line.split(None, 1)
            if n.strip() in manifest:
                raise AssertionError(f'{n.strip()} listed twice in the manifest')
            manifest[n.strip()] = h
    return manifest


class Experiment:
    def __init__(self, folder, stamped_sha):
        from qiskit import qpy
        self.dir = folder
        self.sha = lambda n: hashlib.sha256(open(os.path.join(folder, n), 'rb').read()).hexdigest()
        self.manifest_sha = self.sha('FREEZE_MANIFEST.txt')
        assert stamped_sha.strip().lower() == self.manifest_sha, \
            f'this FREEZE_MANIFEST.txt ({self.manifest_sha}) is not the timestamped one'
        self.manifest = read_manifest(folder)
        unlisted = [n for n in NEEDED + [NOTEBOOK] if n not in self.manifest]
        assert not unlisted, f'not in FREEZE_MANIFEST.txt: {unlisted}'
        self.check_hashes()
        ts = self.load('TIMESTAMP_CHECK.json')
        assert ts and ts.get('manifest_sha256') == self.manifest_sha and ts.get('attested_utc'), \
            'TIMESTAMP_CHECK.json missing or not for this manifest: verify the OpenTimestamps proof first'
        self.attested = datetime.datetime.fromisoformat(ts['attested_utc'].replace('Z', '+00:00'))
        assert self.attested < now(), 'attestation time is in the future'
        self.timestamp = ts
        self.exp_id = 'p111-' + self.manifest_sha[:12]
        self.tags = {s: f'{self.exp_id}-{s}' for s in STEPS}

        self.plan = json.load(open(os.path.join(folder, 'p111_plan.json')))
        self.verify = json.load(open(os.path.join(folder, 'p111_verify_report.json')))
        self.power = json.load(open(os.path.join(folder, 'p111_power.json')))
        qsha = self.sha('p111_circuits.qpy')
        assert self.plan['live_target'], 'p111_plan.json is a development build; only the live build is run'
        assert self.verify['pass'] and self.verify['qpy_sha256'] == self.plan['qpy_sha256'] == qsha and \
            self.verify['plan_sha256'] == self.sha('p111_plan.json') and \
            self.verify['target_sha256'] == self.plan['target']['sha256'] == self.sha('p111_target.pkl'), \
            'verification report does not pass or is not for these files'
        assert self.power['qpy_sha256'] == qsha and self.power['shots_rule'] and \
            self.power['noise_source']['target_sha256'] == self.plan['target']['sha256'], \
            'power check missing, failed, or not for these circuits and calibration'
        b_circs = qpy.load(open(os.path.join(folder, 'p111_circuits.qpy'), 'rb'))
        b_names = self.plan['circuit_names']
        assert [c.name for c in b_circs] == b_names, 'circuit names in the QPY differ from the plan'
        assert [PV.fingerprint(c) for c in b_circs] == [self.verify['fingerprints'][n] for n in b_names]
        a_man = json.load(open(os.path.join(folder, 'testB_manifest.json')))
        a_circs = qpy.load(open(os.path.join(folder, 'testB_transpiled.qpy'), 'rb'))
        a_names = PV.testB_names(a_man)
        assert len(a_circs) == len(a_names) == self.verify['A2_circuits'] and \
            [m['name'] for m in self.verify['A2_mapping']] == a_names and \
            [PV.fingerprint(c) for c in a_circs] == self.verify['A2_fingerprints'], \
            'archived Test B circuits differ from the verified ones'
        pilot_phys = sorted({self.plan['pairs'][p]['physical_pair'] for p in PILOT_PAIRS})
        p_names = [f'{p}/{b}' for p in PILOT_PAIRS for b in ('XX', 'YY', 'ZZ')] + \
                  [f'calib/{pp}/{ab}' for pp in pilot_phys for ab in ('00', '01', '10', '11')]
        by_name = dict(zip(b_names, b_circs))
        self.circs = {'A2': a_circs, 'pilot': [by_name[n] for n in p_names], 'B': b_circs}
        self.names = {'A2': a_names, 'pilot': p_names, 'B': b_names}
        self.fps = {s: [PV.fingerprint(c) for c in self.circs[s]] for s in STEPS}
        self.shots = {'A2': a_man['shots'], 'pilot': PILOT_SHOTS, 'B': self.power['shots_rule']}
        self.width = {'A2': PV.V, 'pilot': 2, 'B': 2}
        self.service = self.backend = self.sampler_cls = None
        print(f'all {len(self.manifest) - 1} frozen files match the timestamped manifest (attested '
              f'{self.attested.isoformat()}) | experiment {self.exp_id}')
        print(' | '.join(f'{s}: {len(self.circs[s])} circuits x {self.shots[s]} shots' for s in STEPS))

    # ------------------------------------------------------------ files
    def check_hashes(self):
        missing = [n for n in self.manifest if n != NOTEBOOK and not os.path.exists(os.path.join(self.dir, n))]
        assert not missing, f'put these in {self.dir}: {missing}'
        wrong = [n for n in self.manifest if n != NOTEBOOK and self.manifest[n] != self.sha(n)]
        assert not wrong, f'do not match FREEZE_MANIFEST.txt: {wrong}'
        assert self.sha('FREEZE_MANIFEST.txt') == self.manifest_sha, 'FREEZE_MANIFEST.txt changed'

    def load(self, name):
        p = os.path.join(self.dir, name)
        if not os.path.exists(p):
            return None
        with open(p) as f:
            return json.load(f)

    def save(self, name, obj):
        """Atomic write, read back; stops (with the content printed) if it cannot be saved."""
        data = json.dumps(obj, indent=1, default=float)
        path = os.path.join(self.dir, name)
        try:
            with open(path + '.tmp', 'w') as f:
                f.write(data)
            os.replace(path + '.tmp', path)
            ok = open(path).read() == data
        except OSError:
            ok = False
        if not ok:
            raise RuntimeError(f'could not save {name}; record it by hand: {data[:500]}')
        print(f'saved {name}')

    # ------------------------------------------------------------ backend
    def connect(self, service, backend_name, sampler_cls):
        import qiskit
        import qiskit_ibm_runtime
        self.versions = {'qiskit': qiskit.__version__, 'qiskit_ibm_runtime': qiskit_ibm_runtime.__version__}
        assert all(self.versions[k] == self.plan['versions'][k] for k in self.versions), \
            f"versions {self.versions} differ from the build {self.plan['versions']}"
        assert backend_name == self.plan['backend'], f"the build is for {self.plan['backend']}"
        self.service, self.sampler_cls = service, sampler_cls
        self.backend = service.backend(backend_name)
        print('using', self.backend.name, '| pending jobs:', self.backend.status().pending_jobs)
        for s in STEPS:
            print(f'existing {s} jobs for this experiment:', [j.job_id() for j in self.tagged(s)])

    def tagged(self, step):
        try:
            return list(self.service.jobs(job_tags=[self.tags[step]], limit=20))
        except Exception as e:
            raise RuntimeError(f'could not query IBM for {self.tags[step]} jobs ({e}); do not submit until this works')

    def target_failures(self, step):
        bad = []
        for i, c in enumerate(self.circs[step]):
            for ins in c.data:
                if ins.operation.name == 'barrier':
                    continue
                qs = tuple(c.find_bit(q).index for q in ins.qubits)
                if not self.backend.target.instruction_supported(ins.operation.name, qs,
                                                                  parameters=list(ins.operation.params)):
                    bad.append((self.names[step][i], ins.operation.name, qs))
                    break
        return bad

    def target_report(self):
        for s in STEPS:
            bad = self.target_failures(s)
            print(f'{s}:', 'supported by today\'s target' if not bad else f'NOT SUPPORTED {bad[:3]}')

    # ------------------------------------------------------------ submission
    def record(self, step, jid, recovered):
        return {'job_id': jid, 'experiment_id': self.exp_id, 'step': step, 'tag': self.tags[step],
                'backend': self.backend.name, 'shots': self.shots[step], 'circuit_names': self.names[step],
                'fingerprints': self.fps[step], 'manifest_sha256': self.manifest_sha, 'versions': self.versions,
                'timestamp_check': self.timestamp, 'recovered_from_tag': recovered,
                'recorded_utc': now().isoformat()}

    def submit_once(self, step):
        """Submit the step once; returns the job ID (existing or new). Every check runs here, just before submitting."""
        assert self.backend is not None, 'connect first'
        rec, existing = self.load(f'{step}_jobid.json'), self.tagged(step)
        if rec:
            assert rec.get('experiment_id') == self.exp_id, f'{step}_jobid.json belongs to another experiment: move it'
        ids = sorted({j.job_id() for j in existing} | ({rec['job_id']} if rec else set()))
        if len(ids) > 1:
            raise RuntimeError(f'{step}: several jobs for this experiment {ids}; resolve by hand before going on')
        if ids:
            if not rec:
                self.save(f'{step}_jobid.json', self.record(step, ids[0], recovered=True))
            print(f'{step}: already submitted, job {ids[0]}; not resubmitting')
            return ids[0]
        prior = {'pilot': 'A2', 'B': 'pilot'}.get(step)
        assert prior is None or self.load(f'{prior}_raw_counts.json'), f'retrieve {prior} first'
        self.check_hashes()
        bad = self.target_failures(step)
        assert not bad, f'{step}: not supported by today\'s target {bad[:3]}: stop (DESIGN.md, Part C step 6)'
        try:
            print('IBM usage report:', self.service.usage())
        except Exception as e:
            print('usage report unavailable:', e, '- check the remaining time on the IBM Quantum dashboard')
        n = len(self.circs[step])
        if input(f'{step}: {n} circuits x {self.shots[step]} shots, up to {BUDGET[step]} QPU seconds. '
                 'Type YES to submit: ').strip() != 'YES':
            raise RuntimeError('not confirmed; nothing submitted')
        assert now() > self.attested
        s = self.sampler_cls(mode=self.backend)
        s.options.max_execution_time = BUDGET[step]
        s.options.environment.job_tags = [self.tags[step]]
        s.options.dynamical_decoupling.enable = False
        s.options.twirling.enable_gates = False
        s.options.twirling.enable_measure = False
        job = s.run(self.circs[step], shots=self.shots[step])
        print(f'{step.upper()} JOB ID (write it down):', job.job_id())
        self.save(f'{step}_jobid.json', self.record(step, job.job_id(), recovered=False))
        return job.job_id()

    # ------------------------------------------------------------ retrieval
    def submitted_fingerprints(self, job):
        """Fingerprints of the circuits the service recorded as this job's input, or None if it does not return them."""
        try:
            pubs = job.inputs['pubs']
            return [PV.fingerprint(p.circuit if hasattr(p, 'circuit') else p[0]) for p in pubs]
        except Exception as e:
            print('the service did not return the job inputs:', repr(e)[:200])
            return None

    def retrieve(self, step):
        """Raw counts by circuit name, after the identity and technical checks. Saves them."""
        names = self.names[step]
        rec = self.load(f'{step}_jobid.json')
        assert rec and rec['experiment_id'] == self.exp_id and rec['circuit_names'] == names and \
            rec['fingerprints'] == self.fps[step] and rec['shots'] == self.shots[step], \
            f'{step}: job record missing or not for these frozen circuits'
        job = self.service.job(rec['job_id'])
        assert self.tags[step] in (job.tags or []), f"{step}: job does not carry the tag {self.tags[step]}"
        assert job.backend().name == self.plan['backend'], f'{step}: job ran on {job.backend().name}'
        fps = self.submitted_fingerprints(job)
        if fps is None:
            assert not rec['recovered_from_tag'], \
                f'{step}: job adopted from its tag and its circuits cannot be checked; resolve by hand'
            payload = 'not returned by the service; job submitted by this procedure (record saved at submission)'
        else:
            assert fps == self.fps[step], f'{step}: the submitted circuits differ from the frozen ones'
            payload = 'submitted circuits match the frozen fingerprints'
        print(f'{step}: {payload}')
        res = job.result()
        assert len(res) == len(names), f'{step}: {len(res)} results for {len(names)} circuits'
        counts = {}
        for i, n in enumerate(names):
            m = res[i].data.meas
            assert m.num_shots == self.shots[step], f'{step}: {n} has {m.num_shots} shots'
            c = m.get_counts()
            width = len(next(iter(c)).replace(' ', ''))
            assert width == self.width[step], f'{step}: unexpected bit width {width} in {n}'
            counts[n] = c
        try:
            usage = job.usage()
        except Exception:
            usage = job.metrics().get('usage')
        self.save(f'{step}_raw_counts.json', {'job_id': rec['job_id'], 'experiment_id': self.exp_id,
                                              'backend': self.backend.name, 'usage': usage, 'payload_check': payload,
                                              'counts': counts})
        print(f'{step}: retrieved, usage {usage}')
        return counts

    # ------------------------------------------------------------ analysis
    def analyse(self, step):
        counts = self.retrieve(step)
        if step == 'A2':
            out = AN.analyse_A2(counts, json.load(open(os.path.join(self.dir, 'testB_plan.json'))))
            for k, v in out['witnesses'].items():
                print(f"  {k:16s} E {v['E']:.3f} (se {v['E_se']:.3f})")
        else:
            out = AN.analyse_B(counts, self.plan)
            for n, r in out['pairs'].items():
                fc = '' if r['F_corrected'] is None else f"  corrected {r['F_corrected']:.3f}"
                print(f"  {n:16s} F {r['F']:.3f} (se {r['F_se']:.3f}, ideal {r['ideal_F']}){fc}")
            for p, b in out['baselines'].items():
                print(f"  baseline {p}: F {b['F']:.3f} (se {b['F_se']:.3f})")
        if step != 'pilot':                           # the pilot is reported, not tested
            for k, t in out['tests'].items():
                print(f"  {k}: {t['value']:+.3f}, z {t['z']:.1f}, {'holds' if t['holds'] else 'FAILS'}")
            print({'A2': 'PART A2: ', 'B': 'PART B: '}[step] +
                  {('A2', True): 'REPLICATED', ('A2', False): 'NOT REPLICATED', ('B', True): 'PASS',
                   ('B', False): 'FAIL'}[(step, out['pass'])])
        self.save(f'{step}_result.json', out)
        return out

    def job_record(self, since_days=1):
        jobs = self.service.jobs(limit=100, created_after=now() - datetime.timedelta(days=since_days))
        rec = [{'job_id': j.job_id(), 'status': str(j.status()), 'tags': getattr(j, 'tags', None),
                'created': str(getattr(j, 'creation_date', ''))} for j in jobs]
        for r in rec:
            print(r)
        self.save('job_record.json', rec)
