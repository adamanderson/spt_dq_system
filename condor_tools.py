import os.path
import subprocess

def condor_submit(errdir, outdir, logdir, proxycert):
    submit = [
        'executable = autoproc_job.sh',
        'universe = vanilla',
        'error = {}'.format(os.path.join(errdir, 'error.log')),
        'output = {}'.format(os.path.join(outdir, 'out.log')),
        'log = {}'.format(os.path.join(logdir, 'condor.log')),
        '+WANT_RCC_ciconnect = True',
        '+ProjectName = "spt-all"',
        'request_cpus = 1',
        'request_memory = 2GB',
        'request_disk = 10GB',
        'transfer_output_files = ""',
        'should_transfer_files = YES',
        '(OpSysAndVer == "CentOS6" || OpSysAndVer == "RedHat6" || OpSysAndVer == "SL6")',
        'x509userproxy = {}'.format(proxycert),
        'when_to_transfer_output = ON_EXIT',
        'transfer_input_files = plotter.sh, plotter.py, {}'.format(proxycert),
        'transfer_executable = True',
        'queue']
    try:
        result = subprocess.check_output('condor_submit', input='\n'.join(submit).encode('ascii'))
        return int(result.split()[-1][:-1])
    except subprocess.CalledProcessError as e:
        print('Error submitting job: %s' % e.output)
        return None
