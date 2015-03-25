import subprocess
import sys
from itertools import *
from os.path import abspath, join, dirname, isfile, basename
from os import listdir

from GoTools import *
from GeoUtils.IO import InputFile


def summary_patch(p, out):
    out.write('%s\n' % type(p))
    curve = type(p) is Curve

    kts = p.GetKnots(with_multiplicities=True)
    if curve:
        kts = [kts]
    for i, ks in enumerate(kts):
        out.write('kmult dir %i: ' % (i+1))
        out.write(' '.join(str(len(list(b))) for _, b in groupby(ks)) + '\n')


def summary_g2(f, out):
    patches = ReadG2(f)
    if type(patches) is not list:
        patches = [patches]

    for i, patch in enumerate(patches):
        out.write('patch %i\n' % (i+1))
        summary_patch(patch, out)


def summary_xinp(f, out):
    xinp = InputFile(f)
    for ts_name in xinp.GetTopologySets():
        out.write('topset: %s\n' % ts_name)
        ts = xinp.GetTopologySet(ts_name, convention='ifem')
        keys = sorted(ts.keys())
        for k in keys:
            out.write('%i vertices: %s edges: %s faces: %s\n'
                      % (k, ' '.join(str(v) for v in ts[k].vertex),
                         ' '.join(str(e) for e in ts[k].edge),
                         ' '.join(str(f) for f in ts[k].face)))


def summary_dir(path, out):
    files = [f for f in listdir(path) if f.endswith('.g2')]
    files.sort()

    for f in files:
        out.write(f + '\n')
        summary_g2(join(path, f), out)

        xinp = join(path, f.split('.')[0] + '.xinp')
        if isfile(xinp):
            out.write(basename(xinp) + '\n')
            summary_xinp(xinp, out)
                

if __name__ == '__main__':
    root = dirname(abspath(__file__))

    testfiles = [f for f in listdir(join(root, 'tests')) if f.endswith('.yaml')]
    testfiles.sort()

    script = join(root, 'nrel.py')

    try:
        nprocs = int(sys.argv[1])
    except IndexError:
        nprocs = 4

    total_ok = True

    for testfile in testfiles:
        sys.stdout.write(testfile + '... ')
        sys.stdout.flush()

        proc = subprocess.Popen(['geoModeler', script,
                                 'paramfile=' + join(root, 'tests', testfile),
                                 'out=testres', 'nprocs=1', 'nprocs_mg=%i' % nprocs],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        success = stderr == '' and 'WARNING' not in stdout and 'ERROR' not in stdout

        if not success:
            print 'FAILED'
            total_ok = False
            continue

        summary = join(root, 'testres', 'summary')
        compare = join(root, 'tests', testfile.split('.')[0] + '.out')
        with open(summary, 'w') as f:
            summary_dir(join(root, 'testres'), f)
        proc = subprocess.Popen(['diff', summary, compare],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        success = stdout == '' and stderr == ''

        print 'OK' if success else 'FAILED'

        if not success:
            print stdout if stderr == '' else stderr
        total_ok = total_ok and success

    sys.exit(0 if total_ok else 1)
