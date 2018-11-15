import unittest
import random
import os
import subprocess
import numpy

_nuc = numpy.array(["A", "T", "C", "G"])


def randSeq(length):
    return ''.join(_nuc[[int(random.random()*4) for i in xrange(length)]])


class TestFaidx(unittest.TestCase):

    def setUp(self):
        # Create tiny genome fasta
        ref = open('temp.fa', 'w')
        ref.write('>x\n')
        ref.write(randSeq(20) + '\n')
        ref.write(randSeq(20) + '\n')
        ref.close()

    def test_faidx(self):
        # Call function
        cmd = ['python', 'MHcut/run_mhcut.py', '-ref', 'temp.fa']
        # dump = open('/dev/null')
        cmd_out = subprocess.check_output(cmd)
        # dump.close()
        self.assertTrue('Indexing completed' in cmd_out)
        self.assertTrue(os.path.isfile('temp.fa.fai'))

    def tearDown(self):
        # Remove input and output files
        os.remove('temp.fa')
        os.remove('temp.fa.fai')


class TestMHcutNoPam(unittest.TestCase):

    def setUp(self):
        # Create tiny genome fasta
        ref = open('temp.fa', 'w')
        ref.write('>x\n')
        ref.write(randSeq(20) + '\n')
        ref.write(randSeq(15) + 'ATCTA\n')
        ref.write('AGCTAG' + 'ATATATAT' + 'AGCTAG\n')
        ref.write('TTACT' + randSeq(15) + '\n')
        ref.write(randSeq(20) + '\n')
        ref.close()
        # Create deletions
        var = open('temp.tsv', 'w')
        var.write('chr\tstart\tend\tinfo\n')
        # Deletion with MH
        var.write('x\t41\t54\twith_mh\n')
        # Random deletion
        var.write('x\t30\t35\tno_mh\n')
        var.close()

    def test_mhcut(self):
        # Call function
        cmd = ['python', 'MHcut/run_mhcut.py', '-ref', 'temp.fa', '-var',
               'temp.tsv', '-out', 'temp']
        # dump = open('/dev/null')
        cmd_out = subprocess.check_output(cmd)
        # dump.close()
        self.assertTrue('Done' in cmd_out)
        self.assertTrue(os.path.isfile('temp-variants.tsv'))

    def tearDown(self):
        # Remove input and output files
        os.remove('temp.fa')
        os.remove('temp.fa.fai')
        os.remove('temp.tsv')
        os.remove('temp-cartoons.tsv')
        os.remove('temp-variants.tsv')
        os.remove('temp-guides.tsv')


if __name__ == '__main__':
    unittest.main()
