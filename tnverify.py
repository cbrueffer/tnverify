#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Christian Brueffer (ORCID: 0000-0002-3826-0989)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.

import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import random
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


col_format = 8
col_sample_start = 9


def vcf2flag(x):
    if x == ".":
        return "NA"
    if not ":" in x:
        print "Could not parse VCF entry", x
        raise NotImplementedError
    x = x.split(":")[0]
    if x == ".":
        return "0"
    if x == "0/0":
        return ("0")
    # if not ":" in x
    #     print "Could not parse VCF entry", x
    #     raise NotImplementedError
    if x.startswith("0/"):
        return ("1")
    if len(x) > 2 and x[0] == x[2]:
        return ("2")
    print x
    raise NotImplementedError


def init_logger(level=None, logfile=None):
    """Initializes a logger.  If logfile is given, the logger will log to the
    file, otherwise to the console."""

    import logging

    if level is None:
        level = logging.INFO

    if logfile:
        handler = logging.FileHandler(logfile)
    else:
        handler = logging.StreamHandler()

    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger = logging.getLogger("")
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger


def exec_and_log(command, logger=None):
    """Executes a given commandline and logs the command output to the given logger."""
    import subprocess

    try:
        s = subprocess.Popen(command.split(),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        while True:
            line = s.stdout.readline()
            if not line:
                break
            if logger:
                logger.info(line.rstrip())
    except Exception as e:
        raise e


def get_file_dims(f):
    """Determines the number of rows, columns and comments
    in the given file."""
    rows = columns = comments = 0
    for line in f:
        if line.startswith("#"):
            comments += 1
        else:
            rows += 1
            c = line.rstrip().split("\t")
            if columns == 0:
                columns = len(c[col_sample_start:])

    f.seek(0)

    return rows, columns, comments


class tnverify:

    def __init__(self, snpfile, samplefile, logger=None):
        if logger:
            self.logger = logger
        else:
            self.logger = init_logger()

        self.sample_paths, self.sample_labels = self.read_samplefile(samplefile)
        #outfile = ""
        #call_snps(samplelist, reference, regions, outfile)
    
        self.flagmatrix, self.samplenames = self.vcf2ndarray(snpfile,
                                                             include_random=False)

        self.clusterplot(self.flagmatrix)

    def call_snps(self, samples, reference, regions, outfile):
        # samtools mpileup -IguDB -f $ref -l $regions ${BamDir}/${samples} | bcftools view -vcg - > ${OutDir}/result-contralat.txt"
        samtools_cmd = "samtools mpileup -IguDB -f %s -l %s %s | bcftools view -vcg - > %s" % (reference, regions, " ".join(samples), outfile)

        exec_and_log(samtools_cmd, this.logger)

    def generate_random_sample(self, length):
        """Return a sample vector with random variant calls."""
        return [random.randint(0, 2) for x in range(length)]

    def vcf2ndarray(self, vcffile, include_random=True):
        vcfmatrix = None
        samplenames = None
        valid_count = 0
        with open(vcffile) as vcf:
            nrows, ncols, ncomments = get_file_dims(vcf)

            # +1 to leave space for a sample with random variant calls
            if include_random:
                vcfmatrix = np.ndarray((nrows, ncols+1))
            else:
                vcfmatrix = np.ndarray((nrows, ncols))
            for k, line in enumerate(vcf):
                if line.startswith("##"):
                    continue
                cols = line.rstrip().split("\t")
                samples = cols[col_sample_start:]
                format = cols[col_format]

                if line.startswith("#"):
                    # print samples
                    samplenames = samples
                #    print "\t".join(sampleNames)
                    continue

                # make sure that GT is always first
                if not format[:2] == "GT":
                    #print >> sys.stderr, "Skipping entry with format %s (no genotype found)" % format
                    continue
    
                # convert vcf entries to simple flags
                flags = [vcf2flag(x) for x in samples]
                # captures also tri-allelic SNPs (e.g. s738333)
                # samples = ["2" if x.startswith("1/1") else x for x in samples]

                sflags = set(flags)
                if len(sflags) == 1 and sflags == set(["NA"]):
                    #print "Skipping", flags
                    continue
    
                valid_count += 1
                vcfmatrix[k-ncomments, :ncols] = np.asarray(flags)

        self.logger.info("%i valid variations found" % valid_count)

        if include_random:
            vcfmatrix[:, ncols] = self.generate_random_sample(nrows)
            samplenames.append("random")

        return vcfmatrix, samplenames

    def read_samplefile(self, sfile):
        sample_paths = []
        sample_labels = []
        with open(sfile, "r") as s:
            for line in s:
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                sample_paths.append(cols[0])
                sample_labels.append(cols[1])

        return sample_paths, sample_labels

    def clusterplot(self, snpmatrix, distmeth="canberra", linkmeth="single",
                    filename="foo", fileformat="png"):
        mat = snpmatrix.transpose()

        dist_matrix = pdist(mat, distmeth)
        linkage_matrix = linkage(dist_matrix, linkmeth)

        fig = plt.figure()
        plt.clf()
        dendrogram(linkage_matrix, labels=self.samplenames,
                   leaf_label_rotation=45)

        fig.savefig(".".join([filename, fileformat]), format=fileformat)


if __name__ == "__main__":
    #import argparse
    import sys

    #parser = argparse.ArgumentParser(description='Verify tumor-normal pair identities')
    #args = parser.parse_args()

    reference = "/runarea1/references/b37d5-dbSNP135/human_g1k_v37_decoy_dbSNP135_10M.fasta"
    regions = "/casa3/project_archive/v2exome-analysis-tools/ucsc-snp135Common-b37-hapmap-maf5allpops.bed"

    try:
        file = sys.argv[1]
        samplefile = sys.argv[2]

        tnv = tnverify(file, samplefile)
    except KeyboardInterrupt:
        print "Program interrupted by user, exiting..."
    except Exception as e:
        import traceback
        print >>sys.stderr, traceback.format_exc()
