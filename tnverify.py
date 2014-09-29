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

import logging
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
        raise NotImplementedError("Could not parse VCF entry %s" % x)
    x = x.split(":")[0]
    if x == ".":
        return "0"
    if x == "0/0":
        return ("0")
    # if not ":" in x
    #     print "Could not parse VCF entry", x
    #     raise NotImplementedError("Could not parse VCF entry", x)
    if x.startswith("0/"):
        return ("1")
    if len(x) > 2 and x[0] == x[2]:
        return ("2")
    raise NotImplementedError("Don't know how to deal with this VCF entry: %s" % x)


def init_logger(level=logging.INFO, logfile=None):
    """Initializes a logger.  Logs to the console and to a logfile, if specified."""
    logger = logging.getLogger("")
    logger.setLevel(level)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")

    cons_handler = logging.StreamHandler()
    cons_handler.setLevel(level)
    cons_handler.setFormatter(formatter)
    logger.addHandler(cons_handler)

    if logfile:
        file_handler = logging.FileHandler(logfile)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

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

    def __init__(self, workdir, vcffile, samplefile, reference, logger=None):
        if logger:
            self.logger = logger
        else:
            self.logger = init_logger()

        self.logger.info("Specified parameters:\n" +
                        "Work directory: %s\n" % os.path.abspath(workdir) +
                        "VCF file: %s\n" % os.path.abspath(vcffile) +
                        "Sample map file: %s\n" % os.path.abspath(samplefile) +
                        "Reference file: %s\n" % os.path.abspath(reference))

        self.sample_paths, self.sample_labels = self.read_samplefile(samplefile)

        bcftools_out = "bcftools_temp.vcf"

        #call_snps(samplelist, reference, regions, outfile)
    
        self.flagmatrix, self.vcfsamplenames = self.vcf2ndarray(vcffile,
                                                             add_random_sample=False)
        self.filter_uninformative_snps()

        if vcffile:
            leaf_labels = self.vcfsamplenames
        else:
            leaf_labels = self.sample_labels

        self.clusterplot(self.flagmatrix, leaf_labels)

    def call_snps(self, samples, reference, regions, outfile):
        """Run samtools and bcftools to call SNPs."""
        # samtools mpileup -IguDB -f $ref -l $regions ${BamDir}/${samples} | bcftools view -vcg - > ${OutDir}/result-contralat.txt"
        samtools_cmd = "samtools mpileup -IguDB -f %s -l %s %s | bcftools view -vcg - > %s" % (reference, regions, " ".join(samples), outfile)

        exec_and_log(samtools_cmd, this.logger)

    def generate_random_sample(self, length):
        """Return a sample vector with random variant calls."""
        return [random.randint(0, 2) for x in range(length)]

    def vcf2ndarray(self, vcffile, add_random_sample=True):
        """Converts a VCF file into a NumPy ndarray matrix of values
        0 (homozygous reference), 1 (heterozygous) and 2 (homozygous
        SNP).

        Returns the matrix and a list of corresponding column names.
        """
        vcfmatrix = None
        samplenames = None
        valid_count = 0
        with open(vcffile) as vcf_input:
            nrows, ncols, ncomments = get_file_dims(vcf_input)

            # +1 to leave space for a sample with random variant calls
            if add_random_sample:
                vcfmatrix = np.ndarray((nrows, ncols+1))
            else:
                vcfmatrix = np.ndarray((nrows, ncols))
            for k, line in enumerate(vcf_input):
                if line.startswith("##"):
                    continue
                cols = line.rstrip().split("\t")
                entries = cols[col_sample_start:]
                format = cols[col_format]

                if line.startswith("#"):
                    samplenames = entries
                    continue

                # make sure that GT is always first
                if not format[:2] == "GT":
                    self.logger.debug("Skipping entry with format %s (no genotype found)" % format)
                    continue
    
                # convert vcf entries to simple flags
                flags = []
                try:
                    for i in range(len(entries)):
                        f = vcf2flag(entries[i])
                        flags.append(f)
                except NotImplementedError as e:
                    self.logger.debug(str(e))
                # captures also tri-allelic SNPs (e.g. s738333)
                # samples = ["2" if x.startswith("1/1") else x for x in samples]

                sflags = set(flags)
                if len(sflags) == 1 and sflags == set(["NA"]):
                    #print "Skipping", flags
                    continue
    
                valid_count += 1
                vcfmatrix[k-ncomments, :ncols] = np.asarray(flags)

        self.logger.info("%i valid variations found" % valid_count)

        if add_random_sample:
            self.logger.info("Adding a random sample to the SNP matrix.")
            vcfmatrix[:, ncols] = self.generate_random_sample(nrows)
            samplenames.append("random")

        return vcfmatrix, samplenames

    def filter_uninformative_snps(self):
        """Removes uninformative SNPs from the SNP matrix.  A SNP is
        informative, if it has at least two different states (1/2/3)
        across all samples.
        """
        uninf_rows = [x for x in range(self.flagmatrix.shape[0]) if
                 len(np.unique(self.flagmatrix[x, :])) == 1]
        self.logger.debug("SNP matrix before filtering: %i rows, %i cols" %
                          self.flagmatrix.shape)
        self.flagmatrix = np.delete(self.flagmatrix, uninf_rows, 0)
        self.logger.info("Removed %i uninformative SNPs." % len(uninf_rows))
        self.logger.debug("SNP matrix after filtering: %i rows, %i cols" %
                          self.flagmatrix.shape)

    def read_samplefile(self, sfile):
        """Read and parse the sample map file.
        
        Format: bampath<TAB>label
        """
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

    def clusterplot(self, snpmatrix, leaf_labels, distmeth="canberra", linkmeth="single",
                    filename="tnverify_hierarchical_clustering", fileformat="png"):
        """Writes a hierarchical clustering plot."""
        mat = snpmatrix.transpose()

        dist_matrix = pdist(mat, distmeth)
        linkage_matrix = linkage(dist_matrix, linkmeth)

        fig = plt.figure()
        plt.clf()
        dendrogram(linkage_matrix, labels=leaf_labels,
                   leaf_rotation=45)

        fig.savefig(".".join([filename, fileformat]), format=fileformat)


def verbosity_to_loglevel(verbosity, skip=3):
    """Map verbosity level to loglevel.  Skip determines which
    levels will be skipped."""
    levels = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]
    select = min(verbosity + skip, len(levels) - 1)
    return levels[select]


if __name__ == "__main__":
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(description='Verify tumor-normal pair identities')
    parser.add_argument("workdir", help="Directory for intermediary files")
    parser.add_argument("samplemap", help="Map of BAM file to label")
    parser.add_argument("reference", help="Reference FASTA sequence")
    parser.add_argument("regions", help="SNP regions in BED format")
    parser.add_argument("-m", "--matrix", help="Existing SNP matrix")
    parser.add_argument("-v", "--verbosity", help="Increase logging verbosity",
                        action="count", default=0)
    args = parser.parse_args()

    if not args.matrix:
        args.matrix = "tempfile"

    if not args.workdir:
        home = os.path.expanduser("~")
        args.workdir = os.path.join(home, "tnverify_run")

    logfile = os.path.join(args.workdir, "tnverify_log.txt")
    loglevel = verbosity_to_loglevel(args.verbosity)
    logger = init_logger(loglevel, logfile)
    logger.debug("Setting logging verbosity: %s" % loglevel)

    try:
        tnv = tnverify(args.workdir, args.matrix, args.samplemap,
                       args.reference, logger=logger)
    except KeyboardInterrupt:
        logger.info("Program interrupted by user, exiting.")
    except Exception as e:
        import traceback
        logger.debug(traceback.format_exc())
