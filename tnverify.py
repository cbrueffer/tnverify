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

    def __init__(self, workdir, regions, reference, bcftools_out="bcftools_temp.vcf",
                 vcffile=None, samplefile=None, logger=None):
        if logger:
            self.logger = logger
        else:
            self.logger = init_logger()

        self.workdir = workdir
        self.bcftools_out = os.path.join(self.workdir, bcftools_out)
        self.regions = regions
        self.vcffile = vcffile
        self.samplefile = samplefile
        self.reference = reference

        self.logger.info("Specified parameters:")
        self.logger.info("Work directory: %s" % os.path.abspath(self.workdir))
        self.logger.info("Regions file: %s" % os.path.abspath(self.regions))
        self.logger.info("VCF file: %s" % os.path.abspath(self.vcffile))
        self.logger.info("Sample map file: %s" % os.path.abspath(self.samplefile))
        self.logger.info("Reference file: %s" % os.path.abspath(self.reference))

        if self.samplefile:
            self.sample_paths, self.sample_labels = self.read_samplefile()
            self.call_snps()

        self.flagmatrix, self.vcfsamplenames = self.vcf2ndarray(self.vcffile)

        self.overall_flagmtx = self.newcalled_vcf
        self.overall_leaf_labels = self.sample_labels

        self.filter_uninformative_snps()
        self.add_random_sample()
        self.clusterplot()

    def call_snps(self):
        """Run samtools and bcftools to call SNPs."""
        # -I           do not perform indel calling
        # -g           generate BCF output (genotype likelihoods)
        # -u           generate uncompressed BCF output
        # -D           output per-sample DP in BCF (require -g/-u)
        # -B           disable BAQ computation
        samtools_cmd = "samtools mpileup -IguDB -f %s -l %s %s" % (self.reference,
                                                                   self.regions,
                                                                   " ".join(self.sample_paths))
        # -v        output potential variant sites only (force -c)
        # -c        SNP calling (force -e)
        # -g        call genotypes at variant sites (force -c)
        bcftools_cmd = "bcftools view -vcg - > %s" % self.bcftools_out

        exec_cmd = "|".join([samtools_cmd, bcftools_cmd])
        exec_and_log(exec_cmd, self.logger)

    def add_random_sample(self):
        """Adds a sample consisting of random variant calls to the flag
        matrix, and the "random" name to the leaf labels."""
        self.logger.info("Adding a random sample to the SNP matrix.")
        length = self.overall_flagmtx.shape[0]  # number of rows
        rsamp = np.random.randint(0, 3, size=(length, 1))
        self.overall_flagmtx = np.append(self.overall_flagmtx, rsamp, axis=1)
        self.overall_leaf_labels.append("random")

    def vcf2ndarray(self, vcffile):
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

        return vcfmatrix, samplenames

    def filter_uninformative_snps(self):
        """Removes uninformative SNPs from the SNP matrix.  A SNP is
        informative, if it has at least two different states (1/2/3)
        across all samples.
        """
        uninf_rows = [x for x in range(self.overall_flagmtx.shape[0]) if
                 len(np.unique(self.overall_flagmtx[x, :])) == 1]
        self.logger.debug("SNP matrix before filtering: %i rows, %i cols" %
                          self.overall_flagmtx.shape)
        self.overall_flagmtx = np.delete(self.overall_flagmtx, uninf_rows, 0)
        self.logger.info("Removed %i uninformative SNPs." % len(uninf_rows))
        self.logger.debug("SNP matrix after filtering: %i rows, %i cols" %
                          self.overall_flagmtx.shape)

    def read_samplefile(self):
        """Read and parse the sample map file.

        Format: bampath<TAB>label
        """
        sample_paths = []
        sample_labels = []
        with open(self.samplefile, "r") as s:
            for line in s:
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                sample_paths.append(cols[0])
                sample_labels.append(cols[1])

        return sample_paths, sample_labels

    def clusterplot(self, distmeth="canberra", linkmeth="single",
                    filename="tnverify_hierarchical_clustering", fileformat="png"):
        """Writes a hierarchical clustering plot."""
        mat = self.overall_flagmtx.transpose()

        dist_matrix = pdist(mat, distmeth)
        linkage_matrix = linkage(dist_matrix, linkmeth)

        fig = plt.figure()
        plt.clf()
        dendrogram(linkage_matrix, labels=self.overall_leaf_labels,
                   leaf_rotation=45)

        fig.savefig(".".join([filename, fileformat]), format=fileformat)


def verbosity_to_loglevel(verbosity, skip=3):
    """Map verbosity level to loglevel.  Skip determines which
    levels will be skipped."""
    levels = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]
    select = min(verbosity + skip, len(levels) - 1)
    return levels[select]


def is_valid_file(path):
    """Tests whether or not a file can be opened read-only.

    Either returns the file name, or raises a parser exception."""
    try:
        testfile = open(path, "r")
        testfile.close()
        return path
    except:
        msg = "Cannot open file %s!" % path
        raise argparse.ArgumentTypeError(msg)


if __name__ == "__main__":
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(description='Verify tumor-normal pair identities')
    parser.add_argument("workdir", help="Directory for intermediary files")
    parser.add_argument("samplemap", help="Map of BAM file to label", type=is_valid_file)
    parser.add_argument("reference", help="Reference FASTA sequence", type=is_valid_file)
    parser.add_argument("regions", help="SNP regions in BED format", type=is_valid_file)
    parser.add_argument("-f", "--vcffile", help="VCF file", type=is_valid_file)
    parser.add_argument("-v", "--verbosity", help="Increase logging verbosity",
                        action="count", default=0)
    parser.add_argument("--version", action="version", version="%(prog)s 0.1")
    args = parser.parse_args()

    if not args.vcffile:
        args.vcffile = None

    if not args.workdir:
        home = os.path.expanduser("~")
        args.workdir = os.path.join(home, "tnverify_run")

    logfile = os.path.join(args.workdir, "tnverify_run.log")
    loglevel = verbosity_to_loglevel(args.verbosity)
    logger = init_logger(loglevel, logfile)
    logger.debug("Setting logging verbosity: %s" % loglevel)

    try:
        tnv = tnverify(args.workdir, args.regions, args.reference, vcffile=args.vcffile,
                       samplefile=args.samplemap, logger=logger)
    except KeyboardInterrupt:
        logger.info("Program interrupted by user, exiting.")
    except Exception as e:
        import traceback
        logger.debug(traceback.format_exc())
