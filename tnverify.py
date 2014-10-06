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


VCF_COL_CHROM = 0
VCF_COL_POS = 1
VCF_COL_FORMAT = 8
VCF_COL_FIRST_SAMPLE = 9


def vcf2flag(x):
    """Converts a VCF genotype (GT) entry into a state 0/1/2; 0 (homozygous
    reference), 1 (heterozygous) and 2 (homozygous SNP)."""
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


def exec_variant_calling_pipeline(samtools_cmd, bcftools_cmd, outfile, logger):
    """Executes a combination of samtools and bcftools to call variants and
    produce a VCF file."""
    import subprocess
    import threading
    from time import sleep
    from Queue import Queue, Empty

    def logq_enqueue(subp, pipe, q):
        while subp.poll() is None:
            line = pipe.readline()
            q.put(line)

    def logq_drain(t1, t2, q, logger):
        while t1.is_alive() or t2.is_alive() and not q.empty():
            try:
                line = q.get_nowait()
            except Empty:
                sleep(5)
            else:
                if line:
                    logger.info(line.rstrip())
                q.task_done()

    with open(outfile, "w") as vcffile:
        try:
            logger.info("Executing cmdline: %s | %s" % (samtools_cmd,
                                                        bcftools_cmd))
            samtools = subprocess.Popen(samtools_cmd.split(),
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
            bcftools = subprocess.Popen(bcftools_cmd.split(),
                                 stdin=samtools.stdout,
                                 stdout=vcffile,
                                 stderr=subprocess.PIPE)
            # Allow samtools to receive a SIGPIPE if bcftools exits.
            samtools.stdout.close()

            q = Queue()
            t_samtools = threading.Thread(target=logq_enqueue, args=(samtools, samtools.stderr, q))
            t_samtools.daemon = True
            t_bcftools = threading.Thread(target=logq_enqueue, args=(bcftools, bcftools.stderr, q))
            t_bcftools.daemon = True
            t_log = threading.Thread(target=logq_drain, args=(t_samtools, t_bcftools, q, logger))
            t_log.daemon = True

            t_samtools.start()
            t_bcftools.start()
            t_log.start()

            t_samtools.join()
            logger.debug("samtools thread terminated")
            t_bcftools.join()
            logger.debug("bcftools thread terminated")
            t_log.join()
            logger.debug("log thread terminated")

            exitword = "prematurely" if samtools.returncode < 0 else "successfully"
            logger.debug("samtools exited %s (exit code: %s)" % (exitword, samtools.returncode))
            exitword = "prematurely" if bcftools.returncode < 0 else "successfully"
            logger.debug("bcftools exited %s (exit code: %i)" % (exitword, bcftools.returncode))
        except:
            raise


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
                columns = len(c[VCF_COL_FIRST_SAMPLE:])

    f.seek(0)

    return rows, columns, comments


class tnverify:

    def __init__(self, workdir, regions, reference, bcftools_prefix="bcftools_",
                 vcffile=None, samplefile=None, logger=None):
        if logger:
            self.logger = logger
        else:
            self.logger = init_logger()

        self.workdir = workdir
        self.bcftools_prefix = bcftools_prefix
        self.regions = regions
        self.vcffile = vcffile
        self.samplefile = samplefile
        self.reference = reference

        self.logger.info("Specified parameters:")
        if self.workdir is not None:
            self.logger.info("Work directory: %s" % os.path.abspath(self.workdir))
        if self.regions is not None:
            self.logger.info("Regions file: %s" % os.path.abspath(self.regions))
        if self.vcffile is not None:
            self.logger.info("VCF file: %s" % os.path.abspath(self.vcffile))
        if self.samplefile is not None:
            self.logger.info("Sample map file: %s" % os.path.abspath(self.samplefile))
        if self.reference is not None:
            self.logger.info("Reference file: %s" % os.path.abspath(self.reference))

        self.flagmtxall = []
        self.leaflabelsall = []
        self.genomeposall = []

        if self.vcffile is not None:
            flagmtx, vcflabels, genomepos = self.vcf2ndarray(self.vcffile)
            self.flagmtxall.append(flagmtx)
            self.leaflabelsall.extend(vcflabels)
            self.genomeposall.append(genomepos)

        if self.samplefile is not None:
            self.sample_paths, self.sample_labels = self.read_samplefile()

            logger.info("Calling SNPs for %i samples: %s" %
                        (len(self.sample_labels), ", ".join(self.sample_labels)))

            snpcalling_outfiles = self.call_snps()

            # XXX change sample labels in the VCF file

            mtx_mergelist = []
            for i, ofile in enumerate(snpcalling_outfiles):
                flagmtx, vcflabels, genomepos = self.vcf2ndarray(ofile)
                self.flagmtxall.append(flagmtx)
                self.genomeposall.append(genomepos)

            self.leaflabelsall.extend(self.sample_labels)

        # determine the SNP positions common to all samples
        gpos_common = self.get_common_genome_positions(self.genomeposall)

        self.overall_flagmtx = self.merge_snpmtx(gpos_common, self.genomeposall,
                                                 self.flagmtxall)

        self.overall_leaf_labels = self.leaflabelsall

        self.filter_uninformative_snps()
        self.add_random_sample()
        self.clusterplot()


    def call_snps(self):
        """Run samtools and bcftools to call SNPs."""
        vcflist = []
        for sample, label in zip(self.sample_paths, self.sample_labels):
            # -I           do not perform indel calling
            # -g           generate BCF output (genotype likelihoods)
            # -u           generate uncompressed BCF output
            # -D           output per-sample DP in BCF (require -g/-u)
            # -B           disable BAQ computation
            samtools_cmd = "samtools mpileup -IguDB -f %s -l %s %s" % (self.reference,
                                                                       self.regions,
                                                                       sample)
            # -v        output potential variant sites only (force -c)
            # -c        SNP calling (force -e)
            # -g        call genotypes at variant sites (force -c)
            bcftools_cmd = "bcftools view -vcg -"

            samplebase = os.path.basename(sample)
            filestem, ext = os.path.splitext(samplebase)
            outfile = ".".join([filestem, "vcf"])

            self.logger.info("Processing sample: %s" % label)
            try:
                exec_variant_calling_pipeline(samtools_cmd, bcftools_cmd, outfile, logger)
            except Exception as e:
                self.logger.warn("Problem while processing sample: %s" % str(e))
            else:
                vcflist.append(outfile)
                self.logger.info("Sample %s processed successfully" % label)

        return vcflist

    def get_common_genome_positions(self, genomepos_list):
        """Returns a set of genome positions common to all SNP matrixes."""
        gpos_setlist = [set(l_gpos) for l_gpos in genomepos_list]
        gpos_common = set.intersection(*gpos_setlist)
        self.logger.info("%i SNPs common to %i samples" % (len(gpos_common),
                                                      len(genomepos_list)))
        return gpos_common

    def merge_snpmtx(self, gpos_common, genomepos_list, mtx_list):
        """Remove SNPs that are not common between all flag matrixes and merge
        the trimmed matrixes together.  Samples (columns) are append on the
        right side of the matrix in order of the list mtx_list."""
        self.logger.info("Creating merged SNP matrix...")
        for k, (mtx, gpos_list) in enumerate(zip(mtx_list, genomepos_list)):
            self.logger.debug("Matrix %i contains %i SNPs and %i samples" % (k, mtx.shape[0], mtx.shape[1]))
            rm_indexes = [i for i, x in enumerate(gpos_list) if x not in gpos_common]

            if len(rm_indexes) > 0:
                self.logger.info("Matrix %i: deleting %i SNPs" % (k, len(rm_indexes)))
                mtx_list[k] = np.delete(mtx, rm_indexes, 0)
                self.logger.info("Matrix %i has new dimensions: %i SNPs and %i samples" % (k, mtx_list[k].shape[0],
                                                 mtx_list[k].shape[1]))
            else:
                self.logger.info("Matrix %i: no changes performed." % k)

        snpmtx_merge = np.concatenate(mtx_list, 1)
        self.logger.info("Merged matrix dimensions: %i SNPs, %i samples" % snpmtx_merge.shape)
        return snpmtx_merge

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
        invalid_count = 0
        with open(vcffile) as vcf_input:
            nrows, ncols, ncomments = get_file_dims(vcf_input)
            self.logger.info("Reading VCF file %s containing %i samples and %i variations" % (vcffile, ncols, nrows))

            genomepos = []
            vcfmatrix = np.ndarray((nrows, ncols))
            for k, line in enumerate(vcf_input):
                if line.startswith("##"):
                    continue
                cols = line.rstrip().split("\t")
                entries = cols[VCF_COL_FIRST_SAMPLE:]
                fmt = cols[VCF_COL_FORMAT]

                if line.startswith("#"):
                    samplenames = entries
                    continue

                # make sure that GT is always first
                if not fmt[:2] == "GT":
                    #self.logger.debug("Skipping entry with format %s (no genotype found)" % fmt)
                    invalid_count += 1
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
                # Identifier for a SNP position. Example: 5:456334 (chrom:pos_on_chrom).
                genomepos.append(":".join([cols[VCF_COL_CHROM], cols[VCF_COL_POS]]))
                vcfmatrix[k-ncomments, :ncols] = np.asarray(flags)

        self.logger.info("Found %i invalid variants in VCF file" % invalid_count)
        self.logger.info("Found %i valid SNPs in VCF file" % valid_count)

        return vcfmatrix, samplenames, genomepos

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
            self.logger.debug("Reading sample map file %s." % self.samplefile)
            for line in s:
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                sample_paths.append(cols[0].strip())
                if len(cols) == 2:
                    sample_labels.append(cols[1].strip())
                elif len(cols) == 1:
                    sample_labels.append(cols[0].strip())

            self.logger.debug("Read %i samples and %i labels from sample map." %
                              (len(sample_paths), len(sample_labels)))

        return sample_paths, sample_labels

    def clusterplot(self, distmeth="canberra", linkmeth="single",
                    filename="tnverify_hierarchical_clustering", fileformat="png"):
        """Writes a hierarchical clustering plot."""
        self.logger.debug("Preparing for hierarchical clustering.")
        mat = self.overall_flagmtx.transpose()

        dist_matrix = pdist(mat, distmeth)
        linkage_matrix = linkage(dist_matrix, linkmeth)

        self.logger.debug("Creating dendrogram.")

        fig = plt.figure()
        plt.clf()
        dendrogram(linkage_matrix, labels=self.overall_leaf_labels,
                   leaf_rotation=45)

        f_out = ".".join([filename, fileformat])
        outpath = os.path.join(self.workdir, f_out)

        fig.savefig(outpath, format=fileformat)
        self.logger.info("Dendrogram plot saved to file %s." % outpath)


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
    parser.add_argument("-w", "--workdir", help="Directory for intermediary files",
                        default=os.path.join(os.path.expanduser("~"), "tnverify_run"))
    parser.add_argument("-s", "--samplemap", help="Map of BAM file to label",
                        type=is_valid_file, default=None)
    parser.add_argument("-r", "--reference", help="Reference FASTA sequence",
                        type=is_valid_file)
    parser.add_argument("-b", "--bed", help="SNP regions in BED format",
                        type=is_valid_file)
    parser.add_argument("-f", "--vcffile", help="VCF file", type=is_valid_file,
                       default=None)
    parser.add_argument("-v", "--verbosity", help="Increase logging verbosity",
                        action="count", default=0)
    parser.add_argument("--version", action="version", version="%(prog)s 0.1")
    args = parser.parse_args()

    if args.vcffile is None and args.samplemap is None:
        parser.error("Either -s/--samplemap or -f/--vcsfile needs to be specified.")

    if args.samplemap is not None:
        if args.reference is None:
            parser.error("--reference is required when -s/--samplemap is specified.")

        if args.bed is None:
            parser.error("--bed is required when -s/--samplemap is specified.")

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
