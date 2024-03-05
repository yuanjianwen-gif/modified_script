__author__ = 'gpratt'
#### original code is from https://github.com/YeoLab/makebigwigfiles, and I did some changes to fit
##1.scale factor
##2.output name

from subprocess import call
import subprocess
import os
import pysam

def genome_coverage_bed(in_bam=None, out_bed_graph=None, genome=None, strand=None, scale=1, five_prime=False):
    with open(out_bed_graph, 'w') as out_bed_graph:
        call = "genomeCoverageBed -ibam {} -bg -strand {} -scale {} -du".format(in_bam,
                                                                                             strand,
                                                                                             scale)


        if five_prime:
            call += " -5 "
        else:
            call += " -split "
        subprocess.check_call(call, shell=True, stdout=out_bed_graph)


def get_norm_constant(in_bam):
    samfile = pysam.Samfile(in_bam)
    mapped_reads = float(samfile.mapped) / 1000000
    norm_constant = 1. / mapped_reads
    return norm_constant

def bed_graph_to_big_wig(in_bed_graph, genome, out_big_wig):
    priming_call = "bedGraphToBigWig {} {} {}".format(in_bed_graph, genome, out_big_wig)

    with open(os.devnull, 'w') as fnull:
        subprocess.call(priming_call, shell=True, stdout=fnull)
def check_for_index(bamfile):

    """

    Checks to make sure a BAM file has an index, if the index does not exist it is created

    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file

    """

    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))

    if os.path.exists(bamfile + ".bai"):
        return

    if not bamfile.endswith(".bam"):
        raise NameError("file %s not of correct type" % (bamfile))
    else:
        process = call(["samtools", "index", str(bamfile)])

        if process == -11:
            raise NameError("file %s not of correct type" % (bamfile))      
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Makes Pretty bed Graph Files!")
    parser.add_argument("--bam", help="bam file to make bedgraphs from", required=True)
    parser.add_argument("--genome", help="chromosome sizes because some things need it", required=True)
    parser.add_argument("--five_prime", help="plot only 5' end of fragment", action="store_true", default=False)
#    parser.add_argument("--bw_pos", help="positive bw file name", required=True) ###changed
#    parser.add_argument("--bw_neg", help="negative bw file name", required=True) ###changed
    parser.add_argument("--pre", help="prefix for  bg and bw file name", required=True)
    parser.add_argument("--scale_factor", help="scale factor which multiply original signal, 0 means use mapping 1/(reads/1000000)", default=0, type=float, required=False)
    args = parser.parse_args()
    bamFile = args.bam
    genome = args.genome
    bw_pos=args.pre + ".pos.bw"
    bw_neg=args.pre + ".neg.bw"

    check_for_index(bamFile)
    if  args.scale_factor == 0 :
      norm_constant = get_norm_constant(bamFile)
    else :
      norm_constant = args.scale_factor
    print("scale_factor=", norm_constant)
    neg_norm_constant = norm_constant * -1

    bedGraphFilePos = args.pre + "pos.bg" #bamFile.replace(".bam", ".norm.pos.bg")
    bedGraphFileNeg = args.pre + "neg.bg" #bamFile.replace(".bam", ".norm.neg.bg")
  
    #Tracks are reversed because orientation of TrueSeq kits
    genome_coverage_bed(in_bam=bamFile,
                        out_bed_graph=bedGraphFilePos,
                        strand="-", genome=genome,
                        scale=norm_constant,
                        five_prime=args.five_prime)
    bed_graph_to_big_wig(bedGraphFilePos, genome, bw_pos)

    genome_coverage_bed(in_bam=bamFile,
                        out_bed_graph=bedGraphFileNeg,
                        strand="+", genome=genome,
                        scale=neg_norm_constant,
                        five_prime=args.five_prime)
    bed_graph_to_big_wig(bedGraphFileNeg, genome, bw_neg)
