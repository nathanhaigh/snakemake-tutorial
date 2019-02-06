ACCESSIONS = [
  "ACBarrie",
#  "Alsen",
#  "Baxter",
#  "Chara",
#  "Drysdale",
#  "Excalibur",
#  "Gladius",
#  "H45",
#  "Kukri",
#  "Pastor",
#  "RAC875",
#  "Volcanii",
#  "Westonia",
#  "Wyalkatchem",
#  "Xiaoyan",
#  "Yitpi",
]

MAX_THREADS = 32
ADAPTERS = "TruSeq3-PE.fa"
CHR = "chr4A"
CHR_START = "688055092"
CHR_END = "688113092"
REFERENCE = "references/" + CHR + ":" + CHR_START + "-" + CHR_END + ".fasta.gz"

N_BENCHMARKS = 3


from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

singularity:
	"docker://continuumio/miniconda3:4.5.12"

############################
# Include other Snakefiles #
############################
include:
	"rules/setup_data.smk"

#######################################
# Convienient rules to define targets #
#######################################
localrules:
	all

rule all:
	input:
		"reports/raw_reads_multiqc.html",
		"reports/qc_reads_multiqc.html",
		expand("mapped/{accession}.bam", accession=ACCESSIONS),


rule setup_data:
	input:
		REFERENCE,
		expand("raw_reads/{accession}_R{read}.fastq.gz", accession=ACCESSIONS, read=[1,2]),


################
# Rules Proper #
################

rule fastqc_raw:
	input:
		"raw_reads/{prefix}.fastq.gz",
	output:
		zip  = "reports/raw_reads/{prefix}_fastqc.zip",
		html = "reports/raw_reads/{prefix}_fastqc.html",
	conda:
		"envs/tutorial.yml"
	threads:
		MAX_THREADS
	benchmark:
		repeat("benchmarks/fastqc_raw/{prefix}.txt", N_BENCHMARKS),
	wrapper:
		"0.31.1/bio/fastqc"
#	shell:
#		"""
#		fastqc --threads {threads} {input}
#		"""

rule multiqc_raw:
	input:
		expand("reports/raw_reads/{accession}_R{read}_fastqc.zip", accession=ACCESSIONS, read=[1,2]),
	output:
		"reports/raw_reads_multiqc.html",
	conda:
		"envs/tutorial.yml"
	benchmark:
		repeat("benchmarks/multiqc_raw/benchmark.txt", N_BENCHMARKS),
	wrapper:
		"0.31.1/bio/multiqc"
#	shell:
#		"""
#		multiqc --filename {output} {input}
#		"""

rule download_trimmomatic_pe_adapters:
	input:
		lambda wildcards: HTTP.remote("raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/" + wildcards.adapters, keep_local=True),
	output:
		"misc/trimmomatic_adapaters/{adapters}"
	conda:
		"envs/tutorial.yml"
	shell:
		"""
		mv {input} {output}
		"""

rule trimmomatic_pe:
	input:
		r1          = "raw_reads/{accession}_R1.fastq.gz",
		r2          = "raw_reads/{accession}_R2.fastq.gz",
		adapaters   = lambda wildcards: "misc/trimmomatic_adapaters/" + ADAPTERS
	output:
		r1          = "qc_reads/{accession}_R1.fastq.gz",
		r2          = "qc_reads/{accession}_R2.fastq.gz",
		# reads where trimming entirely removed the mate
		r1_unpaired = "qc_reads/{accession}_R1.unpaired.fastq.gz",
		r2_unpaired = "qc_reads/{accession}_R2.unpaired.fastq.gz",
	conda:
		"envs/tutorial.yml"
	params:
		trimmer = [
			"ILLUMINACLIP:misc/trimmomatic_adapaters/" + ADAPTERS + ":2:30:10:3:true",
			"LEADING:2",
			"TRAILING:2",
			"SLIDINGWINDOW:4:15",
			"MINLEN:36",
		],
	benchmark:
		repeat("benchmarks/trimmomatic_pe/{accession}.txt", N_BENCHMARKS),
	wrapper:
		"0.31.1/bio/trimmomatic/pe"

rule fastqc_trimmed:
	input:
		"qc_reads/{prefix}.fastq.gz",
	output:
		zip  = "reports/qc_reads/{prefix}_fastqc.zip",
		html = "reports/qc_reads/{prefix}_fastqc.html",
	conda:
		"envs/tutorial.yml"
	benchmark:
		repeat("benchmarks/fastqc_trimmed/{prefix}.txt", N_BENCHMARKS),
	wrapper:
		"0.31.1/bio/fastqc"
#	shell:
#		"""
#		fastqc {input}
#		"""

rule multiqc_trimmed:
	input:
		expand("reports/qc_reads/{accession}_R{read}_fastqc.zip", accession=ACCESSIONS, read=[1,2]),
	output:
		"reports/qc_reads_multiqc.html",
	conda:
		"envs/tutorial.yml"
	benchmark:
		repeat("benchmarks/multiqc_trimmed/benchmark.txt", N_BENCHMARKS),
	wrapper:
		"0.31.1/bio/multiqc"
#	shell:
#		"""
#		multiqc --filename {output} {input}
#		"""

rule bwa_index:
	input:
		"{ref}"
	output:
		"{ref}.amb",
		"{ref}.ann",
		"{ref}.bwt",
		"{ref}.pac",
		"{ref}.sa"
	conda:
		"envs/tutorial.yml"
#	log:
#		"logs/bwa_index/{ref}.log"
	params:
		prefix    = "{ref}",
		algorithm = "bwtsw"
	benchmark:
		repeat("benchmarks/bwa_index/{ref}.txt", N_BENCHMARKS),
	wrapper:
		"0.31.1/bio/bwa/index"
#	shell:
#		"""
#		bwa index \
#		  -p {params.prefix} \
#		  -a {params.algorithm} \
#		  {input} \
#		2> {log}
#		"""

rule bwa_mem:
	input:
		reference = expand(REFERENCE + ".{ext}", ext=["amb","ann","bwt","pac","sa"]),
		reads     = [ "qc_reads/{sample}_R1.fastq.gz", "qc_reads/{sample}_R2.fastq.gz"],
#		r1        = "qc_reads/{sample}_R1.fastq.gz",
#		r2        = "qc_reads/{sample}_R2.fastq.gz",
	output:
		"mapped/{sample}.bam"
	conda:
		"envs/tutorial.yml"
	params:
		index      = REFERENCE,
		extra      = r"-R '@RG\tID:{sample}\tSM:{sample}'",
		sort       = "none",
		sort_order = "queryname",
		sort_extra = ""
	threads:
		MAX_THREADS
	benchmark:
		repeat("benchmarks/bwa_mem/{sample}.txt", N_BENCHMARKS),
	wrapper:
		"0.31.1/bio/bwa/mem"


