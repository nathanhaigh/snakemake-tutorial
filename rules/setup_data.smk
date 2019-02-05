from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule download_chromosome:
	input:
		lambda wildcards: HTTP.remote("urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/iwgsc_refseqv1.0_" + wildcards.chr + ".fsa.zip", keep_local=True),
	output:
		"references/iwgsc_refseqv1.0_{chr}.fsa.zip"
	shell:
		"""
		mv {input} {output}
		"""

rule bgzip_chromosome:
	input:
		"references/iwgsc_refseqv1.0_{chr}.fsa.zip",
	output:
		"references/iwgsc_refseqv1.0_{chr}.fasta.gz",
	conda:
		"../envs/tutorial.yml"
	threads:
		MAX_THREADS
	benchmark:
		repeat("benchmarks/bgzip_chromosome/{chr}.txt", N_BENCHMARKS),
	shell:
		"""
		unzip -p {input} \
		  | bgzip --threads {threads} \
		  > {output}
		"""

rule bgzip_chromosome_subregion:
	input:
		"references/iwgsc_refseqv1.0_{chr}.fasta.gz",
	output:
		"references/{chr}:{start}-{end}.fasta.gz",
	conda:
		"../envs/tutorial.yml"
	threads:
		MAX_THREADS
	benchmark:
		repeat("benchmarks/bgzip_chromosome_subregion/{chr}:{start}-{end}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools faidx {input} {wildcards.chr}:{wildcards.start}-{wildcards.end} \
		  | bgzip --threads {threads} \
		  > {output}
		"""


rule extract_reads:
	output:
		r1 = "raw_reads/{accession}_R1.fastq.gz",
		r2 = "raw_reads/{accession}_R2.fastq.gz",
		index = temp("{accession}.realigned.bam.bai"),
	conda:
		"../envs/tutorial.yml"
	params:
		base_url = "http://crobiad.agwine.adelaide.edu.au/dawn/jbrowse-prod/data/local/by_chr/mapped_reads_merged/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.gz/minimap2_defaults/whole_genome/PE/BPA/",
		chr      = "chr4A_part2",
		start    = "235500000",
		end      = "235558000",
	benchmark:
		repeat("benchmarks/extract_reads/{accession}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools view -hu "{params.base_url}/chr4A_part2/{wildcards.accession}.realigned.bam" {params.chr}:{params.start}-{params.end} \
		  | samtools collate -uO - \
		  | samtools fastq -F 0x900 -1 {output.r1} -2 {output.r2} -s /dev/null -0 /dev/null -
		"""

