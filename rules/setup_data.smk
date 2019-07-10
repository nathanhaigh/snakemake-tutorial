ruleorder:
	extract_chromosome_subregion > get_chromosome

rule get_chromosome:
	output:
		"references/{chr}.fasta.gz",
	conda:
		"../envs/tutorial.yml"
	params:
		reference_url = "http://crobiad.agwine.adelaide.edu.au/dawn/jbrowse-prod/data/wheat_full/references/161010_Chinese_Spring_v1.0_pseudomolecules.fasta.gz",
	threads:
		MAX_THREADS
	benchmark:
		repeat("benchmarks/get_chromosome_subregion/{chr}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools faidx {params.reference_url} {wildcards.chr} \
		  | bgzip --threads {threads} \
		  > {output}
		"""

rule extract_chromosome_subregion:
	input:
		"references/{chr}.fasta.gz",
	output:
		"references/{chr}:{start}-{end}.fasta.gz",
	conda:
		"../envs/tutorial.yml"
	threads:
		MAX_THREADS
	benchmark:
		repeat("benchmarks/extract_chromosome_subregion/{chr}:{start}-{end}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools faidx {input} {wildcards.chr}:{wildcards.start}-{wildcards.end} \
		  | bgzip --threads {threads} \
		  > {output}
		"""

rule extract_reads:
	input:
		reference = "references/{chr}.fasta.gz",
		index     = expand("references/{{chr}}.fasta.gz.{ext}", ext=["fai", "gzi"]),
	output:
		r1 = "raw_reads/{chr}:{start}-{end}/{accession}_R1.fastq.gz",
		r2 = "raw_reads/{chr}:{start}-{end}/{accession}_R2.fastq.gz",
		index = "raw_reads/{chr}:{start}-{end}/{accession}.cram.crai",
	conda:
		"../envs/tutorial.yml"
	params:
		base_url = "http://crobiad.agwine.adelaide.edu.au/dawn/jbrowse-prod/data/wheat_full/minimap2_defaults/whole_genome/PE/BPA",
	benchmark:
		repeat("benchmarks/extract_reads/{chr}:{start}-{end}/{accession}.txt", N_BENCHMARKS),
	shell:
		"""
		# Need to fudge things more than I'd like, since the cari file that samtools downloads for each chromosome will clobber each other
		cd raw_reads/{wildcards.chr}:{wildcards.start}-{wildcards.end}/
		samtools view -hu --reference ../../{input.reference} "{params.base_url}/{wildcards.chr}/{wildcards.accession}.cram" {wildcards.chr}:{wildcards.start}-{wildcards.end} \
		  | samtools collate -uO - \
		  | samtools fastq -F 0x900 -1 {wildcards.accession}_R1.fastq.gz -2 {wildcards.accession}_R2.fastq.gz -s /dev/null -0 /dev/null -
		"""

rule index_fasta:
	input:
		"references/{chr}.fasta.gz",
	output:
		expand("references/{{chr}}.fasta.gz.{ext}", ext=["fai", "gzi"]),
	conda:
		"../envs/tutorial.yml"
	benchmark:
		repeat("benchmarks/index_fasta/{chr}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools faidx {input}
		"""
