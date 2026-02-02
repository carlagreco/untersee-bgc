rule assemble:
	input:
		fastq="{sample}.fastq.gz"
	output:
		"assembly/{sample}/assembly.fasta"
	conda:
		"envs/flye.yaml"
	threads: 32
	shell:
		"flye --meta --threads {threads} --nano-hq {input.fastq} --out-dir assembly/{wildcards.sample}"
		
rule medaka:
	input:
		assembly="assembly/{sample}/assembly.fasta",
		fastq="filt_data/{sample}.fastq.gz"
	output:
		"polish/{sample}/{sample}_polished.fasta"
	threads: 32
	conda: 
		"envs/medaka.yaml"
	shell:
		"medaka_consensus -i {input.fastq} -d {input.assembly} -o polish -t {threads} -m r941_min_sup_g507 && cp polish/consensus.fasta {output}"

rule quast:
	input: "polish/{sample}/{sample}_polished.fasta"
	output: "quast/{sample}_quast.done"
	threads: 16
	conda:
		"envs/quast.yaml"
	shell:
		"""
		quast -o quast/{wildcards.sample} -t {threads} {input}
		touch {output}
		""" 

rule map:
	input:
		fasta="polish/{sample}/{sample}_polished.fasta",
		fastq="filt_data/{sample}.fastq.gz"
	output:
		"mapped/{sample}.sam"
	threads: 32
	conda: "envs/medaka.yaml"
	shell:
		"minimap2 -ax map-ont -t {threads} {input.fasta} {input.fastq} > {output}"

# Define rule to convert SAM to BAM
rule sam_to_bam:
	input:
		sam="mapped/{sample}.sam"
	output:
		bam="mapped/{sample}.sorted.bam",
		bai="mapped/{sample}.sorted.bam.bai"
	threads: 2
	conda: "envs/medaka.yaml"
	shell:		
		"""
		samtools view -bS {input.sam} | samtools sort -o {output.bam} -
		samtools index {output.bam}
		"""

# Define rule to input abundance files into MaxBin2 and MetaBAT2
rule maxbin2:
	input:
		abundance="abundance/{sample}_coverage.txt",
		contig="polish/{sample}/{sample}_polished.fasta"
	output:
		bins="maxbin2/{sample}.summary"
	params: basefile="maxbin2/{sample}"
	threads: 16
	conda: "envs/maxbin2.yaml"
	shell:
		"run_MaxBin.pl -abund {input.abundance} -contig {input.contig} -out {params.basefile} -thread {threads}"

rule format_abundance_file:
	input: "{sample}_polished.fasta.depth.txt"
	output: "abundance/{sample}_coverage.txt"
	shell: 
		"""
		awk 'NR>1 {{print $1"\t"$3}}' {input} > {output}
		"""

rule gene_call:
	input: 
		"polish/{sample}/{sample}_polished.fasta"
	output:
		faa="gene_call/{sample}.faa",
		gff="gene_call/{sample}.gff"
	conda: "envs/pyrodigal.yaml"
	threads: 1
	shell:
		"pyrodigal -i {input} -a {output.faa} -f gff -d gene_call/{wildcards.sample}.fna -o {output.gff} -p meta"
		
rule input_for_antismash:
	input: 
		assembly = "polish/{sample}/{sample}_polished.fasta",
		genes = "gene_call/{sample}.gff"
	output:
		assembly="antismash/{sample}_polished.fasta",
		genes="antismash/{sample}.gff"
	shell: "cp {input.assembly} {output.assembly} && cp {input.genes} {output.genes}"
	
rule antismash:
	input: 
		assembly = "antismash/{sample}_polished.fasta",
		genes = "antismash/{sample}.gff"
	output:
		outdir = directory("antismash/{sample}"),
		htmlfile = "antismash/{sample}/index.html"
	threads: 16
	container: "docker://antismash/standalone:6.1.1"
	shell:
		"antismash --tigrfam --cc-mibig --clusterhmmer --output-dir {output.outdir} --cb-knownclusters --cb-general --cb-general -c {threads} --smcog-trees --genefinding-gff3 {input.genes} {input.assembly}"

rule kraken2:
	input: "polish/{sample}/{sample}_polished.fasta"
	output: outfile="kraken2/{sample}_kraken.out",
		report="kraken2/{sample}_kraken.report"
	threads: 16
	conda: "envs/kraken2.yaml"
	shell:
		"kraken2 --db k2_pluspf_20221209 --threads {threads} --use-names --output {output.outfile} --report {output.report} {input}"


rule final_summary:
	input: kraken = "kraken2/{sample}_kraken.out",
		antismash = "antismash/{sample}/index.html",
		maxbinbins = "maxbin2/{sample}.summary",
		metabatbins = "metabat2/{sample}_metabat.file"
	output: "{sample}.done"
	shell: "touch {output}"