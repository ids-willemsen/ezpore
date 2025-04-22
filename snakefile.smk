#from lib2to3.pgen2.tokenize import group

#from poetry.utils.helpers import directory

configfile: "settingsfile.yaml"

def validate_config(config):
    required_keys = [
        'demultiplex', 'min', 'max', 'quality', 'trim_primers', 
        'primer_error_rate', 'min_abundance', 'cluster_perc', 'clustering',
        'rank', 'threads', 'input_file', 'group', 'barcode_file'
    ]

    # Check for missing keys
    for key in required_keys:
        if key not in config:
            raise KeyError(f"Missing required config parameter: {key}")

    # Check types
    if not isinstance(config["demultiplex"], bool):
        raise TypeError("Parameter 'demultiplex' should be a boolean.")
    if not isinstance(config["clustering"], bool):
        raise TypeError("Parameter 'clustering' should be a boolean.")
    if not isinstance(config["min"], int):
        raise TypeError("Parameter 'min' should be an integer.")
    if not isinstance(config["max"], int):
        raise TypeError("Parameter 'max' should be an integer.")
    if not isinstance(config["quality"], int):
        raise TypeError("Parameter 'quality' should be an integer.")
    if not isinstance(config["trim_primers"], bool):
        raise TypeError("Parameter 'trim_primers' should be a boolean.")
    if not isinstance(config["primer_error_rate"], float):
        raise TypeError("Parameter 'primer_error_rate' should be a float.")
    if not isinstance(config["min_abundance"], float):
        raise TypeError("Parameter 'min_abundance' should be a float.")
    if not isinstance(config["threads"], int):
        raise TypeError("Parameter 'threads' should be an integer.")
    if not isinstance(config["input_file"], str):
        raise TypeError("Parameter 'input_file' should be a string.")
    if not isinstance(config["group"], str):
        raise TypeError("Parameter 'group' should be a string.")
    if not isinstance(config["barcode_file"], str):
        raise TypeError("Parameter 'barcode_file' should be a string.")

    # Check if 'cluster_perc' is either a boolean or a float
    # If cluster_perc is a float, check its range
    if isinstance(config["cluster_perc"], float):
        if not (0 <= config["cluster_perc"] <= 1):
            raise ValueError("Parameter 'cluster_perc' should be between 0 and 1 when it's a float.")

    # Check value ranges
    if config["min"] < 0:
        raise ValueError("Parameter 'min' should be non-negative.")
    if config["max"] <= config["min"]:
        raise ValueError("Parameter 'max' should be greater than 'min'.")
    if not (0 <= config["primer_error_rate"] <= 1):
        raise ValueError("Parameter 'primer_error_rate' should be between 0 and 1.")

    # Check if 'group' is one of the allowed values
    allowed_groups = {"16S_bac", "18S_nem", "ITS_fungi"}
    if config["group"] not in allowed_groups:
        raise ValueError(f"Parameter 'group' must be one of {allowed_groups}, mind that this is case sensitive. "
                         f"Found: {config['group']}")

# Call the validation function
print(config)
validate_config(config)

if not config.get("clustering", True):  # Defaults to True if not specified
    print("Clustering is set to False")

def read_barcodes(barcode_file): #function to parse barcodes from barcode file
    with open(barcode_file) as f:
        return [line.strip() for line in f if line.strip()]  # Remove empty lines

barcode_file = config["barcode_file"]
# Load barcodes
barcodes = read_barcodes(barcode_file)

#only add rules when emu/vsearch is selected
rule_all_classifier = list()

if config["classifier"] == "emu":
    rule_all_classifier.append("results/emu-combined-{}-counts.tsv".format(config["rank"]))
    #rule_all_classifier.append(expand("results/{barcode}_rel-abundance.tsv",barcode=barcodes))
#if config["classifier"] == vsearch:
#   etc

rule all:
    input:
        #"{}".format(config["group"]),
        #expand("selected_demux_barcodes/{barcode}.fastq",barcode=barcodes),
        #expand("filtered/filtered_{barcode}.fastq", barcode=barcodes),
        #expand("classifier_input/{barcode}.fasta",barcode=barcodes),
        rule_all_classifier

#select correct database url
url = ""

if config["group"] == "18S_nem":
    url = "https://www.dropbox.com/scl/fi/r5llvu4s6x32pcr3znsxc/18S_nem.zip?rlkey=1f8g7qbsyvk9oo447ha5l6m97&st\
    =ujkcvtae&dl=1"
elif config["group"] == "16S_bac":
    url = "https://www.dropbox.com/scl/fi/zrh0ujr4bkj4uqbq2c7bo/16S_bac.zip?rlkey=mwr6njaxjvd6yes9mmpxg779m&st\
    =focsovmh&dl=1"
elif config["group"] == "ITS_fungi":
    url = "https://www.dropbox.com/scl/fi/2mjmt34z0zmr204b3ed2s/ITS_fun.zip?rlkey=ew99c0jzq9ujmlcf1b67vemch&st\
    =tehvlq24&dl=1"



rule download_database:
    output:
        directory("{}".format(config["group"]))
    conda:
        "ezpore_conda.yaml"
    params:
        url = f"{url}",
        group = config["group"]
    log: "logs/download_database.log"
    shell:
        """
        wget -O {params.group}.zip "{params.url}"
        unzip {params.group}.zip
        """

if config.get("demultiplex", None) is True:  # Only run if demultiplexing is enabled
    rule demultiplex:
        input:
            fastq=config["input_file"]
        output:
            fastq_files = expand("demux/unknown_run_id_EXP-NBD196_{barcode}.fastq",barcode=barcodes),
            temp1 = temp("dorado-0.8.3-linux-x64.tar.gz"),
            temp2 = directory(temp("dorado-0.8.3-linux-x64/"))
        params:
            demux_dir = "demux"
        threads: config["threads"]
        log: "logs/demultiplex.log"
        shell:
            """
            wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.3-linux-x64.tar.gz
            tar -xvzf dorado-0.8.3-linux-x64.tar.gz
            ./dorado-0.8.3-linux-x64/bin/dorado demux \
                --kit-name EXP-NBD196 \
                --emit-fastq \
                --output-dir {params.demux_dir} \
                {input.fastq}
            """

    rule select_files: #copies only the selected files based on barcode_file.txt
        input:
            "demux/unknown_run_id_EXP-NBD196_{barcode}.fastq"
        output:
            "selected_demux_barcodes/{barcode}.fastq"
        run:
            import shutil
            src = input[0]
            dst = output[0]
            shutil.copy(src, dst)  # Copy file to destination

else:
    rule select_files_nodemux: #copies only the selected files based on barcode_file.txt
        input:
            "demux/{barcode}.fastq"
        output:
            "selected_demux_barcodes/{barcode}.fastq"
        run:
            import shutil
            src = input[0]
            dst = output[0]
            shutil.copy(src, dst)  # Copy file to destination
            os.remove(src)         # Delete the original file


rule quality_filtering: #runs nanofilt to filter for quality and length
    input:
        fastq="selected_demux_barcodes/{barcode}.fastq"
    output:
        "filtered/filtered_{barcode}.fastq"
    conda:
        "ezpore_conda.yaml"
    log: "logs/nanofilt.log"
    params:
        min_length = config["min"],
        max_length = config["max"],
        quality = config["quality"]
    shell:
        """
        NanoFilt --length {params.min_length} --maxlength {params.max_length} -q {params.quality} {input.fastq} > {output}
        """


if config["group"] == "16S_bac" or config["group"] == "18S_nem":
    if config.get("trim_primers", None) is False and config.get("clustering", None) is False:
        rule move_emu_filter:  #temporarily moves files as input for vsearch
            input:
                "filtered/filtered_{barcode}.fastq"
            output:
                "classifier_input/filtered_{barcode}.fasta"
            conda:
                "ezpore_conda.yaml"
            shell:
                """
                cp {input} {output}
                """

    if config.get("trim_primers", None) is True:
        rule primer_trimming: #trimps primers with cutadapt
            input:
                fastq="filtered/filtered_{barcode}.fastq"
            output:
                "trimmed/trimmed_{barcode}.fastq"
            conda:
                "ezpore_conda.yaml"
            log: "logs/cutadapt.log"
            params:
                primer_error_rate = config["primer_error_rate"],
                threads = config["threads"],
                fw_primer = config["forward_primer"],
                rv_primer = config["reverse_primer"]
            shell:
                """
                cutadapt --error-rate {params.primer_error_rate} --match-read-wildcards --revcomp --cores {params.threads} \
                -g {params.fw_primer} -a {params.rv_primer} --times 2 --quiet {input.fastq} > {output}
                """

        if config.get("clustering", None) is True:
            rule move_vsearch: #temporarily moves files as input for vsearch
                input:
                    "trimmed/trimmed_{barcode}.fastq"
                output:
                    temp("vsearch_input/{barcode}.fastq")
                conda:
                    "ezpore_conda.yaml"
                shell:
                    """
                    cp {input} {output}
                    """
        else:
            rule move_emu_trim:  #temporarily moves files as input for vsearch
                input:
                    "trimmed/trimmed_{barcode}.fastq"
                output:
                    "classifier_input/{barcode}.fasta"
                conda:
                    "ezpore_conda.yaml"
                shell:
                    """
                    cp {input} {output}
                    """


if config["group"] == "ITS_fungi": #extracts ITS for fungi
    rule ITSexpress:
        input:
            fastq = "filtered/filtered_{barcode}.fastq"
        output:
            "ITS_extract/ITS_{barcode}.fastq"
        params:
            threads = config["threads"]
        log: "logs/itsexpress.log"
        shell:
            "itsxpress --fastq {input.fastq} --single_end --outfile {output} --region ALL --threads {params.threads}"

    if config.get("clustering",None) is True:
        rule move_vsearch: #temporarily moves files as input for vsearch
            input:
                "ITS_extract/ITS_{barcode}.fastq"
            output:
                temp("vsearch_input/{barcode}.fastq")
            shell:
                """
                cp {input} {output}
                """
    else:
        rule move_emu_ITS:  #temporarily moves files as input for vsearch
            input:
                "ITS_extract/ITS_{barcode}.fastq"
            output:
                "classifier_input/{barcode}.fasta"
            conda:
                "ezpore_conda.yaml"
            shell:
                """
                cp {input} {output}
                """


if config.get("clustering", None) is True: #!= FALSE as cluster_perc can range between 0-1,
    rule vsearch_clustering: #clustering with vsearch
        input:
            "vsearch_input/{barcode}.fastq"
        output:
            "clustered/clustered_{barcode}.fasta"
        params:
            cluster_perc = config["cluster_perc"],
            threads = config["threads"]
        conda:
            "ezpore_conda.yaml"
        log: "logs/vsearch_cluster.log"
        shell:
            """
            vsearch --cluster_fast {input} --id {params.cluster_perc} --sizeout --threads {params.threads} \
                --consout {output}
            """

    rule vsearch_rereplicate: #rereplicates with vsearch
        input:
            "clustered/clustered_{barcode}.fasta"
        output:
            "clustered_rerep/rerep_clustered_{barcode}.fasta"
        conda:
            "ezpore_conda.yaml"
        log: "logs/vsearch_rereplicate.log"
        shell:
            """
            vsearch --rereplicate {input} --relabel sequence --output {output}
            """

    rule move_emu_cluster:  #temporarily moves files as input for vsearch
        input:
            "clustered_rerep/rerep_clustered_{barcode}.fasta"
        output:
            "classifier_input/{barcode}.fasta"
        conda:
            "ezpore_conda.yaml"
        shell:
            """
            cp {input} {output}
            """


#options: 1. vsearch - ezpz; 2. no vsearch, ITS extract; 3. no vsearch, trim primers; 4. no vsearch, no primer trimming


if config["classifier"] == "emu":
    rule emu: #runs emu
        input:
            fasta = "classifier_input/{barcode}.fasta",
            db_path = config["group"]
        output:
            "results/{barcode}_rel-abundance.tsv",
        params:
            threads = config["threads"],
            min_abundance = config["min_abundance"]
        conda:
            "ezpore_conda.yaml"
        log: "logs/emu.log"
        shell:
            """
            emu abundance {input.fasta} --db {input.db_path} --keep-counts --min-abundance {params.min_abundance} \
                --type map-ont --threads {params.threads}
            """

    rule emu_combine: #combines results into the final OTU table
        input:
            expand("results/{barcode}_rel-abundance.tsv",barcode=barcodes),
        output:
            "results/emu-combined-{}-counts.tsv".format(config["rank"]),
            "results/emu-combined-{}.tsv".format(config["rank"])
        params:
            rank = config["rank"]
        conda:
            "ezpore_conda.yaml"
        shell:
            """
            emu combine-outputs results {params.rank} --counts
            emu combine-outputs results {params.rank}
            """

