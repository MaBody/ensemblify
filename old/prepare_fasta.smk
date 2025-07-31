# Remove gaps from already aligned data
rule prepare_sequences:
    input:
        msa = DATA_DIR / "{dataset}"
    output:
        sequences = OUT_DIR / "{dataset}" / "sequences.fasta"
    params:
        dataset = lambda wildcards: wildcards.dataset
    run:
        # Only process if dataset is included in sample
        if params.dataset in datasets:
            msa_obj = Alignment(input.msa).msa

            with open(output.sequences, "w") as outfile:
                for sequence in msa_obj:
                    outfile.write(f">{sequence.id}\n")
                    outfile.write(str(sequence.seq).replace("-", ""))
                    outfile.write("\n")
