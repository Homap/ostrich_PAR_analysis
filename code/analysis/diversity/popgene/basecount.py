def seq_length(sequence):
    """Return heterozygosity from
            any sequence as a string """
    # The first thing is to test the input
    if not isinstance(sequence, str):
            raise Exception("Sequence is not a string")
    homozygote='ATCG'
    SNPs='RYSWKMBDHV'
    homozygote_count = len([base.upper() for base in sequence if base.upper() in homozygote])
    SNPs_count = len([base.upper() for base in sequence if base.upper() in SNPs])
    return(homozygote_count+SNPs_count)