process iqtree {
    publishDir "${params.outDir}/results", mode: 'copy'

    container 'thetinhornfish/iqtree2:2.0.7'

    input:
    path alignment_file
    path partition_file

    output:
    path "*.treefile"
    path "*.log"
    path "*.iqtree"
    path "*.nwk"

    script:
    """
    echo "Running IQ-TREE with alignment: ${alignment_file} and partitions: ${partition_file}"
    
    iqtree2 -s "${alignment_file}" \
        -spp "${partition_file}" \
        -m MFP \
        -bb 1000 \
        -nt AUTO \
        -pre iqtree_output
    
    # Copy treefile to .nwk format
    cp iqtree_output.treefile iqtree_output.nwk
    """
}
