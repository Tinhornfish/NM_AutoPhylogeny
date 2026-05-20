process tree_image_creation {
    publishDir "${params.outDir}/results", mode: 'copy'

    // Use the ete3 image built for this project
    container 'thetinhornfish/ete3:3.1.3'

    input:
    path treefile

    output:
    path "tree.png"

    script:
    """
    echo "Creating tree image from ${treefile}"
    xvfb-run -a ete3 view -t ${treefile} --image tree.png
    """
}
