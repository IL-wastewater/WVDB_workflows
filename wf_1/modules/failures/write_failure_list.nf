process WRITE_FAILURE_LIST {
    tag "failure_list"
    executor 'local'
    cpus 1
    time '10m'

    input:
    path failure_files
    val batch_output_dir

    output:
    path "failure_list.txt"

    script:
    """
    set -euo pipefail

    mkdir -p ${batch_output_dir}

    echo -e "sample\tprocess\treason\texit_status\tfailed_file" > failure_list.txt

    for f in ${failure_files}; do
        sample=\$(awk -F': ' '/^sample:/ {print \$2; exit}' "\$f")
        process=\$(awk -F': ' '/^process:/ {print \$2; exit}' "\$f")
        reason=\$(awk -F': ' '/^reason:/ {print \$2; exit}' "\$f")
        exit_status=\$(awk -F': ' '/^exit_status:/ {print \$2; exit}' "\$f")

        printf "%s\t%s\t%s\t%s\t%s\n" \\
            "\${sample:-NA}" \\
            "\${process:-NA}" \\
            "\${reason:-NA}" \\
            "\${exit_status:-NA}" \\
            "\$f" >> failure_list.txt
    done

    cp failure_list.txt ${batch_output_dir}/failure_list.txt
    """
}
