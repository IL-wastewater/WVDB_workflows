process CHMOD_OUTPUT {
    tag "chmod_output"
    executor 'local'
    cpus 1

    input:
    val output_dirs
    val stats_done

    output:
    path "chmod_complete.txt"

    script:
    def dirs = output_dirs.collect { "'${it}'" }.join(' ')

    """
    set -euo pipefail

    for d in ${dirs}; do
        if [ -d "\$d" ]; then
            find "\$d" -type f -exec chmod a-w {} +
        else
            echo "[WARN] Output directory does not exist: \$d" >&2
        fi
    done

    echo "chmod a-w complete" > chmod_complete.txt
    """
}
