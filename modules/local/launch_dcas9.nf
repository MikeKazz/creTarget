process LAUNCH_DCAS9 {
    tag "${chunk_id}"

    input:
    tuple val(chunk_id), path(annotated)
    path launch_script
    path ckpt_dir

    output:
    tuple val(chunk_id), path("${chunk_id}.with_launch.tsv"), emit: annotated
    path "versions.yml", emit: versions

    script:
    def python_arg = params.launch_python ? "--python-bin ${params.launch_python}" : ''

    /*
     * GPU assignment.
     *
     * Under a scheduler with GPU resources requested (--gres=gpu:1 on SLURM,
     * and the equivalent elsewhere), the scheduler sets CUDA_VISIBLE_DEVICES
     * itself and we must not touch it -- overriding it would hand the task a
     * GPU that was allocated to somebody else.
     *
     * With the local executor there is nothing doing that allocation, so every
     * concurrent task would pile onto GPU 0 and run it out of memory. When
     * --gpu_count is set, each task is pinned to one device round-robin by
     * task index instead. Pair it with maxForks = gpu_count so no two tasks
     * share a device.
     */
    def gpu_pin = params.gpu_count
        ? "export CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-${(task.index - 1) % (params.gpu_count as int)}}"
        : "# no --gpu_count set; leaving CUDA_VISIBLE_DEVICES as the scheduler left it"
    """
    ${gpu_pin}
    echo "CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-<unset>}" >&2

    cretarget score-launch-dcas9 \\
        --annotated ${annotated} \\
        --launch-script ${launch_script} \\
        --ckpt-dir ${ckpt_dir} \\
        ${python_arg} \\
        --workdir . \\
        --out ${chunk_id}.with_launch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cretarget: \$(cretarget --version | sed 's/creTarget //')
    END_VERSIONS
    """
}
