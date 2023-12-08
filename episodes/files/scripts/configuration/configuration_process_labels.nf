nextflow.enable.dsl = 2

process P1 {

    label "bigmem"

    script:
    """
    echo P1: Using $task.cpus cpus and $task.memory memory.
    """
}

process P2 {

    label "bigmem"

    script:
    """
    echo P2: Using $task.cpus cpus and $task.memory memory.
    """
}

workflow {
    P1()
    P2()
}
